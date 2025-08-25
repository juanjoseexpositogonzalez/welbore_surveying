"""
Survey Data Trajectory Calculator
A Streamlit application for calculating wellbore trajectory parameters
using the minimum curvature method with DDI calculations.
"""

from __future__ import annotations

import math
from io import StringIO
from typing import Any, Dict, List, Optional

import plotly.graph_objects as go
import polars as pl
import streamlit as st
from plotly.graph_objects import Figure

# Configure Streamlit page
st.set_page_config(
    page_title="Survey Data Trajectory Calculator",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded",
)


def degrees_to_radians(degrees: float) -> float:
    """
    Convert degrees to radians.

    Args:
        degrees: Angle in degrees

    Returns:
        Angle in radians
    """
    return degrees * math.pi / 180


def interpolate_survey_point(
    df_results: pl.DataFrame, target_md: float
) -> Optional[Dict[str, float]]:
    """
    Interpolate survey point at target MD using minimum curvature method.

    Args:
        df_results: DataFrame with calculated trajectory parameters
        target_md: Target measured depth for interpolation

    Returns:
        Dictionary with interpolated survey parameters or None if target_md is out of range
    """
    md_values = df_results["MD"].to_list()

    # Check if target MD is within range
    if target_md < min(md_values) or target_md > max(md_values):
        return None

    # If target MD exactly matches an existing survey point, return that point
    for i, md in enumerate(md_values):
        if abs(md - target_md) < 0.1:  # Within 0.1 ft tolerance
            return {
                "MD": float(df_results["MD"][i]),
                "Inc": float(df_results["Inc"][i]),
                "Azi": float(df_results["Azi"][i]),
                "TVD": float(df_results["TVD"][i]),
                "North": float(df_results["North"][i]),
                "East": float(df_results["East"][i]),
                "DLS": float(df_results["DLS"][i]),
                "DDI": float(df_results["DDI"][i]),
            }

    # Find bracketing survey stations
    upper_idx = None
    for i, md in enumerate(md_values):
        if md > target_md:
            upper_idx = i
            break

    if upper_idx is None or upper_idx == 0:
        return None

    lower_idx = upper_idx - 1

    # Get bracketing survey data
    md1, md2 = md_values[lower_idx], md_values[upper_idx]
    inc1, inc2 = (
        float(df_results["Inc"][lower_idx]),
        float(df_results["Inc"][upper_idx]),
    )
    azi1, azi2 = (
        float(df_results["Azi"][lower_idx]),
        float(df_results["Azi"][upper_idx]),
    )

    # Previous cumulative values
    tvd1 = float(df_results["TVD"][lower_idx])
    north1 = float(df_results["North"][lower_idx])
    east1 = float(df_results["East"][lower_idx])

    # Calculate interpolation factor
    total_course = md2 - md1
    target_course = target_md - md1
    factor = target_course / total_course

    # Linear interpolation for inclination and azimuth
    target_inc = inc1 + (inc2 - inc1) * factor
    target_azi = azi1 + (azi2 - azi1) * factor

    # Handle azimuth wrap-around (e.g., 350° to 10°)
    azi_diff = azi2 - azi1
    if azi_diff > 180:
        azi_diff -= 360
    elif azi_diff < -180:
        azi_diff += 360
    target_azi = azi1 + azi_diff * factor
    if target_azi < 0:
        target_azi += 360
    elif target_azi >= 360:
        target_azi -= 360

    # Apply minimum curvature method for position calculation
    # Convert to radians
    inc1_rad = degrees_to_radians(inc1)
    target_inc_rad = degrees_to_radians(target_inc)
    azi1_rad = degrees_to_radians(azi1)
    target_azi_rad = degrees_to_radians(target_azi)

    # Calculate dogleg angle
    cos_beta = math.cos(target_inc_rad - inc1_rad) - math.sin(inc1_rad) * math.sin(
        target_inc_rad
    ) * (1 - math.cos(target_azi_rad - azi1_rad))

    cos_beta = max(-1.0, min(1.0, cos_beta))
    beta = math.acos(cos_beta)

    # Calculate ratio factor
    if abs(beta) < 1e-10:
        rf = 1.0
    else:
        rf = (2 / beta) * math.tan(beta / 2)

    # Calculate incremental displacements
    delta_e = (
        (target_course / 2)
        * (
            math.sin(inc1_rad) * math.sin(azi1_rad)
            + math.sin(target_inc_rad) * math.sin(target_azi_rad)
        )
        * rf
    )

    delta_n = (
        (target_course / 2)
        * (
            math.sin(inc1_rad) * math.cos(azi1_rad)
            + math.sin(target_inc_rad) * math.cos(target_azi_rad)
        )
        * rf
    )

    delta_v = (target_course / 2) * (math.cos(inc1_rad) + math.cos(target_inc_rad)) * rf

    # Calculate interpolated position
    target_tvd = tvd1 + delta_v
    target_north = north1 + delta_n
    target_east = east1 + delta_e

    # Calculate DLS for this segment
    target_dls = (
        (beta * 180 / math.pi) / target_course * 100 if target_course > 0 else 0.0
    )

    # Simplified DDI calculation (approximate)
    vsec = math.sqrt(target_east**2 + target_north**2)
    if target_tvd > 0:
        target_ddi = math.log10((vsec * target_md) / target_tvd) if vsec > 0 else 0.0
    else:
        target_ddi = 0.0

    return {
        "MD": target_md,
        "Inc": target_inc,
        "Azi": target_azi,
        "TVD": target_tvd,
        "North": target_north,
        "East": target_east,
        "DLS": target_dls,
        "DDI": target_ddi,
    }


def calculate_trajectory_parameters(df: pl.DataFrame) -> pl.DataFrame:
    """
    Calculate trajectory parameters using minimum curvature method.
    Following the DDI calculation steps from the provided methodology.

    Args:
        df: Input DataFrame with columns MD, Inc, Azi

    Returns:
        DataFrame with calculated trajectory parameters

    Raises:
        ValueError: If required columns are missing
    """
    # Validate required columns
    required_columns = {"MD", "Inc", "Azi"}
    if not required_columns.issubset(set(df.columns)):
        missing_cols = required_columns - set(df.columns)
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Initialize result storage
    results: List[Dict[str, float]] = []

    # First row - surface location (assumed values)
    first_row: Dict[str, float] = {
        "MD": float(df["MD"][0]),
        "Inc": float(df["Inc"][0]),
        "Azi": float(df["Azi"][0]),
        "TVD": float(df["MD"][0]) * math.cos(degrees_to_radians(float(df["Inc"][0]))),
        "North": 0.0,
        "East": 0.0,
        "DLS": 0.0,
        "Build_Rate": 0.0,
        "Turn_Rate": 0.0,
        "Vsec": 0.0,
        "Closure_Distance": 0.0,
        "Closure_Azimuth": 0.0,
        "DDI": 0.0,
    }
    results.append(first_row)

    # Calculate for subsequent survey points
    for i in range(1, len(df)):
        # Current and previous survey data
        md1, inc1, azi1 = (
            float(df["MD"][i - 1]),
            float(df["Inc"][i - 1]),
            float(df["Azi"][i - 1]),
        )
        md2, inc2, azi2 = float(df["MD"][i]), float(df["Inc"][i]), float(df["Azi"][i])

        # Step 2.1: Convert to radians
        inc1_rad: float = degrees_to_radians(inc1)
        inc2_rad: float = degrees_to_radians(inc2)
        azi1_rad: float = degrees_to_radians(azi1)
        azi2_rad: float = degrees_to_radians(azi2)

        # Step 2.2: Calculate β (dogleg angle)
        cos_beta: float = math.cos(inc2_rad - inc1_rad) - math.sin(inc1_rad) * math.sin(
            inc2_rad
        ) * (1 - math.cos(azi2_rad - azi1_rad))

        # Ensure cos_beta is within valid range for acos
        cos_beta = max(-1.0, min(1.0, cos_beta))
        beta: float = math.acos(cos_beta)

        # Step 2.3: Calculate RF (Ratio Factor)
        if abs(beta) < 1e-10:  # Essentially zero
            rf: float = 1.0
        else:
            rf = (2 / beta) * math.tan(beta / 2)

        # Course length
        delta_md: float = md2 - md1

        # Step 2.4: Calculate DLS
        dls: float = (beta * 180 / math.pi) / delta_md * 100 if delta_md > 0 else 0.0

        # Step 2.5: Calculate ΔE, ΔN, ΔV
        delta_e: float = (
            (delta_md / 2)
            * (
                math.sin(inc1_rad) * math.sin(azi1_rad)
                + math.sin(inc2_rad) * math.sin(azi2_rad)
            )
            * rf
        )

        delta_n: float = (
            (delta_md / 2)
            * (
                math.sin(inc1_rad) * math.cos(azi1_rad)
                + math.sin(inc2_rad) * math.cos(azi2_rad)
            )
            * rf
        )

        delta_v: float = (delta_md / 2) * (math.cos(inc1_rad) + math.cos(inc2_rad)) * rf

        # Cumulative values
        north: float = results[i - 1]["North"] + delta_n
        east: float = results[i - 1]["East"] + delta_e
        tvd: float = results[i - 1]["TVD"] + delta_v

        # Step 2.6: Calculate Vsec (horizontal displacement)
        vsec: float = math.sqrt(east**2 + north**2)

        # Step 2.7: Calculate ∑DLS
        sum_dls: float = sum(r["DLS"] for r in results) + dls

        # Step 2.8: Calculate Tortuosity
        tortuosity: float = sum_dls / i if i > 0 else 0.0

        # Step 2.9: Calculate AHD (Actual Horizontal Displacement)
        ahd: float = math.sqrt(delta_e**2 + delta_n**2)

        # Step 2.11: Calculate DDI
        if tvd > 0 and ahd > 0 and tortuosity > 0:
            ddi: float = math.log10((ahd * md2 * tortuosity) / tvd)
        else:
            ddi = 0.0

        # Build and turn rates
        if delta_md > 0:
            build_rate: float = (inc2 - inc1) / delta_md * 100  # degrees per 100ft
            turn_rate: float = (azi2 - azi1) / delta_md * 100  # degrees per 100ft
        else:
            build_rate = 0.0
            turn_rate = 0.0

        # Closure distance and azimuth
        closure_distance: float = math.sqrt(north**2 + east**2)
        if abs(north) > 1e-10 or abs(east) > 1e-10:
            closure_azimuth: float = math.atan2(east, north) * 180 / math.pi
            if closure_azimuth < 0:
                closure_azimuth += 360
        else:
            closure_azimuth = 0.0

        # Store results
        row: Dict[str, float] = {
            "MD": md2,
            "Inc": inc2,
            "Azi": azi2,
            "TVD": tvd,
            "North": north,
            "East": east,
            "DLS": dls,
            "Build_Rate": build_rate,
            "Turn_Rate": turn_rate,
            "Vsec": vsec,
            "Closure_Distance": closure_distance,
            "Closure_Azimuth": closure_azimuth,
            "DDI": ddi,
        }
        results.append(row)

    return pl.DataFrame(results)


def create_3d_trajectory_plot(df: pl.DataFrame) -> Figure:
    """
    Create 3D trajectory visualization.

    Args:
        df: DataFrame with trajectory parameters

    Returns:
        Plotly 3D scatter plot figure
    """
    fig = go.Figure()

    # Add trajectory line
    fig.add_trace(
        go.Scatter3d(
            x=df["East"].to_list(),
            y=df["North"].to_list(),
            z=df["TVD"].to_list(),
            mode="lines+markers",
            name="Wellbore Trajectory",
            line=dict(color="blue", width=4),
            marker=dict(size=3, color="red"),
        )
    )

    # Customize layout
    fig.update_layout(
        title="3D Wellbore Trajectory",
        scene=dict(
            xaxis_title="East (ft)",
            yaxis_title="North (ft)",
            zaxis_title="TVD (ft)",
            zaxis=dict(autorange="reversed"),  # TVD increases downward
        ),
        width=800,
        height=600,
    )

    return fig


def create_plan_view_plot(df: pl.DataFrame) -> Figure:
    """
    Create plan view (top-down) plot.

    Args:
        df: DataFrame with trajectory parameters

    Returns:
        Plotly 2D scatter plot figure
    """
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df["East"].to_list(),
            y=df["North"].to_list(),
            mode="lines+markers",
            name="Wellbore Path",
            line=dict(color="blue", width=3),
            marker=dict(size=6, color="red"),
        )
    )

    fig.update_layout(
        title="Plan View (Top-Down)",
        xaxis_title="East (ft)",
        yaxis_title="North (ft)",
        width=400,
        height=400,
        showlegend=False,
    )

    return fig


def create_vertical_section_plot(df: pl.DataFrame) -> Figure:
    """
    Create vertical section view.

    Args:
        df: DataFrame with trajectory parameters

    Returns:
        Plotly 2D scatter plot figure
    """
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df["Vsec"].to_list(),
            y=df["TVD"].to_list(),
            mode="lines+markers",
            name="Vertical Section",
            line=dict(color="green", width=3),
            marker=dict(size=6, color="red"),
        )
    )

    fig.update_layout(
        title="Vertical Section View",
        xaxis_title="Horizontal Displacement (ft)",
        yaxis_title="TVD (ft)",
        yaxis=dict(autorange="reversed"),  # TVD increases downward
        width=400,
        height=400,
        showlegend=False,
    )

    return fig


def create_ddi_plot(df: pl.DataFrame) -> Figure:
    """
    Create DDI vs MD plot.

    Args:
        df: DataFrame with trajectory parameters

    Returns:
        Plotly 2D scatter plot figure
    """
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df["MD"].to_list(),
            y=df["DDI"].to_list(),
            mode="lines+markers",
            name="Dogleg Deviation Index",
            line=dict(color="purple", width=3),
            marker=dict(size=4, color="red"),
        )
    )

    fig.update_layout(
        title="DDI (Dogleg Deviation Index) vs Measured Depth",
        xaxis_title="Measured Depth (ft)",
        yaxis_title="DDI",
        width=400,
        height=300,
        showlegend=False,
    )

    return fig


def create_dls_plot(df: pl.DataFrame) -> Figure:
    """
    Create DLS vs MD plot.

    Args:
        df: DataFrame with trajectory parameters

    Returns:
        Plotly 2D scatter plot figure
    """
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df["MD"].to_list(),
            y=df["DLS"].to_list(),
            mode="lines+markers",
            name="Dogleg Severity",
            line=dict(color="orange", width=3),
            marker=dict(size=4, color="red"),
        )
    )

    fig.update_layout(
        title="Dogleg Severity vs Measured Depth",
        xaxis_title="Measured Depth (ft)",
        yaxis_title="DLS (°/100ft)",
        width=400,
        height=300,
        showlegend=False,
    )

    return fig


def load_sample_data() -> str:
    """
    Generate sample survey data for testing.

    Returns:
        CSV string with sample data
    """
    return """MD,Inc,Azi
0,0,0
500,2,45
1000,15,50
1500,25,55
2000,30,60
2500,32,65
3000,35,70"""


def process_uploaded_file(uploaded_file: Any) -> Optional[str]:
    """
    Process uploaded file and return its content as string.

    Args:
        uploaded_file: Streamlit uploaded file object

    Returns:
        File content as string, or None if processing failed
    """
    try:
        return StringIO(uploaded_file.getvalue().decode("utf-8")).read()
    except UnicodeDecodeError as e:
        st.error(f"Error reading uploaded file: {str(e)}")
        return None


def clean_and_parse_csv(csv_content: str) -> Optional[pl.DataFrame]:
    """
    Clean and parse CSV content, handling various formats and delimiters.

    Args:
        csv_content: Raw CSV content as string

    Returns:
        Parsed DataFrame or None if parsing failed
    """
    try:
        # First, try standard comma-separated format
        try:
            df = pl.read_csv(
                StringIO(csv_content),
                separator=",",
                truncate_ragged_lines=True,
                ignore_errors=True,
            )
            # Check if we have the required columns
            if all(col in df.columns for col in ["MD", "Inc", "Azi"]):
                return df
        except pl.errors.ComputeError:
            # Ignore parsing errors for this format
            pass

        # Try semicolon-separated format (European style)
        try:
            df = pl.read_csv(
                StringIO(csv_content),
                separator=";",
                truncate_ragged_lines=True,
                ignore_errors=True,
            )

            # Handle European decimal format (comma as decimal separator)
            # and clean up column names
            if len(df.columns) >= 3:
                # Take first 3 columns and rename them
                df = df.select(df.columns[:3])
                df.columns = ["MD", "Inc", "Azi"]

                # Convert European decimals (comma) to standard format (dot)
                for col in ["MD", "Inc", "Azi"]:
                    try:
                        # Convert to string, replace comma with dot, then to float
                        df = df.with_columns(
                            pl.col(col)
                            .cast(pl.Utf8)
                            .str.replace(",", ".")
                            .cast(pl.Float64)
                        )
                    except ValueError:
                        # If conversion fails, try direct casting
                        df = df.with_columns(pl.col(col).cast(pl.Float64))

                return df
        except pl.errors.ComputeError as e:
            st.error(f"Error parsing semicolon-separated CSV: {str(e)}")

        # Try tab-separated format
        try:
            df = pl.read_csv(
                StringIO(csv_content),
                separator="\t",
                truncate_ragged_lines=True,
                ignore_errors=True,
            )
            if all(col in df.columns for col in ["MD", "Inc", "Azi"]):
                return df
        except pl.errors.ComputeError:
            pass

        return None

    except pl.errors.ComputeError as e:
        st.error(f"Error processing CSV: {str(e)}")
        return None


def display_interpolation_tool(df_results: pl.DataFrame) -> None:
    """
    Display interpolation tool UI and handle interpolation requests.

    Args:
        df_results: DataFrame with calculated trajectory parameters
    """
    st.subheader("🔧 Survey Point Interpolation")
    st.markdown(
        "Calculate survey parameters at any measured depth using minimum curvature method"
    )

    # Get MD range for input validation
    md_min = float(df_results["MD"].min())  # type: ignore
    md_max = float(df_results["MD"].max())  # type: ignore

    col1, col2 = st.columns([2, 1])

    with col1:
        target_md = st.number_input(
            "Target Measured Depth (ft)",
            min_value=md_min,
            max_value=md_max,
            value=(md_min + md_max) / 2,
            step=10.0,
            help=f"Enter MD between {md_min:.1f} and {md_max:.1f} ft",
            key="interpolation_target_md",
        )

    with col2:
        if st.button("📍 Interpolate", type="primary", key="interpolate_button"):
            # Perform interpolation
            result = interpolate_survey_point(df_results, target_md)

            if result:
                st.session_state.interpolation_result = result
                st.session_state.show_interpolated_point = True
            else:
                st.error(
                    "❌ Could not interpolate at this MD. Check if MD is within survey range."
                )

    # Display interpolation results
    if (
        hasattr(st.session_state, "interpolation_result")
        and st.session_state.interpolation_result
    ):
        result = st.session_state.interpolation_result

        st.success(f"✅ Interpolated survey point at MD {result['MD']:.1f} ft")

        # Display results in a nice format
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Inclination", f"{result['Inc']:.2f}°")
            st.metric("TVD", f"{result['TVD']:.1f} ft")

        with col2:
            st.metric("Azimuth", f"{result['Azi']:.2f}°")
            st.metric("North", f"{result['North']:.1f} ft")

        with col3:
            st.metric("DLS", f"{result['DLS']:.2f}°/100ft")
            st.metric("East", f"{result['East']:.1f} ft")

        with col4:
            st.metric("DDI", f"{result['DDI']:.3f}")
            closure_dist = math.sqrt(result["North"] ** 2 + result["East"] ** 2)
            st.metric("Closure", f"{closure_dist:.1f} ft")

        # Option to add to visualization
        if st.button("📊 Add to 3D Plot", key="add_to_plot_button"):
            st.session_state.show_interpolated_point = True
            st.rerun()


def display_key_metrics(df_results: pl.DataFrame) -> None:
    """
    Display key trajectory metrics in columns.

    Args:
        df_results: DataFrame with calculated results
    """
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total MD", f"{df_results['MD'].max():.0f} ft")
    with col2:
        st.metric("Final TVD", f"{df_results['TVD'].max():.0f} ft")
    with col3:
        st.metric("Max DLS", f"{df_results['DLS'].max():.1f} °/100ft")
    with col4:
        st.metric("Final DDI", f"{df_results['DDI'].max():.2f}")


def display_instructions() -> None:
    """Display usage instructions in an expandable section."""
    with st.expander("📖 Instructions", expanded=True):
        st.markdown("""
        ### How to use this calculator:

        1. **Prepare your data**: CSV file with columns MD, Inc, Azi
           - **MD**: Measured Depth (feet)
           - **Inc**: Inclination (degrees)
           - **Azi**: Azimuth (degrees)

        2. **Upload or use sample**: Click 'Load Sample Data' for a quick demo

        3. **View results**:
           - 3D trajectory visualization
           - Calculated parameters (TVD, North, East, DLS, DDI, etc.)
           - Plan and vertical section views
           - Export results to CSV

        ### Calculations performed:
        - ✅ Minimum Curvature Method
        - ✅ True Vertical Depth (TVD)
        - ✅ North/East coordinates
        - ✅ Dogleg Severity (DLS)
        - ✅ Build/Turn rates
        - ✅ DDI (Dogleg Deviation Index)
        - ✅ Closure distance and azimuth
        """)


def main() -> None:
    """Main Streamlit application function."""
    st.title("🛢️ Survey Data Trajectory Calculator")
    st.markdown(
        "**Upload CSV survey data (MD, Inc, Azi) and get instant trajectory calculations and visualizations**"
    )

    # Sidebar for file upload
    with st.sidebar:
        st.header("📁 Data Input")

        # Sample data option
        if st.button("Load Sample Data"):
            st.session_state.uploaded_data = load_sample_data()

        # File upload
        uploaded_file = st.file_uploader(
            "Upload CSV file",
            type=["csv"],
            help="CSV with columns: MD, Inc, Azi. Supports comma, semicolon, or tab separators. European decimal format (comma) supported.",
        )

        if uploaded_file is not None:
            file_content = process_uploaded_file(uploaded_file)
            if file_content:
                st.session_state.uploaded_data = file_content

    # Process data if available
    if hasattr(st.session_state, "uploaded_data"):
        try:
            # Clean and parse CSV data
            df_input: Optional[pl.DataFrame] = clean_and_parse_csv(
                st.session_state.uploaded_data
            )

            if df_input is None:
                st.error("❌ Could not parse CSV file. Please check the format.")
                st.info(
                    "Supported formats: comma-separated, semicolon-separated, or tab-separated CSV files with columns MD, Inc, Azi"
                )
                return

            st.success(f"✅ Loaded {len(df_input)} survey points")

            # Calculate trajectory parameters
            with st.spinner("Calculating trajectory parameters..."):
                df_results: pl.DataFrame = calculate_trajectory_parameters(df_input)

            # Display input data with interpolation tool
            with st.expander("📊 Input Survey Data", expanded=False):
                st.dataframe(df_input.to_pandas(), use_container_width=True)

                # Add interpolation tool
                display_interpolation_tool(df_results)

            # Calculate trajectory parameters
            with st.spinner("Calculating trajectory parameters..."):
                df_results: pl.DataFrame = calculate_trajectory_parameters(df_input)

            # Display results in tabs
            tab1, tab2, tab3, tab4 = st.tabs(
                ["📈 3D Trajectory", "📊 Results Table", "📉 Charts", "💾 Export"]
            )

            with tab1:
                st.subheader("3D Wellbore Trajectory")
                fig_3d: Figure = create_3d_trajectory_plot(df_results)
                st.plotly_chart(fig_3d, use_container_width=True)

            with tab2:
                st.subheader("Calculated Trajectory Parameters")

                # Format numbers for better display
                df_display = df_results.to_pandas().round(2)
                st.dataframe(df_display, use_container_width=True)

                # Key statistics
                display_key_metrics(df_results)

            with tab3:
                st.subheader("Trajectory Analysis Charts")

                col1, col2 = st.columns(2)
                with col1:
                    fig_plan: Figure = create_plan_view_plot(df_results)
                    st.plotly_chart(fig_plan, use_container_width=True)

                    fig_dls: Figure = create_dls_plot(df_results)
                    st.plotly_chart(fig_dls, use_container_width=True)

                with col2:
                    fig_vertical: Figure = create_vertical_section_plot(df_results)
                    st.plotly_chart(fig_vertical, use_container_width=True)

                    fig_ddi: Figure = create_ddi_plot(df_results)
                    st.plotly_chart(fig_ddi, use_container_width=True)

            with tab4:
                st.subheader("Export Results")

                # Convert to CSV
                csv_output: str = df_results.to_pandas().to_csv(index=False)

                st.download_button(
                    label="📥 Download Results as CSV",
                    data=csv_output,
                    file_name="trajectory_results.csv",
                    mime="text/csv",
                )

                st.info(
                    "💡 **Tip**: Use the downloaded CSV for further analysis or import into other software"
                )

        except (
            pl.errors.ComputeError,
            ValueError,
            KeyError,
        ) as e:
            st.error(f"❌ Error processing data: {str(e)}")
            st.info("Please ensure your CSV has columns: MD, Inc, Azi (or equivalent)")
            st.info(
                "Supported formats: comma-separated, semicolon-separated, or tab-separated"
            )

    else:
        st.info("👆 Upload a CSV file or click 'Load Sample Data' to get started")
        display_instructions()


if __name__ == "__main__":
    main()
