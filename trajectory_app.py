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
            name="Drilling Difficulty Index",
            line={"color": "purple", "width": 3},
            marker={"size": 4, "color": "red"},
        )
    )

    fig.update_layout(
        title="DDI (Drilling Difficulty Index) vs Measured Depth",
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
    except (UnicodeDecodeError, AttributeError) as e:
        st.error(f"Error reading uploaded file: {str(e)}")
        return None


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
            help="CSV should contain columns: MD (Measured Depth), Inc (Inclination), Azi (Azimuth)",
        )

        if uploaded_file is not None:
            file_content = process_uploaded_file(uploaded_file)
            if file_content:
                st.session_state.uploaded_data = file_content

    # Process data if available
    if hasattr(st.session_state, "uploaded_data"):
        try:
            # Read CSV data
            df_input: pl.DataFrame = pl.read_csv(
                StringIO(st.session_state.uploaded_data)
            )

            st.success(f"✅ Loaded {len(df_input)} survey points")

            # Display input data
            with st.expander("📊 Input Survey Data", expanded=False):
                st.dataframe(df_input.to_pandas(), use_container_width=True)

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

        except Exception as e:
            st.error(f"❌ Error processing data: {str(e)}")
            st.info("Please ensure your CSV has columns: MD, Inc, Azi")

    else:
        st.info("👆 Upload a CSV file or click 'Load Sample Data' to get started")
        display_instructions()


if __name__ == "__main__":
    main()
