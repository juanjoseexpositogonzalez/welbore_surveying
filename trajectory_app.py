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
                "T": float(df_results["T"][i]),
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

    # Calculate T parameter (interpolated)
    t1 = float(df_results["T"][lower_idx])
    beta_degrees = beta * 180 / math.pi
    target_t = t1 + beta_degrees

    # Calculate DDI using T parameter
    vsec = math.sqrt(target_east**2 + target_north**2)
    if target_tvd > 0 and vsec > 0 and target_t > 0:
        target_ddi = math.log10((vsec * target_md * target_t) / target_tvd)
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
        "T": target_t,
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
        "AHD": 0.0,
        "Closure_Distance": 0.0,
        "Closure_Azimuth": 0.0,
        "Tortuosity": 0.0,
        "T": 0.0,
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

        # Step 2.4b: Calculate T (cumulative sum of beta in degrees)
        beta_degrees: float = beta * 180 / math.pi
        previous_t: float = results[i - 1]["T"] if i > 0 else 0.0
        t_parameter: float = previous_t + beta_degrees

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

        # Step 2.6: Calculate AHD (Along Hole Depth) - cumulative horizontal displacement
        ahd: float = math.sqrt(east**2 + north**2)

        # Step 2.6b: Calculate Vsec (Vertical Section) - AHD × cos(azimuth increment)
        # Azimuth increment from previous point
        if i > 0:
            previous_azi = float(df["Azi"][i - 1])
            azi_increment = (
                abs(azi2 - previous_azi) * math.pi / 180
            )  # Convert to radians
        else:
            azi_increment = 0.0  # No increment for first point

        vsec: float = ahd * math.cos(azi_increment)

        # Step 2.7: Calculate Tortuosity (corrected formula)
        # Tortuosity = cumulative sum of |ΔMD × ΔDLS| / 100
        # Where ΔDLS is the change in DLS from previous point
        if i > 1:
            previous_dls = results[i - 1]["DLS"]
            delta_dls = abs(dls - previous_dls)
        else:
            delta_dls = abs(dls)  # For first calculation, use current DLS

        delta_tortuosity: float = abs(delta_md * delta_dls) / 100
        previous_tortuosity: float = results[i - 1]["Tortuosity"] if i > 0 else 0.0
        tortuosity: float = previous_tortuosity + delta_tortuosity

        # Step 2.11: Calculate DDI (corrected formula using T parameter)
        # DDI is based on the T parameter (cumulative sum of beta in degrees)
        if tvd > 0 and ahd > 0 and t_parameter > 0:
            ddi: float = max(math.log10((ahd * md2 * t_parameter) / tvd), 0.0)
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
            "AHD": ahd,
            "Closure_Distance": closure_distance,
            "Closure_Azimuth": closure_azimuth,
            "Tortuosity": tortuosity,
            "T": t_parameter,
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
            line={"color": "green", "width": 3},
            marker={"size": 6, "color": "red"},
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
    Create DDI and Tortuosity vs MD plot with dual y-axes.

    Args:
        df: DataFrame with trajectory parameters

    Returns:
        Plotly 2D scatter plot figure with dual y-axes
    """
    fig = go.Figure()

    # Add DDI trace (primary y-axis)
    fig.add_trace(
        go.Scatter(
            x=df["MD"].to_list(),
            y=df["DDI"].to_list(),
            mode="lines+markers",
            name="DDI",
            line={"color": "purple", "width": 3},
            marker={"size": 4, "color": "purple"},
            yaxis="y",
        )
    )

    # Add Tortuosity trace (secondary y-axis)
    fig.add_trace(
        go.Scatter(
            x=df["MD"].to_list(),
            y=df["Tortuosity"].to_list(),
            mode="lines+markers",
            name="Tortuosity",
            line={"color": "orange", "width": 3},
            marker={"size": 4, "color": "orange"},
            yaxis="y2",
        )
    )

    # Update layout with dual y-axes
    fig.update_layout(
        title="DDI & Tortuosity vs Measured Depth",
        xaxis_title="Measured Depth (ft)",
        yaxis=dict(
            title="DDI",
            titlefont=dict(color="purple"),
            tickfont=dict(color="purple"),
            side="left",
        ),
        yaxis2=dict(
            title="Tortuosity",
            titlefont=dict(color="orange"),
            tickfont=dict(color="orange"),
            overlaying="y",
            side="right",
        ),
        width=400,
        height=300,
        showlegend=True,
        legend=dict(
            x=0.02,
            y=0.98,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="Black",
            borderwidth=1,
        ),
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
            line={"color": "orange", "width": 3},
            marker={"size": 4, "color": "red"},
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


def clean_and_parse_csv(csv_content: str) -> Optional[pl.DataFrame]:
    """
    Clean and parse CSV content, handling various formats and delimiters.

    Enhanced robust parsing that supports:
    - Standard comma-separated CSV with dot decimal separator
    - European format CSV with comma decimal separator (quoted or unquoted)
    - Semicolon-separated CSV (prioritized for European files)
    - Tab-separated CSV
    - Files with extra empty columns
    - Mixed number formats and malformed data
    - Automatic column detection and cleaning

    Args:
        csv_content: Raw CSV content as string

    Returns:
        Parsed DataFrame with columns MD, Inc, Azi or None if parsing failed
    """
    parsing_errors: List[str] = []

    def safe_numeric_conversion(df: pl.DataFrame, columns: List[str]) -> pl.DataFrame:
        """Safely convert columns to numeric, handling various formats."""
        for col in columns:
            if col in df.columns:
                try:
                    # Clean the column: remove quotes, replace comma with dot, strip whitespace
                    df = df.with_columns(
                        pl.col(col)
                        .cast(pl.Utf8)  # Ensure string type first
                        .str.strip_chars()  # Remove leading/trailing whitespace
                        .str.strip_chars("\"'")  # Remove quotes
                        .str.replace(
                            ",", "."
                        )  # Replace European decimal comma with dot
                        .str.replace_all(
                            r"[^\d\.\-\+]", "", literal=False
                        )  # Remove non-numeric chars
                    )

                    # Try to convert to float, handling nulls and empty strings
                    df = df.with_columns(
                        pl.col(col).map_elements(
                            lambda x: None
                            if (
                                x is None
                                or str(x).strip() == ""
                                or str(x).strip() == "nan"
                            )
                            else float(str(x)),
                            return_dtype=pl.Float64,
                        )
                    )

                except (
                    pl.exceptions.ComputeError,
                    pl.exceptions.SchemaError,
                    ValueError,
                ) as e:
                    parsing_errors.append(f"Error converting column {col}: {str(e)}")
                    continue
        return df

    def detect_and_clean_columns(df: pl.DataFrame) -> Optional[pl.DataFrame]:
        """Detect MD, Inc, Azi columns and clean the dataframe."""
        if len(df.columns) < 3:
            return None

        # Find the first three columns that might contain data
        data_columns: List[str] = []
        for col in df.columns:
            # Check if column contains mostly non-empty data
            non_null_count = df.select(pl.col(col).is_not_null().sum()).item()
            if non_null_count > len(df) * 0.1:  # At least 10% non-null data
                data_columns.append(col)
            if len(data_columns) == 3:
                break

        if len(data_columns) < 3:
            # Fallback: take first 3 columns regardless
            data_columns = df.columns[:3]

        # Select only the data columns
        df = df.select(data_columns[:3])

        # Rename columns to standard names
        df.columns = ["MD", "Inc", "Azi"]

        # Remove rows where all values are null or empty
        df = df.filter(
            pl.col("MD").is_not_null()
            & pl.col("Inc").is_not_null()
            & pl.col("Azi").is_not_null()
        )

        return df

    try:
        # Strategy 1: Try semicolon-separated format (European style) FIRST
        try:
            df = pl.read_csv(
                StringIO(csv_content),
                separator=";",
                truncate_ragged_lines=True,
                ignore_errors=True,
                null_values=["", "NULL", "null", "N/A", "n/a", "nan"],
            )

            # Clean and detect columns
            df = detect_and_clean_columns(df)
            if df is not None and len(df) > 0:
                # Convert to numeric
                df = safe_numeric_conversion(df, ["MD", "Inc", "Azi"])

                # Validate that we have some valid numeric data
                valid_md = df.select(pl.col("MD").is_not_null().sum()).item() > 0
                valid_inc = df.select(pl.col("Inc").is_not_null().sum()).item() > 0
                valid_azi = df.select(pl.col("Azi").is_not_null().sum()).item() > 0
                if valid_md and valid_inc and valid_azi:
                    return df

        except (pl.exceptions.ComputeError, ValueError) as e:
            parsing_errors.append(f"Semicolon-separated parsing failed: {str(e)}")

        # Strategy 2: Try comma-separated with European decimal handling
        try:
            # Try with explicit schema first
            try:
                df = pl.read_csv(
                    StringIO(csv_content),
                    separator=",",
                    schema={"MD (ft)": pl.Utf8, "Inc (°)": pl.Utf8, "Azi (°)": pl.Utf8},
                    truncate_ragged_lines=True,
                    null_values=["", "NULL", "null", "N/A", "n/a"],
                )
            except (pl.exceptions.ComputeError, pl.exceptions.SchemaError):
                # Fallback: read without schema
                df = pl.read_csv(
                    StringIO(csv_content),
                    separator=",",
                    truncate_ragged_lines=True,
                    null_values=["", "NULL", "null", "N/A", "n/a"],
                )

            # Clean column names - remove units and special characters
            new_columns: List[str] = []
            for col in df.columns:
                clean_col = col.strip()
                if any(
                    keyword in clean_col.upper()
                    for keyword in ["MD", "DEPTH", "MEASURED"]
                ):
                    new_columns.append("MD")
                elif any(
                    keyword in clean_col.upper() for keyword in ["INC", "INCLINATION"]
                ):
                    new_columns.append("Inc")
                elif any(
                    keyword in clean_col.upper() for keyword in ["AZI", "AZIMUTH", "AZ"]
                ):
                    new_columns.append("Azi")
                else:
                    new_columns.append(clean_col)

            df.columns = new_columns

            # Check if we have the required columns
            if all(col in df.columns for col in ["MD", "Inc", "Azi"]):
                df = safe_numeric_conversion(df, ["MD", "Inc", "Azi"])

                # Validate conversion
                valid_md = df.select(pl.col("MD").is_not_null().sum()).item() > 0
                valid_inc = df.select(pl.col("Inc").is_not_null().sum()).item() > 0
                valid_azi = df.select(pl.col("Azi").is_not_null().sum()).item() > 0
                if valid_md and valid_inc and valid_azi:
                    return df

        except (pl.exceptions.ComputeError, ValueError) as e:
            parsing_errors.append(f"Comma-separated parsing failed: {str(e)}")

        # Strategy 3: Try tab-separated format
        try:
            df = pl.read_csv(
                StringIO(csv_content),
                separator="\t",
                truncate_ragged_lines=True,
                ignore_errors=True,
                null_values=["", "NULL", "null", "N/A", "n/a"],
            )

            df = detect_and_clean_columns(df)
            if df is not None and len(df) > 0:
                df = safe_numeric_conversion(df, ["MD", "Inc", "Azi"])

                valid_md = df.select(pl.col("MD").is_not_null().sum()).item() > 0
                valid_inc = df.select(pl.col("Inc").is_not_null().sum()).item() > 0
                valid_azi = df.select(pl.col("Azi").is_not_null().sum()).item() > 0
                if valid_md and valid_inc and valid_azi:
                    return df

        except (pl.exceptions.ComputeError, ValueError) as e:
            parsing_errors.append(f"Tab-separated parsing failed: {str(e)}")

        # If all strategies failed, show detailed error information
        st.error("❌ Could not parse CSV file. Please check the format.")

        with st.expander("🔍 Detailed parsing errors", expanded=False):
            st.write("**Attempted parsing strategies:**")
            for i, error in enumerate(parsing_errors, 1):
                st.write(f"{i}. {error}")

            st.write("\n**CSV content preview:**")
            lines = csv_content.split("\n")[:5]
            for i, line in enumerate(lines):
                st.code(f"Line {i+1}: {line}")

        st.info("""
        **Expected format:** CSV with columns MD, Inc, Azi (or similar)

        **Supported formats:**
        - Standard: `MD,Inc,Azi` with dot decimal separator
        - European: `MD (ft),Inc (°),Azi (°)` with comma decimal separator
        - Semicolon: `MD;Inc;Azi` with comma decimal separator
        - Tab-separated with any decimal format

        **Example of valid European format:**
        ```
        MD (ft),Inc (°),Azi (°)
        0,0,0
        "23,3",0,"69,97"
        ```
        """)

        return None

    except (pl.exceptions.ComputeError, ValueError, UnicodeDecodeError) as e:
        st.error(f"❌ Unexpected error processing CSV: {str(e)}")
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
        - ✅ DDI (Drilling Difficulty Index)
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

        except (pl.exceptions.ComputeError, ValueError, AttributeError, TypeError) as e:
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
