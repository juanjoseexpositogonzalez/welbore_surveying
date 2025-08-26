# 🛢️ Survey Data Trajectory Calculator

A comprehensive web-based application for analyzing wellbore trajectory data using industry-standard minimum curvature method calculations. Built with Streamlit and Polars for high-performance data processing and interactive visualizations.

![Python](https://img.shields.io/badge/python-v3.13+-blue.svg)
![Streamlit](https://img.shields.io/badge/streamlit-v1.48+-red.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

## 📋 Table of Contents

- [Overview](#-overview)
- [Features](#-features)
- [Installation](#-installation)
- [Usage](#-usage)
- [Mathematical Background](#-mathematical-background)
- [Calculation Methodology](#-calculation-methodology)
- [Input Data Format](#-input-data-format)
- [Output Parameters](#-output-parameters)
- [Visualizations](#-visualizations)
- [Technical Architecture](#-technical-architecture)
- [API Reference](#-api-reference)
- [Roadmap](#-roadmap)
- [Contributing](#-contributing)
- [License](#-license)

## 🎯 Overview

The Survey Data Trajectory Calculator is designed to bridge the gap between expensive enterprise drilling software (costing $50K+ annually) and basic spreadsheet tools. It provides professional-grade trajectory analysis capabilities in an intuitive, web-based interface accessible to drilling engineers, consultants, and students.

This application implements industry-standard minimum curvature calculations with enhanced DDI (Dogleg Deviation Index) analysis using corrected mathematical formulas that distinguish between different trajectory complexity measures for improved accuracy.

### Key Benefits
- **Cost-effective**: Free alternative to expensive drilling software
- **Modern interface**: Web-based, no installation required
- **Industry standard**: Uses minimum curvature method (ISCWSA approved)
- **Enhanced accuracy**: Corrected DDI and tortuosity calculations using T parameter
- **Comprehensive analysis**: Calculate 13+ trajectory parameters including AHD and Vsec
- **Interactive visualizations**: 3D plots, dual y-axis DDI/Tortuosity analysis
- **Format flexibility**: Supports European CSV formats (semicolon, comma decimals)
- **Export capabilities**: Results downloadable as CSV

## ✨ Features

### Core Functionality
- ✅ **CSV Data Import**: Drag-and-drop file upload with validation
- ✅ **European CSV Support**: Automatic detection of semicolon delimiters and comma decimals
- ✅ **Minimum Curvature Calculations**: Industry-standard trajectory analysis
- ✅ **Real-time Processing**: Instant calculations and visualizations
- ✅ **Multiple Views**: 3D, plan view, vertical section, and parameter plots
- ✅ **Enhanced DDI Analysis**: Dual y-axis plots with DDI and Tortuosity
- ✅ **Export Results**: Download calculated parameters as CSV
- ✅ **Sample Data**: Built-in test data for immediate evaluation

### Calculated Parameters
| Parameter | Description | Units |
|-----------|-------------|-------|
| **TVD** | True Vertical Depth | ft |
| **North** | North coordinate displacement | ft |
| **East** | East coordinate displacement | ft |
| **DLS** | Dogleg Severity | °/100ft |
| **DDI** | Dogleg Deviation Index | dimensionless |
| **Build Rate** | Inclination change rate | °/100ft |
| **Turn Rate** | Azimuth change rate | °/100ft |
| **Vsec** | Horizontal displacement with azimuth correction | ft |
| **AHD** | Actual Horizontal Displacement (cumulative) | ft |
| **Closure Distance** | Total horizontal offset | ft |
| **Closure Azimuth** | Direction to target | degrees |
| **Tortuosity** | Cumulative trajectory complexity measure | dimensionless |
| **T Parameter** | Cumulative angular deviation for DDI calculation | degrees |

### Visualizations
- 🎯 **3D Trajectory Plot**: Interactive wellbore path visualization
- 📊 **Plan View**: Top-down horizontal projection
- 📈 **Vertical Section**: Side view with TVD
- 📉 **DLS vs MD**: Dogleg severity analysis
- 🟣 **DDI & Tortuosity vs MD**: Dual y-axis plot for trajectory complexity assessment

## 🚀 Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Quick Setup

```bash
# Clone the repository
git clone https://github.com/juanjoseexpositogonzalez/wellbore-surveying.git
cd wellbore-surveying

# Using uv (recommended for fastest setup)
uv sync
uv run poe gui

# OR using traditional pip
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
streamlit run trajectory_app.py
```

### Dependencies

```txt
streamlit>=1.48.0
polars>=1.32.0
plotly>=6.3.0
numpy>=2.3.0
```

## 💻 Usage

### 1. Start the Application
```bash
streamlit run trajectory_app.py
```
The app will open in your default browser at `http://localhost:8501`

### 2. Load Data
**Option A: Sample Data**
- Click "Load Sample Data" for instant testing
- Provides realistic trajectory with 7 survey points

**Option B: Upload CSV**
- Drag and drop your CSV file
- File must contain columns: MD, Inc, Azi
- Supports standard survey data formats

### 3. Analyze Results
- **3D Trajectory**: Interactive 3D wellbore visualization
- **Results Table**: Comprehensive calculated parameters including Tortuosity and T parameter
- **Charts**: Multiple analytical views (plan, vertical, DLS, DDI & Tortuosity dual-axis)
- **Export**: Download results as CSV for further analysis

### 4. Sample Data Format

```csv
MD,Inc,Azi
0,0,0
500,2,45
1000,15,50
1500,25,55
2000,30,60
2500,32,65
3000,35,70
```

## 📐 Mathematical Background

### Minimum Curvature Method

The minimum curvature method is the industry-standard approach for calculating wellbore trajectories, recommended by the ISCWSA (Industry Steering Committee for Wellbore Survey Accuracy). It provides the most accurate representation of the wellbore path between survey stations.

**Key Principles:**
1. **Smooth curvature**: Assumes the wellbore follows a smooth circular arc between survey points
2. **Three-dimensional**: Accounts for both inclination and azimuth changes
3. **Ratio factor**: Applies smoothing to prevent unrealistic sharp turns
4. **Industry acceptance**: Widely adopted across the oil and gas industry

### Coordinate System

The application uses a **right-handed coordinate system**:
- **North**: Positive Y-axis (0° azimuth)
- **East**: Positive X-axis (90° azimuth)
- **TVD**: Positive Z-axis (downward, increasing with depth)
- **Reference**: Surface location at (0, 0, 0)

## 🔢 Calculation Methodology

The application implements the complete minimum curvature workflow following industry best practices:

### Step 1: Data Preparation
```python
# Convert degrees to radians for trigonometric calculations
inc1_rad = inc1 * π / 180
inc2_rad = inc2 * π / 180
azi1_rad = azi1 * π / 180
azi2_rad = azi2 * π / 180
```

### Step 2: Dogleg Angle Calculation (β)
The dogleg angle represents the total angular change between survey stations:

```python
cos_β = cos(inc2 - inc1) - sin(inc1) × sin(inc2) × [1 - cos(azi2 - azi1)]
β = arccos(cos_β)
```

**Physical Meaning**: β quantifies the total 3D angular change, combining both inclination and azimuth changes.

### Step 3: Ratio Factor (RF)
The ratio factor smooths the trajectory calculation:

```python
if β ≈ 0:
    RF = 1.0
else:
    RF = (2/β) × tan(β/2)
```

**Purpose**: Prevents unrealistic sharp angles and provides smooth trajectory curves.

### Step 4: Displacement Calculations
Calculate the incremental displacements between survey stations:

```python
# Course length
ΔMD = MD2 - MD1

# Horizontal displacements
ΔE = (ΔMD/2) × [sin(inc1)×sin(azi1) + sin(inc2)×sin(azi2)] × RF
ΔN = (ΔMD/2) × [sin(inc1)×cos(azi1) + sin(inc2)×cos(azi2)] × RF

# Vertical displacement
ΔV = (ΔMD/2) × [cos(inc1) + cos(inc2)] × RF
```

### Step 5: Cumulative Coordinates
```python
# Build cumulative position
North_i = North_{i-1} + ΔN
East_i = East_{i-1} + ΔE
TVD_i = TVD_{i-1} + ΔV
```

### Step 6: Derived Parameters

**Dogleg Severity (DLS)**
```python
DLS = β × (180/π) × (100/ΔMD)  # degrees per 100ft
```

**Build and Turn Rates**
```python
Build_Rate = (inc2 - inc1) × (100/ΔMD)  # degrees per 100ft
Turn_Rate = (azi2 - azi1) × (100/ΔMD)   # degrees per 100ft
```

**Horizontal Displacement (Vsec)**
```python
# Vsec includes azimuth increment correction
azi_increment = azi2 - azi1
Vsec = AHD * cos(azi_increment * π/180)
```

**Actual Horizontal Displacement (AHD)**
```python
AHD = √(East² + North²)  # Cumulative horizontal displacement
```

**Closure Parameters**
```python
Closure_Distance = √(North² + East²)
Closure_Azimuth = arctan2(East, North) × (180/π)
```

**Dogleg Deviation Index (DDI)**
```python
# Tortuosity calculation (corrected formula)
Tortuosity = Σ(|ΔMD × ΔDLS|) / 100

# T parameter for DDI calculation
T = Σ(β_degrees)  # Cumulative sum of dogleg angles in degrees

# DDI calculation using T parameter
if TVD > 0 and AHD > 0 and T > 0:
    DDI = log₁₀((AHD × MD × T) / TVD)
else:
    DDI = 0
```

**DDI Interpretation:**
- **DDI < 1**: Relatively simple trajectory
- **DDI 1-2**: Moderate complexity
- **DDI > 2**: Complex, tortuous wellbore

## 📊 Input Data Format

### Required Columns

| Column | Description | Units | Range | Notes |
|--------|-------------|-------|-------|-------|
| **MD** | Measured Depth | feet | ≥ 0 | Monotonically increasing |
| **Inc** | Inclination | degrees | 0-180 | 0° = vertical, 90° = horizontal |
| **Azi** | Azimuth | degrees | 0-360 | 0° = North, 90° = East |

### Data Quality Requirements
- **Sorted by MD**: Data must be in ascending measured depth order
- **No gaps**: Consecutive survey points without missing data
- **Realistic values**: Inclination 0-90° for most applications
- **Consistent units**: All measurements in the same unit system
- **Format flexibility**: Supports comma-separated, semicolon-separated (European), and tab-separated files

### CSV Format Example

**Standard Format (Comma-separated):**
```csv
MD,Inc,Azi
0.0,0.0,0.0
500.0,2.1,45.3
1000.0,15.7,50.2
1500.0,25.4,55.1
2000.0,30.8,60.0
```

**European Format (Semicolon-separated with comma decimals):**
```csv
MD;Inc;Azi
0,0;0,0;0,0
500,0;2,1;45,3
1000,0;15,7;50,2
1500,0;25,4;55,1
2000,0;30,8;60,0
```

The application automatically detects and handles both formats, along with tab-separated files.

## 📈 Output Parameters

### Primary Trajectory Data

**TVD (True Vertical Depth)**
- Actual vertical distance from surface
- Critical for geological correlation
- Used in pressure calculations

**North/East Coordinates**
- Horizontal position relative to surface location
- Essential for anti-collision analysis
- Required for well spacing calculations

**DLS (Dogleg Severity)**
- Rate of wellbore direction change
- Critical for drilling equipment limits
- Typical limits: 3-6°/100ft for conventional drilling

**DDI (Dogleg Deviation Index)**
- Quantifies overall trajectory complexity
- Useful for comparing wellbore designs
- Correlates with drilling difficulty and costs

### Secondary Parameters

**AHD (Actual Horizontal Displacement)**
- Cumulative horizontal distance from surface location
- Essential for anti-collision analysis
- Used in DDI calculations

**Vsec (Horizontal Displacement with Azimuth Correction)**
- Azimuth-corrected horizontal displacement
- Accounts for azimuth increment changes
- Different from simple cumulative displacement

**Tortuosity**
- Quantifies overall trajectory complexity using |ΔMD × ΔDLS| accumulation
- More accurate than simple DLS averaging
- Critical component for DDI calculation

**T Parameter**
- Cumulative sum of dogleg angles in degrees
- Required for accurate DDI calculation
- Replaces tortuosity in DDI formula for better accuracy

## 🎨 Visualizations

### 1. 3D Trajectory Plot
- **Interactive rotation and zoom**
- **True spatial representation**
- **Survey point markers**
- **Axis labels with units**

**Features:**
- Mouse controls for 3D navigation
- TVD axis inverted (increases downward)
- Color-coded trajectory line
- Survey station markers

### 2. Plan View (Top-Down)
- **Horizontal projection**
- **Surface facility layout compatibility**
- **Anti-collision reference**

### 3. Vertical Section
- **Side view projection**
- **Geological correlation**
- **Drilling target visualization**

### 4. DLS Analysis
- **Quality control tool**
- **Equipment limitation checking**
- **Drilling performance assessment**

### 5. DDI & Tortuosity Trend
- **Dual y-axis visualization**
- **Complexity assessment**
- **Comparative analysis tool**
- **Drilling difficulty prediction**

## 🏗️ Technical Architecture

### Technology Stack

**Backend Processing:**
- **Polars**: High-performance DataFrame operations
- **NumPy**: Mathematical calculations (minimal usage)
- **Python Math**: Trigonometric functions

**Frontend Interface:**
- **Streamlit**: Web application framework
- **Plotly**: Interactive visualizations
- **HTML/CSS**: UI styling through Streamlit

**Data Flow:**
```
CSV Upload → Polars DataFrame → Calculations → Plotly Visualizations → Export
```

### Performance Characteristics
- **Processing Speed**: <1 second for 100 survey points
- **Memory Usage**: <50MB for typical datasets
- **Scalability**: Tested up to 1000+ survey points
- **Browser Support**: All modern browsers

### Code Quality
- **Type Hints**: Full type annotation coverage
- **Documentation**: Comprehensive docstrings
- **Error Handling**: Robust input validation
- **Testing**: Validated against industry examples

## 📚 API Reference

### Core Functions

#### `calculate_trajectory_parameters(df: pl.DataFrame) -> pl.DataFrame`
Main calculation engine implementing minimum curvature method.

**Parameters:**
- `df`: Input DataFrame with MD, Inc, Azi columns

**Returns:**
- DataFrame with all calculated trajectory parameters

**Raises:**
- `ValueError`: If required columns are missing

#### `create_3d_trajectory_plot(df: pl.DataFrame) -> Figure`
Generate interactive 3D trajectory visualization.

**Parameters:**
- `df`: DataFrame with calculated trajectory parameters

**Returns:**
- Plotly Figure object for 3D visualization

### Utility Functions

#### `degrees_to_radians(degrees: float) -> float`
Convert angle units for trigonometric calculations.

#### `load_sample_data() -> str`
Generate sample CSV data for testing and demonstration.

## 🗺️ Roadmap

### Version 2.0 - Enhanced Analysis
**Target: Q4 2025**

- [ ] **Multi-well Analysis**
  - Import and compare multiple wellbore trajectories
  - Anti-collision analysis with offset wells
  - Batch processing capabilities

- [ ] **Advanced Calculations**
  - Uncertainty ellipsoids (ISCWSA error models)
  - Formation intersection calculations
  - Torque and drag analysis

- [ ] **Data Sources**
  - Direct import from LAS files
  - API integration with drilling data providers
  - Real-time MWD data streaming

### Version 2.1 - Visualization Enhancements
**Target: Q1 2026**

- [ ] **Advanced Plotting**
  - Heat maps for parameter distributions
  - Time-based trajectory animations
  - Geological cross-section overlays

- [ ] **Export Options**
  - PDF report generation
  - PowerPoint slide exports
  - Industry-standard plot templates

- [ ] **Customization**
  - User-defined coordinate systems
  - Custom unit conversions
  - Branding and theming options

### Version 3.0 - Enterprise Features
**Target: Q2 2026**

- [ ] **Database Integration**
  - PostgreSQL/MySQL support
  - Well database management
  - Historical trajectory library

- [ ] **API Development**
  - RESTful API for external integrations
  - Python SDK for programmatic access
  - Webhook support for automated workflows

- [ ] **Advanced Analytics**
  - Machine learning for trajectory optimization
  - Predictive drilling performance models
  - Automated quality control alerts

### Version 3.1 - Collaboration Tools
**Target: Q3 2026**

- [ ] **Team Features**
  - Multi-user project sharing
  - Comment and annotation system
  - Version control for trajectory designs

- [ ] **Integration Platform**
  - Plugin system for custom calculations
  - Third-party software connectors
  - Cloud deployment options

### Long-term Vision
**Target: 2027+**

- [ ] **AI-Powered Analysis**
  - Automated trajectory optimization
  - Intelligent drilling recommendations
  - Predictive maintenance alerts

- [ ] **Mobile Applications**
  - iOS/Android native apps
  - Offline calculation capabilities
  - Field data collection tools

- [ ] **Industry Partnerships**
  - Integration with major drilling contractors
  - Certified calculation algorithms
  - Professional training programs

## 🤝 Contributing

We welcome contributions from the drilling and software engineering communities!

### Development Setup

1. **Fork the repository**
2. **Create feature branch**: `git checkout -b feature/amazing-feature`
3. **Install development dependencies**: `uv sync` or `pip install -r requirements-dev.txt`
4. **Run tests**: `pytest tests/` or `poe test`
5. **Submit pull request**

### Contribution Areas

**High Priority:**
- Additional calculation methods (radius of curvature, average angle)
- Performance optimizations for large datasets
- Enhanced error handling and validation

**Medium Priority:**
- UI/UX improvements and accessibility
- Additional export formats
- Documentation improvements

**Future Features:**
- Integration with external APIs
- Advanced visualization options
- Mobile responsiveness enhancements

### Code Standards
- **Type hints**: All functions must have type annotations
- **Documentation**: Comprehensive docstrings required
- **Testing**: Unit tests for all calculation functions
- **Formatting**: Follow PEP 8 style guidelines

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Third-Party Licenses
- Streamlit: Apache 2.0 License
- Polars: MIT License
- Plotly: MIT License

## 📞 Support

### Getting Help
- **Documentation**: This README and inline code documentation
- **Issues**: GitHub Issues for bug reports and feature requests
- **Discussions**: GitHub Discussions for questions and community support

### Professional Support
For commercial applications and professional consulting:
- **Email**: support@trajectory-calculator.com
- **LinkedIn**: [Professional Profile]
- **Website**: [Project Website]

### Acknowledgments

**Industry Guidance:**
- ISCWSA (Industry Steering Committee for Wellbore Survey Accuracy)
- SPE (Society of Petroleum Engineers) guidelines
- IADC (International Association of Drilling Contractors) standards

**Technical Inspiration:**
- Minimum curvature method implementations
- Industry best practices from drilling professionals
- Open source computational tools community

---

**Built with ❤️ for the drilling engineering community**

*Making professional trajectory analysis accessible to everyone*
