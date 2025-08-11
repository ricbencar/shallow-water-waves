# Shallow Water Waves Model (CLI & GUI)

## Overview

This project presents a robust solution for accurately computing local shallow-foreshore wave-height distribution parameters, critical for coastal engineering applications. It offers both a **Command-Line Interface (CLI)** and a **Graphical User Interface (GUI)**, providing flexibility for users to interact with the model. The core of this project lies in its implementation of the **Composite Weibull distribution**, a sophisticated model specifically designed for shallow water environments.

Unlike deep-water conditions where wave behavior is often well-described by simpler distributions like Rayleigh, shallow foreshores introduce complex phenomena such as depth-induced breaking, non-linear wave-wave interactions, and bottom friction. These processes significantly alter wave height distributions, making accurate prediction challenging (Groenendijk, 1998; Battjes & Groenendijk, 2000). The Composite Weibull distribution addresses these complexities by employing a two-part structure that effectively captures the distinct physical regimes of unbroken and breaking waves.

Through these applications, users can input essential wave parameters to perform local wave heights distribution calculations. The output includes a detailed report of the computed values, offering insights for the design of marine structures, coastal protection systems, and the assessment of wave-related phenomena like overtopping.
<img width="1083" height="463" alt="shallow-water-waves" src="https://github.com/user-attachments/assets/839a6c12-46d9-4470-b1f9-92394d305c43" />

## Composite Weibull Distribution

The Composite Weibull distribution is a two-part distribution specifically designed for shallow water environments, capturing the distinct physical regimes of unbroken and breaking waves. It is described by the following cumulative distribution function (CDF) for normalized wave heights ($\tilde{h} = h / H_{rms}$):

```math
F(\tilde{h}) = \begin{cases}
1 - \exp\left[ - \left( \frac{\tilde{h}}{\tilde{H}_1} \right)^{k_1} \right], & \tilde{h} < \tilde{H}_{tr} \\
1 - \exp\left[ - \left( \frac{\tilde{h}}{\tilde{H}_2} \right)^{k_2} \right], & \tilde{h} \geq \tilde{H}_{tr}
\end{cases}
```

Where:

* $\tilde{h} = h / H_{rms}$ is the normalized wave height.
* $\tilde{H}_1 = H_1 / H_{rms}$ is the normalized scale parameter for the first part of the distribution (unbroken waves).
* $\tilde{H}_2 = H_2 / H_{rms}$ is the normalized scale parameter for the second part of the distribution (breaking waves).
* $\tilde{H}_{tr} = H_{tr} / H_{rms}$ is the dimensionless transitional wave height, marking the boundary between the two parts of the distribution.
* $k_1 = 2.0$ is the exponent (shape parameter) for the first part, which is Rayleigh-shaped.
* $k_2 = 3.6$ is the empirically determined exponent for the second part.

The parameters $\tilde{H}_1$ and $\tilde{H}_2$ are determined by solving a system of non-linear equations to ensure consistency with the normalized $H_{rms}$ and continuity at the transitional wave height.

## Features

* **Dual Interface**: Offers both a command-line interface for quick calculations and scripting, and a graphical user interface for ease of use.
* **Composite Weibull Distribution Model**: Implements a robust model for shallow-foreshore wave-height distribution with empirically determined exponents ($k_1=2.0$, $k_2=3.6$).
* **Key Parameter Calculation**: Takes free-surface variance ($m_0$) as input and computes the mean square wave height ($H_{rms}$), and dimensional/dimensionless transitional wave heights ($H_{tr\_dim}$, $\tilde{H}_{tr}$).
* **Dimensionless Wave-Height Ratios**: Calculates $\tilde{H}_N$ and $\tilde{H}_{1/N}$ for various characteristic wave heights by solving a system of two non-linear equations derived from the Composite Weibull distribution.
* **Dimensional Wave Heights**: Converts dimensionless ratios to actual wave heights in meters.
* **Diagnostic Ratios**: Provides insights into the wave height distribution shape through various diagnostic ratios.
* **Detailed Reporting**: Generates a comprehensive report of input parameters, calculated values, and diagnostic ratios.

## Input Parameters

Both the CLI and GUI applications require the following four input parameters:

* **$H_{m0}$ (in meters)**: The local significant spectral wave height (used for reporting and context).
* **$m_0$ (in m²)**: The free-surface variance. This is the primary input for the energy calculations.
* **$d$ (in meters)**: The local water depth.
* **Beach slope ($1:m$)**: The beach slope expressed as "1:m". For example, entering 100 signifies a slope of $1/100=0.01$.

## Computational Process

Based on the provided inputs, the program performs a series of calculations to determine various wave height parameters and distribution characteristics:

### 1. Mean Square Wave Height ($H_{rms}$)

The mean square wave height ($H_{rms}$) is an important characteristic of the wave field. It is now calculated directly from the input free-surface variance ($m_0$) and local depth ($d$). This formula, used in both the CLI and GUI implementations, incorporates empirical coefficients to better capture shallow-water effects:

```math
H_{rms} = \left(2.69 + 3.24 \cdot \frac{\sqrt{m_0}}{d}\right) \cdot \sqrt{m_0}
```

This relationship demonstrates an increase with the degree of saturation ($\Psi = \sqrt{m_0}/d$), counteracting the bandwidth effect observed in deep water.

### 2. Dimensional Transitional Wave Height ($H_{tr\_dim}$)

The dimensional transitional wave height ($H_{tr\_dim}$) marks the point where the wave height distribution significantly changes due to depth-induced breaking. It is calculated using the local water depth ($d$) and the beach slope ($m$):

```math
\tan(\alpha) = \frac{1}{m}
```

Then, $H_{tr\_dim}$ is computed as:

```math
H_{tr\_dim} = (0.35 + 5.8 \cdot \tan(\alpha)) \cdot d
```

### 3. Subsequent Calculations

The remaining calculations for the dimensionless transitional parameter ($\tilde{H}_{tr}$), dimensionless wave-height ratios, dimensional wave heights, and diagnostic ratios proceed as described in the original documentation, using the newly computed $H_{rms}$ as the basis for normalization.

## Building and Running

### Command-Line Interface (CLI)

The CLI application (`shallow-water-waves_cli.cpp`) can be compiled using g++.

**Compilation Instructions (Windows example with g++):**

```bash
g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic -static -o shallow-water-waves_cli shallow-water-waves_cli.cpp
```

**Usage:**
You can run the CLI application by providing the parameters as command-line arguments or by entering them interactively.

* **With command-line arguments (e.g., Hm0=2.5, m0=0.3906, d=5.0, slopeM=100):**
    ```bash
    shallow-water-waves_cli 2.5 0.3906 5.0 100
    ```

* **Interactive input:**
    ```bash
    shallow-water-waves_cli
    ```
    (The program will then prompt you for the four required inputs.)

### Graphical User Interface (GUI)

The GUI application (`shallow-water-waves_gui.cpp`) is implemented using the native Win32 API.

**Compilation Instructions (Windows example with g++):**

```bash
g++ -O3 -march=native -std=c++17 -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui -mwindows -static
```

**Usage:**
Run the compiled executable (`shallow-water-waves_gui.exe`). A window will appear where you can input the values and click "Compute" to see the results. The report will also be saved to `report.txt`.

## Utility Script: m0_calculator.py

For convenience, a Python script `m0_calculator.py` is included. This utility can be used to estimate the free-surface variance ($m_0$) if it is not known directly. The script uses a more comprehensive physical model, considering wave shoaling and breaking (based on Thornton & Guza, 1983) to derive $m_0$ from offshore wave parameters like peak period and wave height. It serves as a useful pre-processing tool to generate the necessary $m_0$ input for the main C++ applications.

## References

* **Battjes, J. A., & Groenendijk, H. W. (2000).** Wave height distributions on shallow foreshores. *Coastal Engineering*, 40(3), 161-182.
* **Groenendijk, H. W. (1998).** *Shallow foreshore wave height statistics*. M.Sc.-thesis, Delft University of Technology, Department of Civil Engineering, Section Fluid Mechanics, The Netherlands.
* **Longuet-Higgins, M. S. (1980).** On the distribution of the heights of sea waves: Some effects of nonlinearity and finite band width. *Journal of Geophysical Research*, 85(C3), 1519-1523.
* **Thornton, E. B., & Guza, R. T. (1982).** Energy saturation and phase speeds measured on a natural beach. *Journal of Geophysical Research*, 87(C12), 9499-9508.
* **Thornton, E. B., & Guza, R. T. (1983).** Transformation of wave height distribution. *Journal of Geophysical Research*, 88(C10), 5925-5938.
