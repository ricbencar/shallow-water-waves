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

* ` ```math \tilde{h} = h / H_{rms} ``` is the normalized wave height.
* ` ```math \tilde{H}_1 = H_1 / H_{rms} ``` is the normalized scale parameter for the first part of the distribution (unbroken waves).
* ` ```math \tilde{H}_2 = H_2 / H_{rms} ``` is the normalized scale parameter for the second part of the distribution (breaking waves).
* ` ```math \tilde{H}_{tr} = H_{tr} / H_{rms} ``` is the dimensionless transitional wave height, marking the boundary between the two parts of the distribution.
* ` ```math k_1 = 2.0 ``` is the exponent (shape parameter) for the first part, which is Rayleigh-shaped.
* ` ```math k_2 = 3.6 ``` is the empirically determined exponent for the second part.

The parameters ` ```math \tilde{H}_1 ``` and ` ```math \tilde{H}_2 ``` are determined by solving a system of non-linear equations to ensure consistency with the normalized ` ```math H_{rms} ``` and continuity at the transitional wave height. These equations are:

1.  **Normalized ` ```math H_{rms} ``` Constraint (from Groenendijk, 1998, Equation 7.11):**
    ```math
    \tilde{H}_{rms} = 1 = \sqrt{
    \tilde{H}_1^{2} \, \gamma\left( \frac{2}{k_1} + 1, \left( \frac{\tilde{H}_{tr}}{\tilde{H}_1} \right)^{k_1} \right)
    + \tilde{H}_2^{2} \, \Gamma\left( \frac{2}{k_2} + 1, \left( \frac{\tilde{H}_{tr}}{\tilde{H}_2} \right)^{k_2} \right)
    }
    ```
    This equation ensures that the overall root-mean-square of the normalized composite distribution equals one.

2.  **Continuity Condition (from Groenendijk, 1998, Equation 3.4):**
    ```math
    \left( \frac{\tilde{H}_{tr}}{\tilde{H}_1} \right)^{k_1} = \left( \frac{\tilde{H}_{tr}}{\tilde{H}_2} \right)^{k_2}
    ```
    This condition ensures that the cumulative distribution function is continuous at the transitional wave height ` ```math \tilde{H}_{tr} ```.

## Features

* **Dual Interface**: Offers both a command-line interface for quick calculations and scripting, and a graphical user interface for ease of use.
* **Composite Weibull Distribution Model**: Implements a robust model for shallow-foreshore wave-height distribution with empirically determined exponents (` ```math k_1=2.0 ```, ` ```math k_2=3.6 ```). This two-part distribution is designed to reflect the different physical regimes governing smaller (unbroken) and larger (breaking) waves.
* **Key Parameter Calculation**: Computes mean square wave height (` ```math H_{rms} ```), and dimensional/dimensionless transitional wave heights (` ```math H_{tr\_dim} ```, ` ```math \tilde{H}_{tr} ```) from the free-surface variance (` ```math m_0 ```).
* **Dimensionless Wave-Height Ratios**: Calculates ` ```math \tilde{H}_N ``` (wave height with ` ```math 1/N ``` exceedance probability) and ` ```math \tilde{H}_{1/N} ``` (mean of the highest ` ```math 1/N ```-part of wave heights) for various characteristic wave heights. These are determined by solving a system of two non-linear equations derived from the Composite Weibull distribution, ensuring the normalized ` ```math H_{rms} ``` of the distribution equals one. The solution employs a robust numerical strategy using a Newton-Raphson matrix method for simultaneous root-finding, and precise implementations of incomplete gamma functions.
* **Dimensional Wave Heights**: Converts dimensionless ratios to actual wave heights in meters.
* **Diagnostic Ratios**: Provides insights into the wave height distribution shape through various diagnostic ratios.
* **Detailed Reporting**: Generates a comprehensive report of input parameters, calculated values, and diagnostic ratios.

## Input Parameters

Both the CLI and GUI applications require the following four input parameters:

* ` ```math H_{m0} ``` (in meters): The local significant spectral wave height. This value is used for reporting and context but is not a direct input to the core calculation, which now starts from ` ```math m_0 ```.
* ` ```math m_0 ``` (in m²): The free-surface variance (zeroth spectral moment). This is now the primary input representing the total wave energy.
* ` ```math d ``` (in meters): The local water depth.
* **Beach slope (` ```math 1:m ```)**: The beach slope expressed as "` ```math 1:m ```". For example, entering 20 signifies a slope of ` ```math 1/20=0.05 ```. This parameter influences the transitional wave height.

## Computational Process

Based on the provided inputs, the program performs a series of calculations to determine various wave height parameters and distribution characteristics:

### 1. Mean Square Wave Height (` ```math H_{rms} ```)

The mean square wave height (` ```math H_{rms} ```) is an important characteristic of the wave field. Its calculation now starts directly from the input free-surface variance (` ```math m_0 ```) and incorporates empirical coefficients to better capture the shallow-water distribution of extreme waves. This formula, used in both the CLI and GUI implementations, is:

```math
H_{rms} = \left(2.69 + 3.24 \cdot \frac{\sqrt{m_0}}{d}\right) \cdot \sqrt{m_0}
```

This relationship demonstrates an increase with the degree of saturation (` ```math \Psi = \sqrt{m_0}/d ```), counteracting the bandwidth effect observed in deep water.

### 2. Dimensional Transitional Wave Height (` ```math H_{tr\_dim} ```)

The dimensional transitional wave height (` ```math H_{tr\_dim} ```) marks the point where the wave height distribution significantly changes due to depth-induced breaking. It is calculated using the local water depth (` ```math d ```) and the beach slope (` ```math m ```):
The tangent of the beach slope (` ```math \tan(\alpha) ```) is derived from the input ` ```math m ```:

```math
\tan(\alpha) = \frac{1}{m}
```

Then, ` ```math H_{tr\_dim} ``` is computed as:

```math
H_{tr\_dim} = (0.35 + 5.8 \cdot \tan(\alpha)) \cdot d
```

For example, if ` ```math m=20 ```, then ` ```math \tan(\alpha)=1/20=0.05 ```, and ` ```math H_{tr\_dim}=(0.35+5.8 \cdot 0.05) \cdot d=0.64 \cdot d ```. This relationship indicates that steeper slopes tend to result in higher ` ```math H_{tr} ``` values, implying that fewer waves deviate from the Rayleigh distribution on steeper foreshores.

### 3. Dimensionless Transitional Parameter (` ```math \tilde{H}_{tr} ```)

The dimensionless transitional parameter (` ```math \tilde{H}_{tr} ```) normalizes the dimensional transitional wave height by the mean square wave height:

```math
\tilde{H}_{tr} = \frac{H_{tr\_dim}}{H_{rms}}
```

### 4. Dimensionless Wave-Height Ratios (` ```math \tilde{H}_N ``` and ` ```math \tilde{H}_{1/N} ```)

The dimensionless wave-height ratios are critical outputs of the model. The calculation involves solving a system of two non-linear equations derived from the Composite Weibull distribution, ensuring that the normalized ` ```math H_{rms} ``` of the distribution equals one. This is achieved using a Newton-Raphson matrix method for simultaneous root-finding.

The core of this calculation is finding the values of ` ```math \tilde{H}_1 ``` and ` ```math \tilde{H}_2 ``` that satisfy the normalized ` ```math H_{rms} ``` equation (Equation 7.11 from Groenendijk, 1998) and the continuity condition between the two Weibull distributions (Equation 3.4).

### 5. Dimensional Wave Heights (` ```math H ```)

The calculated dimensionless wave-height ratios (` ```math \tilde{H}_N ``` or ` ```math \tilde{H}_{1/N} ```) are then converted back to dimensional wave heights (in meters) by multiplying them by the mean square wave height (` ```math H_{rms} ```):

```math
H = \tilde{H} \cdot H_{rms}
```

### 6. Diagnostic Ratios

Finally, the program computes several diagnostic ratios, which provide insights into the shape of the wave height distribution and the relative significance of extreme waves. These include characteristic wave height ratios ` ```math (H_{1/N})/(H_{1/3}) ``` with ` ```math N = 10, 50, 100, 250, \text{ and } 1000 ```.

## Supporting Mathematical Functions

The core calculations rely on precise implementations of fundamental mathematical functions:

* **Complete Gamma Function (` ```math \Gamma(z) ```):**
    ```math
    \Gamma(a) = \int_0^{\infty} t^{a-1} e^{-t} dt \quad (a > 0)
    ```
* **Unnormalized Lower Incomplete Gamma Function (` ```math \gamma(a,x) ```):**
    ```math
    \gamma(a, x) = \int_0^x t^{a-1} e^{-t} dt
    ```
* **Unnormalized Upper Incomplete Gamma Function (` ```math \Gamma(a,x) ```):**
    ```math
    \Gamma(a, x) = \int_x^{\infty} t^{a-1} e^{-t} dt
    ```

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

For convenience, a Python script `m0_calculator.py` is included. This utility can be used to estimate the free-surface variance (` ```math m_0 ```) if it is not known directly. The script uses a more comprehensive physical model, considering wave shoaling and breaking (based on Thornton & Guza, 1983) to derive ` ```math m_0 ``` from offshore wave parameters like peak period and wave height. It serves as a useful pre-processing tool to generate the necessary ` ```math m_0 ``` input for the main C++ applications.

## References

* **Battjes, J. A., & Groenendijk, H. W. (2000).** Wave height distributions on shallow foreshores. *Coastal Engineering*, 40(3), 161-182.
* **Caires, S., & Van Gent, M. R. A. (2012).** Wave height distribution in constant and finite depths. *Coastal Engineering Proceedings*, 1(33), 15.
* **Goda, Y. (1979).** A review on statistical interpretation of wave data. *Report of the Port and Harbour Research Institute, Japan*, 18(1), 5-32.
* **Groenendijk, H. W. (1998).** *Shallow foreshore wave height statistics*. M.Sc.-thesis, Delft University of Technology, Department of Civil Engineering, Section Fluid Mechanics, The Netherlands.
* **Groenendijk, H. W., & Van Gent, M. R. A. (1998).** *Shallow foreshore wave height statistics; A predictive model for the probability of exceedance of wave heights*. Technical Report H3351, WL | delft hydraulics, The Netherlands.
* **Hasselmann, K., Barnett, T. P., Bouws, E., Carlson, H., Cartwright, D. D. E., Enke, K.,... & Walden, H. (1973).** *Measurements of Wind-Wave Growth and Swell Decay during the Joint North Sea Wave Project (JONSWAP)*. Ergnzungsheft zur Deutschen Hydrographischen Zeitschrift Reihe, A (8), 95.
* **Karmpadakis, I., Swan, C., & Christou, M. (2020).** Assessment of wave height distributions using an extensive field database. *Coastal Engineering*, 157, 103630.
* **Karmpadakis, I., Swan, C., & Christou, M. (2022).** A new wave height distribution for intermediate and shallow water depths. *Coastal Engineering*, 175, 104130.
* **Klopman, G. (1996).** *Extreme wave heights in shallow water*. WL | delft hydraulics, Report H2486, The Netherlands.
* **Klopman, G., & Stive, M. J. F. (1989).** *Extreme waves and wave loading in shallow water*. Paper presented at the E&P Forum Workshop in Paris, Delft Hydraulics, The Netherlands.
* **Longuet-Higgins, M. S. (1952).** On the statistical distribution of heights of sea waves. *Journal of Marine Research*, 11(3), 245-266.
* **Longuet-Higgins, M. S. (1980).** On the distribution of the heights of sea waves: Some effects of nonlinearity and finite band width. *Journal of Geophysical Research*, 85(C3), 1519-1523.
* **Naess, A. (1985).** On the distribution of crest to trough wave heights. *Ocean Engineering*, 12(3), 221-234.
* **Rice, S. O. (1944).** Mathematical analysis of random noise. *Bell System Technical Journal*, 23(3), 282-332.
* **Tayfun, M. A. (1990).** Distribution of large wave heights. *Journal of Waterway, Port, Coastal, and Ocean Engineering*, 116(6), 686-707.
* **Thornton, E. B., & Guza, R. T. (1982).** Energy saturation and phase speeds measured on a natural beach. *Journal of Geophysical Research*, 87(C12), 9499-9508.
* **Thornton, E. B., & Guza, R. T. (1983).** Transformation of wave height distribution. *Journal of Geophysical Research*, 88(C10), 5925-5938.
