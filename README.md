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

The parameters $\tilde{H}_1$ and $\tilde{H}_2$ are determined by solving a system of non-linear equations to ensure consistency with the normalized $H_{rms}$ and continuity at the transitional wave height. These equations are:

1.  **Normalized $H_{rms}$ Constraint (from Groenendijk, 1998, Equation 7.11):**
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
    This condition ensures that the cumulative distribution function is continuous at the transitional wave height $\tilde{H}_{tr}$.

## Features

* **Dual Interface**: Offers both a command-line interface for quick calculations and scripting, and a graphical user interface for ease of use.
* **Composite Weibull Distribution Model**: Implements a robust model for shallow-foreshore wave-height distribution with empirically determined exponents ($k_1=2.0$, $k_2=3.6$). This two-part distribution is designed to reflect the different physical regimes governing smaller (unbroken) and larger (breaking) waves.
* **Key Parameter Calculation**: Takes free-surface variance ($m_0$) as input and computes the mean square wave height ($H_{rms}$), and dimensional/dimensionless transitional wave heights ($H_{tr\_dim}$, $\tilde{H}_{tr}$).
* **Dimensionless Wave-Height Ratios**: Calculates $\tilde{H}_N$ (wave height with $1/N$ exceedance probability) and $\tilde{H}_{1/N}$ (mean of the highest $1/N$-part of wave heights) for various characteristic wave heights ($H_1, H_2, H_{1/3}, H_{1/10}, H_{1/50}, H_{1/100}, H_{1/1000}$). These are determined by solving a system of two non-linear equations derived from the Composite Weibull distribution, ensuring the normalized $H_{rms}$ of the distribution equals one. The solution employs a robust numerical strategy using a Newton-Raphson matrix method for simultaneous root-finding, and precise implementations of incomplete gamma functions.
* **Dimensional Wave Heights**: Converts dimensionless ratios to actual wave heights in meters.
* **Diagnostic Ratios**: Provides insights into the wave height distribution shape through various diagnostic ratios.
* **Detailed Reporting**: Generates a comprehensive report of input parameters, calculated values, and diagnostic ratios.

## Input Parameters

Both the CLI and GUI applications require the following four input parameters:

* $H_{m0}$ (in meters): The local significant spectral wave height. This value is used for reporting and context but is not a direct input to the core calculation, which now starts from $m_0$.
* $m_0$ (in m²): The free-surface variance (zeroth spectral moment). This is now the primary input representing the total wave energy.
* $d$ (in meters): The local water depth.
* Beach slope ($1:m$): The beach slope expressed as "1:m". For example, entering 20 signifies a slope of $1/20=0.05$. This parameter influences the transitional wave height.

## Computational Process

Based on the provided inputs, the program performs a series of calculations to determine various wave height parameters and distribution characteristics:

### 1. Mean Square Wave Height ($H_{rms}$)

The mean square wave height ($H_{rms}$) is an important characteristic of the wave field. Its calculation now starts directly from the input free-surface variance ($m_0$) and incorporates empirical coefficients to better capture the shallow-water distribution of extreme waves. This formula, used in both the CLI and GUI implementations, is:

```math
H_{rms} = \left(2.69 + 3.24 \cdot \frac{\sqrt{m_0}}{d}\right) \cdot \sqrt{m_0}
```

This relationship demonstrates an increase with the degree of saturation ($\Psi = \sqrt{m_0}/d$), counteracting the bandwidth effect observed in deep water.

### 2. Dimensional Transitional Wave Height ($H_{tr\_dim}$)

The dimensional transitional wave height ($H_{tr\_dim}$) marks the point where the wave height distribution significantly changes due to depth-induced breaking. It is calculated using the local water depth ($d$) and the beach slope ($m$):
The tangent of the beach slope ($\tan(\alpha)$) is derived from the input $m$:

```math
\tan(\alpha) = \frac{1}{m}
```

Then, $H_{tr\_dim}$ is computed as:

```math
H_{tr\_dim} = (0.35 + 5.8 \cdot \tan(\alpha)) \cdot d
```

For example, if $m=20$, then $\tan(\alpha)=1/20=0.05$, and $H_{tr\_dim}=(0.35+5.8 \cdot 0.05) \cdot d=0.64 \cdot d$. This relationship indicates that steeper slopes tend to result in higher $H_{tr}$ values, implying that fewer waves deviate from the Rayleigh distribution on steeper foreshores.

### 3. Dimensionless Transitional Parameter ($\tilde{H}_{tr}$)

The dimensionless transitional parameter ($\tilde{H}_{tr}$) normalizes the dimensional transitional wave height by the mean square wave height:

```math
\tilde{H}_{tr} = \frac{H_{tr\_dim}}{H_{rms}}
```

### 4. Dimensionless Wave-Height Ratios ($\tilde{H}\_N$ and $\tilde{H}\_{1/N}$)

The dimensionless wave-height ratios are critical outputs of the model. The calculation involves solving a system of two non-linear equations derived from the Composite Weibull distribution, ensuring that the normalized $H\_{rms}$ of the distribution equals one. This is achieved using a Newton-Raphson matrix method for simultaneous root-finding.

The core of this calculation is finding the values of $\tilde{H}_1$ and $\tilde{H}_2$ that satisfy the normalized $H\_{rms}$ equation (Equation 7.11 from Groenendijk, 1998) and the continuity condition between the two Weibull distributions (Equation 3.4):

```math
f(\tilde{H}_{1_{Hrms}}, \tilde{H}_{2_{Hrms}}, \tilde{H}_{tr}) = \sqrt{
\tilde{H}_{1_{Hrms}}^2 \cdot \gamma\left(\frac{2}{k_1} + 1, \left(\frac{\tilde{H}_{tr}}{\tilde{H}_{1_{Hrms}}}\right)^{k_1}\right) + 
\tilde{H}_{2_{Hrms}}^2 \cdot \Gamma\left(\frac{2}{k_2} + 1, \left(\frac{\tilde{H}_{tr}}{\tilde{H}_{2_{Hrms}}}\right)^{k_2}\right)
} - 1 = 0
```

where $k\_1=2.0$ (representing a Rayleigh-shaped first part of the distribution based on empirical observations for smaller waves) and $k\_2=3.6$ (an empirically determined exponent for the second part, characterizing larger, breaking waves) are global exponents for the Composite Weibull distribution. $H\_{2\_{Hrms}}$ is related to $H\_{1\_{Hrms}}$ and $\tilde{H}\_{tr}$ by the continuity condition between the two Weibull distributions:

```math
H_{2\_{Hrms}} = \tilde{H}_{tr} \cdot \left(\frac{\tilde{H}_{tr}}{H_{1\_{Hrms}}}\right)^{k_1/k_2}
```

Here, $γ(a,x)$ and $Γ(a,x)$ are the unnormalized lower and upper incomplete gamma functions, respectively.

Once $\tilde{H}_1$ (the normalized scale parameter of the first Weibull distribution) and $\tilde{H}_2$ (the normalized scale parameter of the second Weibull distribution) are determined, two types of dimensionless wave heights can be calculated:

* $\tilde{H}\_N$ **(Wave Height with** $1/N$ **Exceedance Probability):** This is the wave height ($H$) such that the probability of a wave exceeding it is $1/N$. It is calculated by first determining a candidate $\tilde{H}\_N$ from the first part of the distribution. If this candidate is less than $\tilde{H}\_{tr}$, then $\tilde{H}\_N$ is taken from the first part. Otherwise, it is taken from the second part of the distribution.

    * If $\tilde{H}\_{N,candidate} < \tilde{H}\_{tr}$: $\tilde{H}\_N = \tilde{H}\_1 \cdot (\ln(N))^{1/k\_1}$

    * If $\tilde{H}\_{N,candidate} \ge \tilde{H}\_{tr}$: $\tilde{H}\_N = \tilde{H}\_2 \cdot (\ln(N))^{1/k\_2}$

* $\tilde{H}\_{1/N}$ **(Mean of the Highest** $1/N$**-part of Wave Heights):** This represents the average height of the highest $N$-th fraction of waves (e.g., $H\_{1/3}$ for significant wave height). The calculation depends on whether $\tilde{H}\_N$ (from the previous step) falls within the first or second part of the Composite Weibull distribution.

    * **Case 1:** $\tilde{H}\_N < \tilde{H}\_{tr}$ (The wave height with $1/N$ exceedance probability is smaller than the transitional wave height). This scenario implies that the integration for $\tilde{H}\_{1/N}$ spans both parts of the Composite Weibull distribution. The formula used is (Groenendijk 1998, Equation A.10):

        $$H_{UN} = NH_1 \left[ \Gamma\left(\frac{1}{k_1}+1, \ln(N)\right) - \Gamma\left(\frac{1}{k_1}+1, \left(\frac{H_{tr}}{H_1}\right)^{k_1}\right) \right] + NH_2 \Gamma\left(\frac{1}{k_2}+1, \left(\frac{H_{tr}}{H_2}\right)^{k_2}\right)$$

        where $Γ(a,x)$ is the unnormalized upper incomplete gamma function.

    * **Case 2:** $\tilde{H}\_N \ge \tilde{H}\_{tr}$ (The wave height with $1/N$ exceedance probability is greater than or equal to the transitional wave height). In this case, the integration for $\tilde{H}\_{1/N}$ only involves the second part of the Composite Weibull distribution. The formula used is (Groenendijk 1998, Equation A.17):

        $$\tilde{H}_{1/N} = N \cdot \tilde{H}_2 \cdot \Gamma\left(\frac{1}{k_2}+1, \ln(N)\right)$$

### 5. Dimensional Wave Heights ($H$)

The calculated dimensionless wave-height ratios ($\tilde{H}_N$ or $\tilde{H}_{1/N}$) are then converted back to dimensional wave heights (in meters) by multiplying them by the mean square wave height ($H_{rms}$):

```math
H = \tilde{H} \cdot H_{rms}
```

### 6. Diagnostic Ratios

Finally, the program computes several diagnostic ratios, which provide insights into the shape of the wave height distribution and the relative significance of extreme waves. These include characteristic wave height ratios $(H_{1/N})/(H_{1/3})$ with $N = 10, 50, 100, 250, \text{ and } 1000$.

## Supporting Mathematical Functions

The core calculations rely on precise implementations of fundamental mathematical functions:

* **Complete Gamma Function (Γ(z)):** This is a generalization of the factorial function to real and complex numbers. In the implementation, `std::tgamma` is used. For calculating the logarithm of the complete gamma function (ln(Γ(a))), `std::lgamma` is employed for improved numerical stability, especially for large values of $a$.

```math
\[
Γ(a) = \int_0^{\infty} t^{a-1} e^{-t} dt \quad (a > 0)
\]
```

* **Unnormalized Lower Incomplete Gamma Function (γ(a,x)):** This function is computed using a hybrid numerical approach for stability and accuracy. For small values of $x$ (specifically, $x < a + 1.0$), a series expansion is used. For larger values of $x$, a continued fraction expansion is employed. This adaptive strategy ensures robust and precise computation across different input ranges.

```math
\[
γ(a, x) = \int_0^x t^{a-1} e^{-t} dt
\]
```

* **Unnormalized Upper Incomplete Gamma Function (Γ(a,x)):** This is calculated as Γ(a) - Γ(a,x).

```math
\[
Γ(a, x) = \int_x^{\infty} t^{a-1} e^{-t} dt
\]
```

## Building and Running

### Command-Line Interface (CLI)

The CLI application (`shallow-water-waves_cli.cpp`) can be compiled using g++ on Windows (or similar compilers on other systems).

**Compilation Instructions (Windows example with g++):**

```bash
g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -o shallow-water-waves_cli shallow-water-waves_cli.cpp
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
    (The program will then prompt you for `Hm0`, `m0`, `d`, and `beach slope m`.)

### Graphical User Interface (GUI)

The GUI application (`shallow-water-waves_gui.cpp`) is implemented using the native Win32 API and standard C++. It can be compiled using g++ on Windows.

**Compilation Instructions (Windows example with g++ and OpenMP):**

```bash
g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui -mwindows -static -static-libgcc -static-libstdc++ -fopenmp
```

**Usage:**
Run the compiled executable (`shallow-water-waves_gui.exe`). A window will appear where you can input the `Hm0`, `m0`, `d`, and `Beach slope m` values in the respective text fields and click "Compute" to see the results. The report will also be saved to `report.txt` in the same directory as the executable.

## Utility Script: m0_calculator.py

For convenience, a Python script `m0_calculator.py` is included. This utility can be used to estimate the free-surface variance ($m_0$) if it is not known directly. The script uses a more comprehensive physical model, considering wave shoaling and breaking (based on Thornton & Guza, 1983) to derive $m_0$ from offshore wave parameters like peak period and wave height. It serves as a useful pre-processing tool to generate the necessary $m_0$ input for the main C++ applications.

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
