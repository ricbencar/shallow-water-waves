# Shallow Water Waves Model (CLI & GUI)

## Overview

This project provides both a **Command-Line Interface (CLI)** and a **Graphical User Interface (GUI)** application for computing local shallow-foreshore wave-height distribution parameters. The underlying model is based on the **Composed Weibull distribution**, as detailed in "Shallow foreshore wave height statistics" by H. Groenendijk (Master's Thesis, Delft University of Technology, 1998).

The applications allow users to input key wave and bathymetric parameters, perform complex wave height distribution calculations, and generate a detailed report of the results.

## Features

* **Dual Interface**: Offers both a command-line interface for quick calculations and scripting, and a graphical user interface for ease of use.

* **Composed Weibull Distribution Model**: Implements a robust model for shallow-foreshore wave-height distribution with empirically determined exponents ($k\_1=2.0$, $k\_2=3.5$).

* **Key Parameter Calculation**: Computes free-surface variance ($m\_0$), mean square wave height ($H\_{rms}$), and dimensional/dimensionless transitional wave heights ($H\_{tr\_dim}$, $\tilde{H}\_{tr}$).

* **Dimensionless Wave-Height Ratios**: Calculates $\tilde{H}\_N$ (wave height with $1/N$ exceedance probability) and $\tilde{H}\_{1/N}$ (mean of the highest $1/N$-part of wave heights) for various characteristic wave heights ($H\_1, H\_2, H\_{1/3}, H\_{1/10}, H\_{1/50}, H\_{1/100}, H\_{1/1000}$). These are determined by solving a system of two non-linear equations derived from the Composite Weibull distribution, ensuring the normalized $H\_{rms}$ of the distribution equals one. The solution employs a robust numerical strategy using a Newton-Raphson matrix method for simultaneous root-finding, and precise implementations of incomplete gamma functions.

* **Dimensional Wave Heights**: Converts dimensionless ratios to actual wave heights in meters.

* **Diagnostic Ratios**: Provides insights into the wave height distribution shape through various diagnostic ratios.

* **Detailed Reporting**: Generates a comprehensive report of input parameters, calculated values, and diagnostic ratios.

## Input Parameters

Both the CLI and GUI applications require the following three input parameters:

* $H\_{m0}$ (in meters): The local significant spectral wave height.

* $d$ (in meters): The local water depth.

* Beach slope ($1:m$): The beach slope expressed as "1:m". For example, entering 20 signifies a slope of $1/20=0.05$.

## Computational Process

Based on the provided inputs, the program performs a series of calculations to determine various wave height parameters and distribution characteristics:

### 1. Free-Surface Variance ($m\_0$)

The free-surface variance, denoted as $m\_0$, is a measure of the total wave energy. It is calculated from the local significant spectral wave height ($H\_{m0}$) using the following formula:

```
m_0 = \left(\frac{H_{m0}}{4}\right)^2
```

### 2. Mean Square Wave Height ($H\_{rms}$)

The mean square wave height ($H\_{rms}$) is an important characteristic of the wave field. Its calculation incorporates empirical coefficients to better capture the shallow-water distribution of extreme waves, deviating from the original deep-water formulas. This formula is derived from Groenendijk (1998), with slightly updated empirical coefficients:

```
H_{rms} = \left(2.69 + 3.24 \cdot \sqrt{\frac{m_0}{d}}\right) \cdot \sqrt{m_0}
```

### 3. Dimensional Transitional Wave Height ($H\_{tr\_dim}$)

The dimensional transitional wave height ($H\_{tr\_dim}$) marks the point where the wave height distribution significantly changes due to depth-induced breaking. It is calculated using the local water depth ($d$) and the beach slope ($m$):
The tangent of the beach slope ($\tan(\alpha)$) is derived from the input $m$:

```
\tan(\alpha) = \frac{1}{m}
```

Then, $H\_{tr\_dim}$ is computed as:

```
H_{tr\_dim} = (0.35 + 5.8 \cdot \tan(\alpha)) \cdot d
```

For example, if $m=20$, then $\tan(\alpha)=1/20=0.05$, and $H\_{tr\_dim}=(0.35+5.8 \cdot 0.05) \cdot d=0.64 \cdot d$.

### 4. Dimensionless Transitional Parameter ($\tilde{H}\_{tr}$)

The dimensionless transitional parameter ($\tilde{H}\_{tr}$) normalizes the dimensional transitional wave height by the mean square wave height:

```
\tilde{H}_{tr} = \frac{H_{tr\_dim}}{H_{rms}}
```

A critical adjustment is applied: if $\tilde{H}\_{tr}$ exceeds 3.5, it is capped at 3.5, and $H\_{tr\_dim}$ is recalculated to maintain consistency, as the model's empirical basis is limited beyond this value:

```
\text{If } \tilde{H}_{tr} > 3.5, \text{ then } \tilde{H}_{tr} = 3.5 \text{ and } H_{tr\_dim} = 3.5 \cdot H_{rms}
```

### 5. Dimensionless Wave-Height Ratios ($\tilde{H}\_N$ and $\tilde{H}\_{1/N}$)

The dimensionless wave-height ratios are critical outputs of the model. The calculation involves solving a system of two non-linear equations derived from the Composite Weibull distribution, ensuring that the normalized $H\_{rms}$ of the distribution equals one. This is achieved using a Newton-Raphson matrix method for simultaneous root-finding.

The core of this calculation is finding the values of $\tilde{H}_1$ and $\tilde{H}_2$ that satisfy the normalized $H_{rms}$ equation (Equation 7.11 from Groenendijk, 1998) and the continuity condition between the two Weibull distributions (Equation 3.4).

Once $\tilde{H}_1$ (the normalized scale parameter of the first Weibull distribution) and $\tilde{H}_2$ (the normalized scale parameter of the second Weibull distribution) are determined, the model calculates two types of dimensionless wave heights:

* **$\tilde{H}\_N$ (Wave Height with $1/N$ Exceedance Probability):** This is the wave height ($H$) such that the probability of a wave exceeding it is $1/N$. It is calculated by first determining a candidate $\tilde{H}\_N$ from the first part of the distribution. If this candidate is less than $\tilde{H}\_{tr}$, then $\tilde{H}\_N$ is taken from the first part. Otherwise, it is taken from the second part of the distribution.

    * If $\tilde{H}\_{N,candidate} < \tilde{H}\_{tr}$: $\tilde{H}\_N = \tilde{H}\_1 \cdot (\ln(N))^{1/k\_1}$

    * If $\tilde{H}\_{N,candidate} \ge \tilde{H}\_{tr}$: $\tilde{H}\_N = \tilde{H}\_2 \cdot (\ln(N))^{1/k\_2}$

* **$\tilde{H}\_{1/N}$ (Mean of the Highest $1/N$-part of Wave Heights):** This represents the average height of the highest $N$-th fraction of waves (e.g., $H\_{1/3}$ for significant wave height). The calculation depends on whether $\tilde{H}\_N$ (from the previous step) falls within the first or second part of the Composite Weibull distribution.

    * **Case 1:** $\tilde{H}\_N < \tilde{H}\_{tr}$ (The wave height with $1/N$ exceedance probability is smaller than the transitional wave height). This scenario implies that the integration for $\tilde{H}\_{1/N}$ spans both parts of the Composite Weibull distribution. The formula used is (Groenendijk 1998, Equation A.10):

        $$\tilde{H}_{1/N} = N \cdot \tilde{H}_1 \cdot \left[\Gamma\left(\frac{1}{k_1}+1, \ln(N)\right) - \Gamma\left(\frac{1}{k_1}+1, \left(\frac{\tilde{H}_{\text{tr}}}{\tilde{H}_1}\right)^{k_1}\right)\right] + N \cdot \tilde{H}_2 \cdot \Gamma\left(\frac{1}{k_2}+1, \left(\frac{\tilde{H}_{\text{tr}}}{\tilde{H}_2}\right)^{k_2}\right)$$

        where Γ(a,x) is the unnormalized upper incomplete gamma function.

    * **Case 2:** $\tilde{H}\_N \ge \tilde{H}\_{tr}$ (The wave height with $1/N$ exceedance probability is greater than or equal to the transitional wave height). In this case, the integration for $\tilde{H}\_{1/N}$ only involves the second part of the Composite Weibull distribution. The formula used is (Groenendijk 1998, Equation A.17):

        $$\tilde{H}_{1/N} = N \cdot \tilde{H}_2 \cdot \Gamma\left(\frac{1}{k_2}+1, \ln(N)\right)$$

### 6. Dimensional Wave Heights ($H$)

The calculated dimensionless wave-height ratios ($\tilde{H}\_N$ or $\tilde{H}\_{1/N}$) are then converted back to dimensional wave heights (in meters) by multiplying them by the mean square wave height ($H\_{rms}$):

```
H = \tilde{H} \cdot H_{rms}
```

### 7. Diagnostic Ratios

Finally, the program computes several diagnostic ratios, which provide insights into the shape of the wave height distribution and the relative significance of extreme waves. These include ratios characteristic wave height ratios $(H\_{1/n})/(H\_{1/3})$ with $n = 10, 50, 100, 250, \text{ and } 1000$.

## Supporting Mathematical Functions

The core calculations rely on precise implementations of fundamental mathematical functions:

* **Complete Gamma Function ($Γ(z)$):** This is a generalization of the factorial function to real and complex numbers. In the implementation, `std::tgamma` is used. For calculating the logarithm of the complete gamma function ($ln(Γ(a)$), `std::lgamma` is employed for improved numerical stability, especially for large values of $a$.

* **Normalized Lower Incomplete Gamma Function ($P(a,x)$ or $Γ(a,x)/Γ(a)$):** This function is computed using a hybrid numerical approach for stability and accuracy. For small values of $x$ (specifically, $x < a + 1.0$), a series expansion is used. For larger values of $x$, a continued fraction expansion is employed. This adaptive strategy ensures robust and precise computation across different input ranges.

* **Normalized Upper Incomplete Gamma Function ($Q(a,x)$ or $Γ(a,x)/Γ(a)$):** This is directly derived as $1 - P(a,x)$.

* **Unnormalized Upper Incomplete Gamma Function ($Γ(a,x)$):** This is calculated as $Q(a,x) \cdot Γ(a)$.

## Building and Running

### Command-Line Interface (CLI)

The CLI application (`shallow-water-waves_cli.cpp`) can be compiled using g++ on Windows (or similar compilers on other systems).

**Compilation Instructions (Windows example with g++):**

```bash
g++ -O2 -Wall shallow-water-waves_cli.cpp -o shallow-water-waves_cli -static -static-libgcc -static-libstdc++
```

**Usage:**
You can run the CLI application by providing the parameters as command-line arguments or by entering them interactively.

* **With command-line arguments (e.g., Hm0=2.5, d=10, slopeM=20):**

    ```
    shallow-water-waves_cli 2.5 10 20
    ```

* **Interactive input:**

    ```
    shallow-water-waves_cli
    ```

    (The program will then prompt you for `Hm0`, `d`, and `beach slope m`.)

### Graphical User Interface (GUI)

The GUI application (`shallow-water-waves_gui.cpp`) is implemented using the native Win32 API and standard C++. It can be compiled using g++ on Windows.

**Compilation Instructions (Windows example with g++ and OpenMP):**

```
g++ -O3 -Wall -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui -mwindows -static -static-libgcc -static-libstdc++ -fopenmp
```

**Usage:**
Run the compiled executable (`shallow-water-waves_gui.exe`). A window will appear where you can input the `Hm0`, `d`, and `Beach slope m` values in the respective text fields and click "Compute" to see the results. The report will also be saved to `report.txt` in the same directory as the executable.

## Credits

This project's underlying model is based on:

* **Groenendijk, H. (1998).** *Shallow foreshore wave height statistics.* Master's Thesis, Delft University of Technology. [[https://repository.tudelft.nl/record/uuid:fe03dda9-40d9-4046-87fb-459f01fcd3d3](https://repository.tudelft.nl/record/uuid:fe03dda9-40d9-4046-87fb-459f01fcd3d3)]
