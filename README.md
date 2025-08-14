# Shallow Water Waves Model (CLI & GUI)

## Overview

This project presents a solution for accurately computing local shallow-foreshore wave-height distribution parameters, which are critical for coastal engineering applications. It offers both a **Command-Line Interface (CLI)** and a **Graphical User Interface (GUI)**, providing flexibility for users to interact with the model. The core of this project lies in its implementation of the **Composite Weibull distribution**, a statistical model specifically designed for shallow water environments.

Unlike deep-water conditions where wave behavior is often well-described by simpler distributions like Rayleigh, shallow foreshores introduce complex phenomena such as depth-induced breaking, non-linear wave-wave interactions, and bottom friction. These processes significantly alter wave height distributions, making accurate prediction challenging (Groenendijk, 1998; Battjes & Groenendijk, 2000). The Composite Weibull distribution addresses these complexities by employing a two-part structure that effectively captures the distinct physical regimes of unbroken and breaking waves.

Through these applications, users can input essential wave parameters to perform local wave heights distribution calculations. The output includes a detailed report of the computed values, offering insights for the design of marine structures, coastal protection systems, and the assessment of wave-related phenomena like overtopping.
<img width="1083" height="463" alt="shallow-water-waves" src="https://github.com/user-attachments/assets/839a6c12-46d9-4470-b1f9-92394d305c43" />

## Composite Weibull Distribution

The Composite Weibull distribution is a two-part distribution specifically designed for shallow water environments, capturing the distinct physical regimes of unbroken and breaking waves. It is described by the following cumulative distribution function (CDF) for normalized wave heights ($\tilde{h} = h / H_{rms}$):

```math
\Large F(\tilde{h}) = \begin{cases} 
1 - \exp\left[ - \left( \frac{\tilde{h}}{\tilde{H}_1} \right)^{k_1} \right], & \tilde{h} < \tilde{H}_{tr} \\
1 - \exp\left[ - \left( \frac{\tilde{h}}{\tilde{H}_2} \right)^{k_2} \right], & \tilde{h} \geq \tilde{H}_{tr} 
\end{cases}
```

Where:

* $\tilde{h} = h / H\_{rms}$ is the normalized wave height.

* $\tilde{H}_1 = H\_1 / H\_{rms}$ is the normalized scale parameter for the first part of the distribution (unbroken waves).

* $\tilde{H}_2 = H\_2 / H\_{rms}$ is the normalized scale parameter for the second part of the distribution (breaking waves).

* $\tilde{H}_{tr} = H\_{tr} / H\_{rms}$ is the dimensionless transitional wave height, marking the boundary between the two parts of the distribution.

* $k\_1 = 2.0$ is the exponent (shape parameter) for the first part, which is Rayleigh-shaped.

* $k\_2 = 3.6$ is the empirically determined exponent for the second part.

The parameters $\tilde{H}_1$ and $\tilde{H}_2$ are determined by solving a system of non-linear equations to ensure consistency with the normalized $H\_{rms}$ and continuity at the transitional wave height. These equations are:

1.  **Normalized** $H\_{rms}$ **Constraint (from Groenendijk, 1998, Equation 7.11):**
   ```math
   \Large \tilde{H}_{rms} = 1 = \sqrt{
   \tilde{H}_1^{2} \, \gamma\left( \frac{2}{k_1} + 1, \left( \frac{\tilde{H}_{tr}}{\tilde{H}_1} \right)^{k_1} \right)
   + \tilde{H}_2^{2} \, \Gamma\left( \frac{2}{k_2} + 1, \left( \frac{\tilde{H}_{tr}}{\tilde{H}_2} \right)^{k_2} \right)
}
   ```
This equation ensures that the overall root-mean-square of the normalized composite distribution equals one.

2.  **Continuity Condition (from Groenendijk, 1998, Equation 3.4):**
   ```math
   \Large \left( \frac{\tilde{H}_{tr}}{\tilde{H}_1} \right)^{k_1} = \left( \frac{\tilde{H}_{tr}}{\tilde{H}_2} \right)^{k_2}
   ```
This condition ensures that the cumulative distribution function is continuous at the transitional wave height $\tilde{H}_{tr}$.

## Features

* **Dual Interface**: Offers both a command-line interface for quick calculations and scripting, and a graphical user interface for ease of use.

* **Composite Weibull Distribution Model**: Implements a robust model for shallow-foreshore wave-height distribution with empirically determined exponents ($k\_1=2.0$, $k\_2=3.6$). This two-part distribution is designed to reflect the different physical regimes governing smaller (unbroken) and larger (breaking) waves.

* **Key Parameter Calculation**: Computes free-surface variance ($m\_0$), mean square wave height ($H\_{rms}$), and dimensional/dimensionless transitional wave heights ($H\_{tr\_dim}$, $\tilde{H}\_{tr}$).

* **Dimensionless Wave-Height Ratios**: Calculates $\tilde{H}\_N$ (wave height with $1/N$ exceedance probability) and $\tilde{H}\_{1/N}$ (mean of the highest $1/N$-part of wave heights) for various characteristic wave heights ($H\_1, H\_2, H\_{1/3}, H\_{1/10}, H\_{1/50}, H\_{1/100}, H\_{1/1000}$). These are determined by solving a system of two non-linear equations derived from the Composite Weibull distribution, ensuring the normalized $H\_{rms}$ of the distribution equals one. The solution employs a robust numerical strategy using a Newton-Raphson matrix method for simultaneous root-finding, and precise implementations of incomplete gamma functions.

* **Dimensional Wave Heights**: Converts dimensionless ratios to actual wave heights in meters.

* **Diagnostic Ratios**: Provides insights into the wave height distribution shape through various diagnostic ratios.

* **Detailed Reporting**: Generates a comprehensive report of input parameters, calculated values, and diagnostic ratios.

## Input Parameters

Both the CLI and GUI applications require the following three input parameters:

* $H_{m0}$ (in meters): The local significant spectral wave height.
* $d$ (in meters): The local water depth.
* Beach slope ($m$): The beach slope expressed as "1:m". For example, entering `100` signifies a slope of $1/100=0.01$.

## Computational Process

Based on the provided inputs, the program performs a series of calculations to determine various wave height parameters and distribution characteristics. The process has been updated to derive the free-surface variance ($m_0$) internally.

### 1. Mean Square Wave Height ($H_{rms}$)

The process begins by calculating the mean square wave height ($H_{rms}$) from the input significant wave height ($H_{m0}$) and water depth ($d$) using the **van Vledder (1991)** distribution. This provides a more physically grounded starting point for shallow water conditions. The formula is:

```math
\Large H_{rms} = H_{m0} \sqrt{ \frac{\Gamma\left(\frac{2}{k} + 1\right)}{\Gamma\left(\frac{1}{k} + 1\right)} }
```

where the parameter $k$ is dependent on the depth-induced saturation:

```math
k = \frac{2}{1 - (H_{m0} / d)}
```

### 2. Free-Surface Variance ($m_0$)

Next, the free-surface variance ($m_0$, the zeroth spectral moment) is back-calculated from the newly computed $H_{rms}$ and the water depth $d$. This is done by numerically solving the following equation for $\sqrt{m_0}$:

```math
\Large H_{rms} = \left(2.69 + \frac{3.24}{d} \sqrt{m_0}\right) \sqrt{m_0}
```

This relationship, which demonstrates an increase with the degree of saturation ($\Psi = \sqrt{m_0}/d$), is inverted to find $m_0$.

### 3. Dimensional Transitional Wave Height ($H_{tr}$)

The dimensional transitional wave height ($H_{tr}$) marks the point where the wave height distribution significantly changes due to depth-induced breaking. It is calculated using the local water depth ($d$) and the beach slope ($m$):
The tangent of the beach slope ($\tan{a}$) is derived from the input $m$:

```math
\Large \tan(\alpha) = \frac{1}{m}
```

Then, $H_{tr}$ is computed as:

```math
\Large H_{tr} = (0.35 + 5.8 \cdot \tan(\alpha)) \cdot d
```

For example, if $m=100$, then $\tan(\alpha)=1/100=0.01$, and $H_{tr}=(0.35+5.8 \cdot 0.01) \cdot d=0.408 \cdot d$. This relationship indicates that steeper slopes tend to result in higher $H_{tr}$ values, implying that fewer waves deviate from the Rayleigh distribution on steeper foreshores.

### 4. Dimensionless Transitional Parameter ($\tilde{H}_{tr}$)

The dimensionless transitional parameter ($\tilde{H}_{tr}$) normalizes the dimensional transitional wave height by the mean square wave height:

```math
\Large \tilde{H}_{tr} = \frac{H_{tr}}{H_{rms}}
```

### 5. Dimensionless Wave-Height Ratios ($\tilde{H}\_N$ and $\tilde{H}\_{1/N}$)

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

```math
\Large \tilde{H}_{1/N} \left[ \Gamma\left(\frac{1}{k_1}+1, \ln(N)\right) - \Gamma\left(\frac{1}{k_1}+1, \left(\frac{H_{tr}}{H_1}\right)^{k_1}\right) \right] + NH_2 \Gamma\left(\frac{1}{k_2}+1, \left(\frac{H_{tr}}{H_2}\right)^{k_2}\right)
```

where $Γ(a,x)$ is the unnormalized upper incomplete gamma function.

    * **Case 2:** $\tilde{H}\_N \ge \tilde{H}\_{tr}$ (The wave height with $1/N$ exceedance probability is greater than or equal to the transitional wave height). In this case, the integration for $\tilde{H}\_{1/N}$ only involves the second part of the Composite Weibull distribution. The formula used is (Groenendijk 1998, Equation A.17):

```math
\Large \tilde{H}_{1/N} = N \cdot \tilde{H}_2 \cdot \Gamma\left(\frac{1}{k_2}+1, \ln(N)\right)
```

### 6. Dimensional Wave Heights ($H$)

The calculated dimensionless wave-height ratios ($\tilde{H}\_N$ or $\tilde{H}\_{1/N}$) are then converted back to dimensional wave heights (in meters) by multiplying them by the mean square wave height ($H\_{rms}$):

```math
\Large H = \tilde{H} \cdot H_{rms}
```

### 7. Diagnostic Ratios

Finally, the program computes several diagnostic ratios, which provide insights into the shape of the wave height distribution and the relative significance of extreme waves. These include characteristic wave height ratios $(H_{1/N})/(H_{1/3})$ with $N = 10, 50, 100, 250, \text{ and } 1000$.

## Supporting Mathematical Functions

The core calculations rely on precise implementations of fundamental mathematical functions:

* **Complete Gamma Function (Γ(z)):** This is a generalization of the factorial function to real and complex numbers. In the implementation, `std::tgamma` is used. For calculating the logarithm of the complete gamma function (ln(Γ(a))), `std::lgamma` is employed for improved numerical stability, especially for large values of $a$.

```math
\Large Γ(a) = \int_0^{\infty} t^{a-1} e^{-t} dt \quad (a > 0)
```

* **Unnormalized Lower Incomplete Gamma Function (γ(a,x)):** This function is computed using a hybrid numerical approach for stability and accuracy. For small values of $x$ (specifically, $x < a + 1.0$), a series expansion is used. For larger values of $x$, a continued fraction expansion is employed. This adaptive strategy ensures robust and precise computation across different input ranges.

```math
\Large γ(a, x) = \int_0^x t^{a-1} e^{-t} dt
```

* **Unnormalized Upper Incomplete Gamma Function (Γ(a,x)):** This is calculated as Γ(a) - Γ(a,x).

```math
\Large Γ(a, x) = \int_x^{\infty} t^{a-1} e^{-t} dt
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

* **With command-line arguments (e.g., Hm0=2.5, d=5.0, slopeM=100):**
    ```bash
    shallow-water-waves_cli 2.5 5.0 100
    ```

* **Interactive input:**
    ```bash
    shallow-water-waves_cli
    ```
    (The program will then prompt you for `Hm0`, `d`, and `beach slope m`.)

### Graphical User Interface (GUI)

The GUI application (`shallow-water-waves_gui.cpp`) is implemented using the native Win32 API and standard C++. It can be compiled using g++ on Windows.

**Compilation Instructions (Windows example with g++):**

```bash
g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui -mwindows -static -static-libgcc -static-libstdc++
```

**Usage:**
Run the compiled executable (`shallow-water-waves_gui.exe`). A window will appear where you can input the `Hm0`, `d`, and `Beach slope m` values in the respective text fields and click "Compute" to see the results. The report will also be saved to text file `report.txt` in the same directory as the executable.

## References

* **Battjes, J. A., & Groenendijk, H. W. (2000).** Wave height distributions on shallow foreshores. *Coastal Engineering*, 40(3), 161-182. [https://data-ww3.ifremer.fr/BIB/Battjes_Groenendijk_CE2000.pdf](https://data-ww3.ifremer.fr/BIB/Battjes_Groenendijk_CE2000.pdf)
* **Caires, S., & Van Gent, M. R. A. (2012).** Wave height distribution in constant and finite depths. *Coastal Engineering Proceedings*, 1(33), 15. [https://doi.org/10.9753/icce.v33.waves.15](https://doi.org/10.9753/icce.v33.waves.15)
* **Goda, Y. (1979).** A review on statistical interpretation of wave data. *Report of the Port and Harbour Research Institute, Japan*, 18(1), 5-32. [https://www.pari.go.jp/PDF/vol018-no01-01.pdf](https://www.pari.go.jp/PDF/vol018-no01-01.pdf)
* **Groenendijk, H. W. (1998).** *Shallow foreshore wave height statistics*. M.Sc.-thesis, Delft University of Technology, Department of Civil Engineering, Section Fluid Mechanics, The Netherlands. [http://resolver.tudelft.nl/uuid:fe03dda9-40d9-4046-87fb-459f01fcd3d3](http://resolver.tudelft.nl/uuid:fe03dda9-40d9-4046-87fb-459f01fcd3d3)
* **Groenendijk, H. W., & Van Gent, M. R. A. (1998).** *Shallow foreshore wave height statistics; A predictive model for the probability of exceedance of wave heights*. Technical Report H3351, WL | delft hydraulics, The Netherlands.
* **Hasselmann, K., Barnett, T. P., Bouws, E., Carlson, H., Cartwright, D. D. E., Enke, K.,... & Walden, H. (1973).** *Measurements of Wind-Wave Growth and Swell Decay during the Joint North Sea Wave Project (JONSWAP)*. Ergnzungsheft zur Deutschen Hydrographischen Zeitschrift Reihe, A (8), 95.
* **James, J. P., & Panchang, V. (2022).** Investigation of Wave Height Distributions and Characteristic Wave Periods in Coastal Environments. *Journal of Geophysical Research: Oceans*, 127, e2021JC018144. [https://doi.org/10.1029/2021JC018144](https://doi.org/10.1029/2021JC018144)
* **Karmpadakis, I., Swan, C., & Christou, M. (2020).** Assessment of wave height distributions using an extensive field database. *Coastal Engineering*, 157, 103630. [https://doi.org/10.1016/j.coastaleng.2019.103630](https://doi.org/10.1016/j.coastaleng.2019.103630)
* **Karmpadakis, I., Swan, C., & Christou, M. (2022).** A new wave height distribution for intermediate and shallow water depths. *Coastal Engineering*, 175, 104130. [https://doi.org/10.1016/j.coastaleng.2022.104130](https://doi.org/10.1016/j.coastaleng.2022.104130)
* **Klopman, G. (1996).** *Extreme wave heights in shallow water*. WL | delft hydraulics, Report H2486, The Netherlands.
* **Klopman, G., & Stive, M. J. F. (1989).** *Extreme waves and wave loading in shallow water*. Paper presented at the E&P Forum Workshop in Paris, Delft Hydraulics, The Netherlands. [https://repository.tudelft.nl/islandora/object/uuid:3d33138b-94d9-46a8-be8c-1ffb39c64a7f/datastream/OBJ/download](https://repository.tudelft.nl/islandora/object/uuid:3d33138b-94d9-46a8-be8c-1ffb39c64a7f/datastream/OBJ/download)
* **Longuet-Higgins, M. S. (1952).** On the statistical distribution of heights of sea waves. *Journal of Marine Research*, 11(3), 245-266. [https://elischolar.library.yale.edu/journal_of_marine_research/774/](https://elischolar.library.yale.edu/journal_of_marine_research/774/)
* **Longuet-Higgins, M. S. (1980).** On the distribution of the heights of sea waves: Some effects of nonlinearity and finite band width. *Journal of Geophysical Research*, 85(C3), 1519-1523. [https://doi.org/10.1029/JC085iC03p01519](https://doi.org/10.1029/JC085iC03p01519)
* **Naess, A. (1985).** On the distribution of crest to trough wave heights. *Ocean Engineering*, 12(3), 221-234. [https://doi.org/10.1016/0029-8018(85)90014-9](https://doi.org/10.1016/0029-8018(85)90014-9)
* **Rice, S. O. (1944).** Mathematical analysis of random noise. *Bell System Technical Journal*, 23(3), 282-332. [https://doi.org/10.1002/j.1538-7305.1944.tb00874.x](https://doi.org/10.1002/j.1538-7305.1944.tb00874.x)
* **Tayfun, M. A. (1990).** Distribution of large wave heights. *Journal of Waterway, Port, Coastal, and Ocean Engineering*, 116(6), 686-707. [https://doi.org/10.1061/(ASCE)0733-950X(1990)116:6(686)](https://doi.org/10.1061/(ASCE)0733-950X(1990)116:6(686))
* **Thornton, E. B., & Guza, R. T. (1982).** Energy saturation and phase speeds measured on a natural beach. *Journal of Geophysical Research*, 87(C12), 9499-9508. [https://doi.org/10.1029/JC087iC12p09499](https://doi.org/10.1029/JC087iC12p09499)
* **Thornton, E. B., & Guza, R. T. (1983).** Transformation of wave height distribution. *Journal of Geophysical Research*, 88(C10), 5925-5938. [https://doi.org/10.1029/JC088iC10p05925](https://doi.org/10.1029/JC088iC10p05925)
* **van Vledder, G. P. (1991).** *A practical method for the computation of the significant wave height in shallow water*. In: *Coastal Engineering 1990*, pp. 653-665. ASCE.
