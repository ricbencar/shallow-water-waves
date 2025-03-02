# Shallow Water Waves calculator

# Program: `shallow-water-waves_gui.cpp`

## Detailed Description

This program computes local shallow-foreshore wave-height distribution parameters using a model based on the Composed Weibull distribution as described in:

> "Shallow foreshore wave height statistics"  
> by H. Groenendijk, Master's Thesis, Delft University of Technology, 1998.

The program is implemented as a Windows GUI application using only the native Win32 API (without any external libraries). It allows the user to input two key parameters:

1. **$H_{m0}$** (in meters): The local significant spectral wave height.  
2. **$d$** (in meters): The local water depth.

Based on these inputs, the program performs the following steps:

---

### 1. Input Acquisition

- The user enters the values for $H_{m0}$ and $d$ into two separate edit controls on the main window.

---

### 2. Parameter Calculation

- **Free-surface variance ($m_0$)** is computed as:

```math
m_0 = \left(\frac{H_{m0}}{4}\right)^2

    Mean square wave height ($H_{rms}$) is computed as:

Hrms=3m0
Hrms​=3m0​
​

    Dimensional transitional wave height ($H_{tr}$) is calculated using the corrected formula:

Htr=(0.12 dm0) Hrms
Htr​=(m0​
​0.12d​)Hrms​

    Dimensionless transitional parameter ($\tilde{H}_{tr}$) is derived as:

H~tr=HtrHrms
H~tr​=Hrms​Htr​​

This parameter is used as the interpolation point and must lie within the table range (2.3–3.5) for proper interpolation.
3. Interpolation of Wave-Height Ratios

    A predefined 25-row table is used:
        tableX: Values ranging from 2.3 to 3.5 representing $\tilde{H}_{tr}$.
        col1 to col7: The table columns contain characteristic ratios relative to $H_{rms}$ (i.e., $\frac{H_i}{H_{rms}}$).

    Natural cubic spline interpolation is performed for each column at the value $\tilde{H}{tr}$, thereby obtaining the dimensionless ratios $\frac{H_i}{H{rms}}$.

4. Conversion to Dimensional Quantities

    The dimensional wave heights (in meters) are then calculated using:

Hi=(HiHrms) Hrms
Hi​=(Hrms​Hi​​)Hrms​

where $\frac{H_i}{H_{rms}}$ are the interpolated values from the table.
5. Report Generation

    A detailed report is generated which includes:
        The input parameters ($H_{m0}$ and $d$).
        The computed intermediate values ($m_0$, $H_{rms}$, $H_{tr}$, $\tilde{H}_{tr}$).
        The dimensionless wave heights ($\frac{H_i}{H_{rms}}$) as directly interpolated.
        The dimensional wave heights (in meters) computed from $H_{rms}$.
        Diagnostic ratios computed from the wave-height values.

6. Graphical User Interface (GUI)

    The main window is designed as a non-resizable window.
    It contains:
        Two edit controls for inputting $H_{m0}$ and $d$.
        A "Compute" button that triggers the parameter computations.
        A multiline, read-only output control (using a 20-pt Courier New font) which displays the detailed report.
    A helper routine ensures that Unix-style newline characters (\n) are converted to the Windows CR-LF (\r\n) format for proper display.

## Compilation Instructions

To compile this application using g++ on a Windows system, use a command similar to the following:

g++ -O2 -Wall -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui \
    -mwindows -static -static-libgcc -static-libstdc++

Explanation of the flags:

    -O2: Enables level 2 optimization for improved performance.
    -Wall: Enables all compiler warning messages to help with debugging.
    -municode: Ensures Unicode support.
    -mwindows: Links against the Windows subsystem rather than the console.
    -static: Links statically to reduce dependency on DLLs.
    -static-libgcc: Links the GCC runtime library statically.
    -static-libstdc++: Links the standard C++ library statically.

## References

    H. Groenendijk, "Shallow foreshore wave height statistics", Master's Thesis, Delft University of Technology, 1998.
