# Shallow Water Waves calculator

## Description:

This program computes local shallow-foreshore wave-height distribution parameters using a model based on the Composed Weibull distribution as described in:

> "Shallow foreshore wave height statistics"  
> by H. Groenendijk, Master's Thesis, Delft University of Technology, 1998.

The program is implemented as a Windows GUI application using only the native Win32 API (without any external libraries). It allows the user to input two key parameters:

1. **Hm0** (in meters) - The local significant spectral wave height.
2. **d** (in meters) - The local water depth.

Based on these inputs, the program performs the following steps:

### 1. Input Acquisition:
- The user enters the values for **Hm0** and **d** into two separate edit controls on the main window.

### 2. Parameter Calculation:
  Free-surface variance (**m0**) is computed as: m0 = (Hm0 / 4)²
  The mean square wave height (Hrms) is computed as: Hrms = 3 * sqrt(m0)
  A dimensional transitional wave height (Htr) is then calculated using the corrected formula: Htr = (0.12 * d / sqrt(m0)) * Hrms
  A dimensionless transitional parameter (H̃_tr) is then derived: H̃_tr = Htr / Hrms

This parameter is used as the interpolation point and must lie within the table range (2.3–3.5) for proper interpolation.

### 3. Interpolation of Wave-Height Ratios:

  A predefined 25-row table is used:
  tableX: Values ranging from 2.3 to 3.5 representing H̃_tr.
  col1 to col7: The table columns contain characteristic ratios relative to Hrms (i.e. Hᵢ/Hrms).
  
  Natural cubic spline interpolation is performed for each column at the value H̃_tr, thereby obtaining the dimensionless ratios Hᵢ/Hrms.

### 4. Conversion to Dimensional Quantities:

  The dimensional wave heights (in meters) are then calculated using: H = (Hᵢ/Hrms) * Hrms
  where (Hᵢ/Hrms) are the interpolated values from the table.

### 5. Report Generation:

  A detailed report is generated which includes:
  The input parameters (Hm0 and d).
  The computed intermediate values (m0, Hrms, Htr, H̃_tr).
  The dimensionless wave heights (H/Hrms) as directly interpolated.
  The dimensional wave heights (in meters) computed from Hrms.
  Diagnostic ratios computed from the wave-height values.

### Compilation Instructions:

  To compile this application using g++ on a Windows system, use a command similar to the following:

```sh
g++ -O2 -Wall -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui \
    -mwindows -static -static-libgcc -static-libstdc++
```

Explanation of the flags:

  -O2: Enables level 2 optimization for improved performance.
  -Wall: Enables all compiler's warning messages to help with debugging.
  -municode: Ensures Unicode support.
  -mwindows: Links against the Windows subsystem rather than the console.
  -static: Links statically to reduce dependency on DLLs.
  -static-libgcc: Links the GCC runtime library statically.
  -static-libstdc++: Links the standard C++ library statically.

### References:

H. Groenendijk, "Shallow foreshore wave height statistics", Master's Thesis, Delft University of Technology, 1998.
