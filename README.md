# Shallow Water Waves Calculator

## Description:

This program computes local shallow-foreshore wave-height distribution parameters using a model based on the Composed Weibull distribution as described in:

> H. Groenendijk, "Shallow foreshore wave height statistics", Master's Thesis, Delft University of Technology, 1998. URL: (https://repository.tudelft.nl/record/uuid:fe03dda9-40d9-4046-87fb-459f01fcd3d3)

![shallow-water-waves](https://github.com/user-attachments/assets/31154777-4b6f-4c90-bb2d-13b83aafc7ba)

It is implemented as a Windows GUI application using only the native Win32 API (without any external libraries). It allows the user to input two key parameters:

1. **Hm0** (in meters) - The local significant spectral wave height.
2. **d** (in meters) - The local water depth.

Based on these inputs, the program performs the following steps:

### 1. Parameter Calculation:
- **Free-Surface Variance (m₀):**  
  Computed as:  
  `m₀ = (0.75 * Hm0 / 4)²`
- **Mean Square Wave Height (H<sub>rms</sub>):**  
  Calculated using the formula:  
  `H₍rms₎ = (2.69 + 3.24 * sqrt(m₀) / d) * sqrt(m₀)`
- **Dimensional Transitional Wave Height (H<sub>tr</sub>):**  
  Determined as:  
  `H₍tr₎ = (0.12 * d / sqrt(m₀)) * H₍rms₎`
- **Dimensionless Transitional Parameter (H̃<sub>tr</sub>):**  
  Derived from:  
  `H̃₍tr₎ = H₍tr₎ / H₍rms₎`  
  This parameter is used as the interpolation point.

### 2. Interpolation of Wave-Height Ratios:
- A predefined **70-row table** is used with values ranging from **0.05 to 3.50** (in increments of 0.05) representing H̃<sub>tr</sub>.
- The table contains **7 columns** (col1 to col7) with characteristic dimensionless wave-height ratios relative to H<sub>rms</sub> (i.e. Hᵢ/H<sub>rms</sub>).
- **Natural cubic spline interpolation** is performed on each column at the computed H̃<sub>tr</sub> to obtain the interpolated dimensionless ratios.

### 3. Conversion to Dimensional Quantities:
- The dimensional wave heights (in meters) are then calculated by:
  
  `H = (Hᵢ/H₍rms₎) * H₍rms₎`
  
  where **Hᵢ/H₍rms₎** are the interpolated values from the table.

### 4. Report Generation:
- A detailed report is generated which includes:
  - The input parameters (**Hm0** and **d**).
  - The computed intermediate values (**m₀**, **H<sub>rms</sub>**, **H<sub>tr</sub>**, **H̃<sub>tr</sub>**).
  - The dimensionless wave heights (H/H<sub>rms</sub>) as directly interpolated.
  - The dimensional wave heights (in meters) computed from H<sub>rms</sub>.
  - Diagnostic ratios computed from the wave-height values.

### Compilation Instructions:
To compile this application using g++ on a Windows system, you can use the following command:

```sh
g++ -O2 -Wall -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui \
    -mwindows -static -static-libgcc -static-libstdc++
```
Explanation of the flags:

- -O2: Enables level 2 optimization for improved performance.
- -Wall: Activates all compiler warning messages to help with debugging.
- -municode: Ensures Unicode support.
- -mwindows: Links against the Windows subsystem rather than the console.
- -static, -static-libgcc, -static-libstdc++: Link statically to reduce dependency on DLLs.

### References:

- Battjes, J. A. and Groenendijk, H. W. (2000). "Wave Height Distributions on Shallow Foreshores," Coastal Engineering, 40, 161–182. DOI: https://doi.org/10.1016/S0378-3839(00)00019-2
- Groenendijk, H. (1998). "Shallow foreshore wave height statistics", Master's Thesis, Delft University of Technology. TU Delft Repository
- Goda, Y. (1975, 2010). "Deformation of Irregular Waves due to Depth-Controlled Wave Breaking" and Random Seas and Design of Maritime Structures, World Scientific – Random Seas and Design of Maritime Structures
