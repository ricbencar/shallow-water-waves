# Shallow Water Waves GUI

## Overview

This program computes local shallow-foreshore wave-height distribution parameters using a model based on the Composed Weibull distribution as described in:

- **"Shallow foreshore wave height statistics"** by H. Groenendijk, Master's Thesis, Delft University of Technology, 1998.

The program is implemented as a Windows GUI application using only the native Win32 API and standard C++ (with OpenMP directives for parallelism). It allows the user to input three key parameters:

1. **Hm0 (in meters)** - The local significant spectral wave height.
2. **d (in meters)** - The local water depth.
3. **Beach slope (1:m)** - The beach slope expressed as “1:m”. (For example, enter 20 for a slope of 1/20 = 0.05.)
![shallow-water-waves](https://github.com/user-attachments/assets/31154777-4b6f-4c90-bb2d-13b83aafc7ba)
## Computations

Based on these inputs, the program computes:

- **Free-surface variance**: 
  ```math
  m0 = (Hm0 / 4)²
  ```

- **Mean square wave height**:
  ```math
  Hrms = (3.00 + 3.50 * sqrt(m0) / d) * sqrt(m0)
  ```
  (Empirical coefficients have been slightly increased relative to the original deep-water formula to better capture the shallow-water distribution of extreme waves.)

- **Dimensional transitional wave height**:
  ```math
  Htr = (0.35 + 5.8 * (1/m)) * d.
  ```
  (For example, if m = 20 then tan(alpha) = 1/20 = 0.05 and 
  ```math
  Htr = (0.35 + 5.8 * 0.05) * d = 0.64 * d.
  ```)

- **Dimensionless transitional parameter**:
  ```math
  H̃_tr = Htr / Hrms.
  ```
 If `H̃_tr` is above 3.5, then it is set to 3.5 and `Htr` is recalculated as Htr = 3.5 * Hrms.

Using a 70-row table (with columns for H1/Hrms, H2/Hrms, etc.), a natural cubic spline interpolation is performed at H̃_tr to obtain dimensionless wave-height ratios. These are then converted to dimensional quantities (in meters) by multiplying with Hrms.

A detailed report is then generated (and written to `report.txt`) with the input parameters, intermediate values, interpolated ratios, and computed dimensional wave heights, as well as diagnostic ratios.

## Compilation Instructions

To compile the program using `g++` on Windows with OpenMP, you can use the following command:

```bash
g++ -O3 -Wall -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui \
    -mwindows -static -static-libgcc -static-libstdc++ -fopenmp
```

## References

1. Groenendijk, H.W. (1998). “Shallow foreshore wave height statistics.” M.Sc. thesis, Delft University of Technology (also Delft Hydraulics Rep. H3245).
2. Groenendijk, H.W. & Van Gent, M.R.A. (1998/1999). “Shallow foreshore wave height statistics: predictive model for probability of exceedance of wave heights.” Delft Hydraulics Rep. H3351. [vliz.be](http://vliz.be).
3. Battjes, J.A. & Groenendijk, H.W. (2000). “Wave height distributions on shallow foreshores.” Coastal Engineering, 40(3):161–182. (Composite Rayleigh–Weibull distribution development) [ygraigarw.github.io](http://ygraigarw.github.io).
4. Van Gent, M.R.A. (2001). “Wave runup on dikes with shallow foreshores.” J. Waterway, Port, Coastal, Ocean Eng., 127(5):254–262. (Notes shallow foreshore effects on wave spectra & distribution) [ascelibrary.org](http://ascelibrary.org).
5. Forristall, G.Z. (2007). “Wave height distribution in shallow water.” Ocean Engineering, 34(11–12):1516–1525. (Weibull distribution with Ursell-dependent parameters) [ygraigarw.github.io](http://ygraigarw.github.io).
6. Mai, S. et al. (2011). “Wave height distributions in shallow waters.” Coastal Engineering Proceedings 1(32), paper 57. (Field validation of BG model; recommends coefficient recalibration) [researchgate.net](http://researchgate.net); [pdfs.semanticscholar.org](http://pdfs.semanticscholar.org).
7. Verhagen, H.J. et al. (2008). “A practical method for the design of coastal structures in shallow water.” (Conference paper; emphasizes using Tm-1,0 and H2% for shallow-water design).