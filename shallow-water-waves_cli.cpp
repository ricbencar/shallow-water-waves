/***********************************************************************
 * Program: shallow-water-waves_cli.cpp
 *
 * Detailed Description:
 *   This program computes local shallow-foreshore wave-height distribution
 *   parameters using a model based on the Composed Weibull distribution.
 *
 *   The command-line application performs the following:
 *     1. If two command-line arguments are provided, they are used as
 *        Hm0 (local significant spectral wave height) and d (local water depth).
 *        Otherwise, the program prompts the user for these values.
 *     2. Computes the following intermediate values:
 *            m0      = (Hm0/4)^2          [free-surface variance]
 *            Hrms    = 3 * sqrt(m0)       [mean square wave height]
 *            Htr_dim = (0.12 * d / sqrt(m0)) * Hrms  
 *                      [dimensional transitional wave height]
 *            H̃_tr   = Htr_dim / Hrms     [dimensionless transitional parameter]
 *     3. Interpolates the dimensionless wave-height ratios (Hᵢ/Hrms)
 *        using natural cubic spline interpolation based on a predefined 25-row table.
 *     4. Computes the dimensional wave heights (in meters) as:
 *            H = (Hᵢ/Hrms) * Hrms.
 *     5. Displays a detailed report on the console that exactly matches
 *        the output produced by the GUI version.
 *
 * (Note: The table values, dimensional heights, and diagnostic ratios
 *  are computed using natural cubic spline interpolation on a fixed 25-row
 *  table. They are included here to exactly mirror the GUI version.)
 *
 * Compilation Instructions (for Windows):
 *   g++ -O2 -Wall shallow-water-waves_cli.cpp -o shallow-water-waves_cli \
 *       -static -static-libgcc -static-libstdc++
 *
 *   To run with command-line arguments (e.g., Hm0=1.5 and d=10):
 *       shallow-water-waves_cli 1.5 10
 ***********************************************************************/

#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <locale>
#include <cstdlib>
#include <string>

//---------------------------------------------------------------------
// Function: computeSplineSecondDerivatives
// Purpose:
//   Computes the second derivatives for a set of (x,y) data points using
//   the natural cubic spline method.
// (These functions are included for full functionality as used in the GUI.)
//---------------------------------------------------------------------
std::vector<double> computeSplineSecondDerivatives(const std::vector<double>& x, const std::vector<double>& y) {
    size_t n = x.size();
    std::vector<double> m(n, 0.0);  // Second derivatives.
    std::vector<double> l(n, 0.0), mu(n, 0.0), z(n, 0.0);
    l[0] = 1.0; mu[0] = 0.0; z[0] = 0.0;
    std::vector<double> h(n - 1);
    for (size_t i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }
    std::vector<double> alpha(n, 0.0);
    for (size_t i = 1; i < n - 1; i++) {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i])
                 - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }
    for (size_t i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n - 1] = 1.0; z[n - 1] = 0.0; m[n - 1] = 0.0;
    for (int j = static_cast<int>(n) - 2; j >= 0; j--) {
        m[j] = z[j] - mu[j] * m[j + 1];
    }
    return m;
}

//---------------------------------------------------------------------
// Function: cubicSplineInterpolation
// Purpose:
//   Performs natural cubic spline interpolation for the query point x0.
//---------------------------------------------------------------------
double cubicSplineInterpolation(const std::vector<double>& x, const std::vector<double>& y, double x0) {
    size_t n = x.size();
    if (n == 0)
        throw std::runtime_error("No data points provided in x.");
    if (n == 1)
        return y[0];
    if (x0 <= x.front()) {
        double t = (x0 - x.front()) / (x[1] - x.front());
        return y.front() + t * (y[1] - y.front());
    }
    if (x0 >= x.back()) {
        double t = (x0 - x[n - 2]) / (x.back() - x[n - 2]);
        return y[n - 2] + t * (y.back() - y[n - 2]);
    }
    size_t i = 0;
    while (i < n - 1 && x0 > x[i + 1])
        i++;
    double h = x[i + 1] - x[i];
    double A = (x[i + 1] - x0) / h;
    double B = (x0 - x[i]) / h;
    std::vector<double> m = computeSplineSecondDerivatives(x, y);
    double y0 = A * y[i] + B * y[i + 1] + ((A * A * A - A) * m[i] + (B * B * B - B) * m[i + 1]) * (h * h) / 6.0;
    return y0;
}

//---------------------------------------------------------------------
// Function: buildReport
// Purpose:
//   Computes all wave parameters using the same formulas as in the GUI
//   version, performs table interpolation, and builds a detailed report
//   as a formatted wide string that exactly matches the GUI output.
//---------------------------------------------------------------------
static std::wstring buildReport(double Hm0, double d) {
    std::wstringstream ss;
    ss << std::fixed << std::setprecision(4);
    if (Hm0 <= 0.0 || d <= 0.0) {
        ss << L"ERROR: Both Hm0 and d must be positive.\n";
        return ss.str();
    }
    double m0 = std::pow(Hm0 / 4.0, 2.0);          // Free-surface variance.
    double Hrms = 3.0 * std::sqrt(m0);               // Mean square wave height.
    double Htr_dim = (0.12 * d / std::sqrt(m0)) * Hrms; // Dimensional transitional wave height.
    double Htr_tilde = (Hrms > 0.0) ? (Htr_dim / Hrms) : 0.0; // Dimensionless transitional parameter.

    // Define the 25-row table with H̃_tr values.
    std::vector<double> tableX = { 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75,
                                   2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.1, 3.15, 3.2, 3.25,
                                   3.3, 3.35, 3.4, 3.45, 3.5 };
    // Table columns: Predefined ratios Hᵢ/Hrms.
    std::vector<double> col1 = {0.988, 0.989, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.997,
                                0.998, 0.998, 0.999, 0.999, 0.999, 0.999, 0.999, 1.000, 1.000, 1.000,
                                1.000, 1.000, 1.000, 1.000, 1.000};
    std::vector<double> col2 = {1.419, 1.433, 1.448, 1.462, 1.475, 1.489, 1.502, 1.515, 1.528, 1.540,
                                1.553, 1.565, 1.577, 1.589, 1.601, 1.612, 1.623, 1.635, 1.646, 1.657,
                                1.668, 1.679, 1.689, 1.700, 1.711};
    std::vector<double> col3 = {1.397, 1.400, 1.402, 1.404, 1.406, 1.407, 1.409, 1.410, 1.411, 1.412,
                                1.413, 1.413, 1.414, 1.414, 1.414, 1.415, 1.415, 1.415, 1.415, 1.415,
                                1.415, 1.416, 1.416, 1.416, 1.416};
    std::vector<double> col4 = {1.774, 1.778, 1.781, 1.784, 1.786, 1.789, 1.791, 1.792, 1.794, 1.795,
                                1.796, 1.797, 1.797, 1.798, 1.798, 1.799, 1.799, 1.799, 1.799, 1.799,
                                1.800, 1.800, 1.800, 1.800, 1.800};
    std::vector<double> col5 = {1.951, 1.954, 1.957, 1.960, 1.962, 1.965, 1.967, 1.969, 1.970, 1.971,
                                1.973, 1.974, 1.974, 1.975, 1.976, 1.976, 1.977, 1.977, 1.977, 1.977,
                                1.977, 1.978, 1.978, 1.978, 1.978};
    std::vector<double> col6 = {2.116, 2.120, 2.123, 2.126, 2.129, 2.132, 2.134, 2.136, 2.137, 2.139,
                                2.140, 2.141, 2.142, 2.143, 2.144, 2.144, 2.144, 2.145, 2.145, 2.145,
                                2.145, 2.146, 2.146, 2.146, 2.146};
    std::vector<double> col7 = {2.465, 2.490, 2.515, 2.539, 2.563, 2.586, 2.609, 2.616, 2.618, 2.620,
                                2.621, 2.623, 2.624, 2.625, 2.625, 2.626, 2.626, 2.627, 2.627, 2.627,
                                2.628, 2.628, 2.628, 2.628, 2.628};
    std::vector<std::vector<double>> tableY_all = {col1, col2, col3, col4, col5, col6, col7};
    std::vector<double> interp_Hrms(7, 0.0);
    for (size_t c = 0; c < 7; c++) {
        interp_Hrms[c] = cubicSplineInterpolation(tableX, tableY_all[c], Htr_tilde);
    }
    std::vector<double> dimensional(7, 0.0);
    for (size_t i = 0; i < 7; i++) {
        dimensional[i] = interp_Hrms[i] * Hrms;
    }
    double ratio_110_13    = (interp_Hrms[2] != 0.0) ? (interp_Hrms[3] / interp_Hrms[2]) : 0.0;
    double ratio_150_13    = (interp_Hrms[2] != 0.0) ? (interp_Hrms[4] / interp_Hrms[2]) : 0.0;
    double ratio_1100_13   = (interp_Hrms[2] != 0.0) ? (interp_Hrms[5] / interp_Hrms[2]) : 0.0;
    double ratio_11000_13  = (interp_Hrms[2] != 0.0) ? (interp_Hrms[6] / interp_Hrms[2]) : 0.0;
    double ratio_11000_110 = (interp_Hrms[3] != 0.0) ? (interp_Hrms[6] / interp_Hrms[3]) : 0.0;

    ss << L"====================\n";
    ss << L"   INPUT PARAMETERS\n";
    ss << L"====================\n";
    ss << L"Hm0 (m) : " << Hm0 << L"\n";
    ss << L"d (m)   : " << d << L"\n\n";
    
    ss << L"===========================\n";
    ss << L"   CALCULATED PARAMETERS\n";
    ss << L"===========================\n";
    ss << L"Free surface variance m0          : " << m0 << L"\n";
    ss << L"Mean square wave height Hrms (m)  : " << Hrms << L"\n";
    ss << L"Dimensional transitional wave height Htr : " << Htr_dim << L"\n";
    ss << L"Dimensionless H~_tr (Htr/Hrms)    : " << Htr_tilde << L"\n\n";
    
    ss << L"===============================================\n";
    ss << L"   DIMENSIONLESS WAVE HEIGHTS (H/Hrms) [Table]\n";
    ss << L"===============================================\n";
    ss << L"H1/Hrms       : " << interp_Hrms[0] << L"\n";
    ss << L"H2/Hrms       : " << interp_Hrms[1] << L"\n";
    ss << L"H1/3 / Hrms   : " << interp_Hrms[2] << L"\n";
    ss << L"H1/10 / Hrms  : " << interp_Hrms[3] << L"\n";
    ss << L"H1/50 / Hrms  : " << interp_Hrms[4] << L"\n";
    ss << L"H1/100 / Hrms : " << interp_Hrms[5] << L"\n";
    ss << L"H1/1000 /Hrms : " << interp_Hrms[6] << L"\n\n";
    
    ss << L"==================================\n";
    ss << L"   DIMENSIONAL WAVE HEIGHTS (m)\n";
    ss << L"==================================\n";
    ss << L"H1 (m)        : " << dimensional[0] << L"\n";
    ss << L"H2 (m)        : " << dimensional[1] << L"\n";
    ss << L"H1/3 (m)      : " << dimensional[2] << L"\n";
    ss << L"H1/10 (m)     : " << dimensional[3] << L"\n";
    ss << L"H1/50 (m)     : " << dimensional[4] << L"\n";
    ss << L"H1/100 (m)    : " << dimensional[5] << L"\n";
    ss << L"H1/1000 (m)   : " << dimensional[6] << L"\n\n";
    
    ss << L"===================\n";
    ss << L"   DIAGNOSTIC RATIOS\n";
    ss << L"===================\n";
    ss << L"(H1/10)/(H1/3)   : " << ratio_110_13 << L"\n";
    ss << L"(H1/50)/(H1/3)   : " << ratio_150_13 << L"\n";
    ss << L"(H1/100)/(H1/3)  : " << ratio_1100_13 << L"\n";
    ss << L"(H1/1000)/(H1/3) : " << ratio_11000_13 << L"\n";
    ss << L"(H1/1000)/(H1/10): " << ratio_11000_110 << L"\n\n";
    
    ss << L"End of Report\n";
    return ss.str();
}

//---------------------------------------------------------------------
// Main: Command-line interface entry point.
//---------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Set locale for proper wide-character output.
    std::wcout.imbue(std::locale(""));
    double Hm0 = 0.0, d = 0.0;
    if (argc >= 3) {
        Hm0 = std::stod(argv[1]);
        d = std::stod(argv[2]);
    } else {
        std::wcout << L"Enter Hm0 (m): ";
        std::wcin >> Hm0;
        std::wcout << L"Enter water depth d (m): ";
        std::wcin >> d;
    }
    std::wstring report = buildReport(Hm0, d);
    std::wcout << L"\n" << report << std::endl;
    return 0;
}
