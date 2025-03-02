/***********************************************************************
 * Program: shallow-water-waves_gui.cpp
 *
 * Detailed Description:
 *   This program computes local shallow-foreshore wave-height distribution
 *   parameters using a model based on the Composed Weibull distribution as
 *   described in:
 *
 *       "Shallow foreshore wave height statistics"
 *        by H. Groenendijk, Master's Thesis, Delft University of Technology, 1998.
 *
 *   The program is implemented as a Windows GUI application using only the
 *   native Win32 API (without any external libraries). It allows the user to 
 *   input two key parameters:
 *
 *       1. Hm0 (in meters) - The local significant spectral wave height.
 *       2. d   (in meters) - The local water depth.
 *
 *   Based on these inputs, the program performs the following steps:
 *
 *   1. Input Acquisition:
 *       - The user enters the values for Hm0 and d into two separate edit
 *         controls on the main window.
 *
 *   2. Parameter Calculation:
 *       - Free-surface variance (m0) is computed as:
 *             m0 = (Hm0 / 4)²
 *       - The mean square wave height (Hrms) is computed as:
 *             Hrms = 3 * sqrt(m0)
 *       - A dimensional transitional wave height (Htr) is then calculated 
 *         using the corrected formula:
 *             Htr = (0.12 * d / sqrt(m0)) * Hrms
 *       - A dimensionless transitional parameter (H̃_tr) is derived:
 *             H̃_tr = Htr / Hrms
 *         This parameter is used as the interpolation point and must lie within
 *         the table range (2.3–3.5) for proper interpolation.
 *
 *   3. Interpolation of Wave-Height Ratios:
 *       - A predefined 25-row table is used:
 *             • tableX: Values ranging from 2.3 to 3.5 representing H̃_tr.
 *             • col1 to col7: The table columns contain characteristic ratios 
 *               relative to Hrms (i.e. Hᵢ/Hrms).
 *       - Natural cubic spline interpolation is performed for each column at 
 *         the value H̃_tr, thereby obtaining the dimensionless ratios Hᵢ/Hrms.
 *
 *   4. Conversion to Dimensional Quantities:
 *       - The dimensional wave heights (in meters) are then calculated using:
 *             H = (Hᵢ/Hrms) * Hrms
 *         where (Hᵢ/Hrms) are the interpolated values from the table.
 *
 *   5. Report Generation:
 *       - A detailed report is generated which includes:
 *             • The input parameters (Hm0 and d).
 *             • The computed intermediate values (m0, Hrms, Htr, H̃_tr).
 *             • The dimensionless wave heights (H/Hrms) as directly interpolated.
 *             • The dimensional wave heights (in meters) computed from Hrms.
 *             • Diagnostic ratios computed from the wave-height values.
 *
 *   6. Graphical User Interface (GUI):
 *       - The main window is designed as a non-resizable window.
 *       - It contains:
 *             • Two edit controls for inputting Hm0 and d.
 *             • A "Compute" button that triggers the parameter computations.
 *             • A multiline, read-only output control (using a 20-pt Courier New font)
 *               which displays the detailed report.
 *       - A helper routine ensures that Unix-style newline characters ("\n")
 *         are converted to the Windows CR-LF ("\r\n") format for proper display.
 *
 *   Compilation Instructions (Detailed):
 *       To compile this application using g++ on a Windows system, use a command
 *       similar to the following:
 *
 *         g++ -O2 -Wall -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui \
 *             -mwindows -static -static-libgcc -static-libstdc++
 *
 *       Explanation of the flags:
 *         -O2                  : Enables level 2 optimization for improved performance.
 *         -Wall                : Enables all compiler's warning messages to help with debugging.
 *         -municode            : Ensures Unicode support.
 *         -mwindows            : Links against the Windows subsystem rather than the console.
 *         -static              : Links statically to reduce dependency on DLLs.
 *         -static-libgcc      : Links the GCC runtime library statically.
 *         -static-libstdc++   : Links the standard C++ library statically.
 *
 *   References:
 *       H. Groenendijk, "Shallow foreshore wave height statistics", Master's Thesis,
 *       Delft University of Technology, 1998.
 *
 ***********************************************************************/

#include <windows.h>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <locale>
#include <codecvt>

/*========================================================================
 * Function: computeSplineSecondDerivatives
 *------------------------------------------------------------------------
 * Purpose:
 *   Given a set of data points (x, y), this function computes the second
 *   derivatives of the function y(x) using the natural cubic spline method.
 *
 * Detailed Process:
 *   - The method uses "natural" boundary conditions by setting the second
 *     derivatives at the endpoints to zero.
 *   - It computes the interval sizes (h[i] = x[i+1] - x[i]) for each adjacent
 *     pair of x values.
 *   - Coefficients (alpha) are computed using differences in y values scaled by
 *     the intervals.
 *   - A tridiagonal system is then solved using forward elimination and back
 *     substitution to determine the second derivatives.
 *
 * Parameters:
 *   const std::vector<double>& x : Vector containing the strictly increasing x values.
 *   const std::vector<double>& y : Vector containing the corresponding y values.
 *
 * Returns:
 *   std::vector<double> : A vector containing the computed second derivatives at each x.
 *
 * Exceptions:
 *   Throws a runtime error if no data points are provided.
 *========================================================================*/
std::vector<double> computeSplineSecondDerivatives(const std::vector<double>& x, const std::vector<double>& y)
{
    size_t n = x.size();
    std::vector<double> m(n, 0.0);  // Second derivatives vector.
    std::vector<double> l(n, 0.0), mu(n, 0.0), z(n, 0.0);

    // Set natural boundary conditions.
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    // Compute interval sizes h[i] = x[i+1] - x[i].
    std::vector<double> h(n - 1);
    for (size_t i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }

    // Compute the alpha coefficients based on differences in y.
    std::vector<double> alpha(n, 0.0);
    for (size_t i = 1; i < n - 1; i++) {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i])
                 - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    // Forward elimination to decompose the tridiagonal system.
    for (size_t i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    // Set boundary condition at the upper end.
    l[n - 1] = 1.0;
    z[n - 1] = 0.0;
    m[n - 1] = 0.0;

    // Back substitution to solve for the second derivatives.
    for (int j = static_cast<int>(n) - 2; j >= 0; j--) {
        m[j] = z[j] - mu[j] * m[j + 1];
    }
    return m;
}

/*========================================================================
 * Function: cubicSplineInterpolation
 *------------------------------------------------------------------------
 * Purpose:
 *   Performs natural cubic spline interpolation to estimate the y value at
 *   a given query point x0, based on provided data points.
 *
 * Detailed Process:
 *   1. If x0 is outside the provided x range, a simple linear extrapolation is
 *      performed using the closest interval.
 *   2. If x0 is within the range, the function locates the interval [x[i], x[i+1]]
 *      that contains x0.
 *   3. The weights A and B for the two endpoints of the interval are computed.
 *   4. The previously computed second derivatives are retrieved, and the
 *      natural cubic spline formula is used to compute the interpolated y value.
 *
 * Parameters:
 *   const std::vector<double>& x : Vector containing x values.
 *   const std::vector<double>& y : Vector containing corresponding y values.
 *   double x0                    : The query point at which interpolation is desired.
 *
 * Returns:
 *   double : The interpolated value y(x0).
 *
 * Exceptions:
 *   Throws a runtime error if no data points are provided.
 *========================================================================*/
double cubicSplineInterpolation(const std::vector<double>& x, const std::vector<double>& y, double x0)
{
    size_t n = x.size();
    if (n == 0)
        throw std::runtime_error("No data points provided in x.");
    if (n == 1)
        return y[0];

    // If x0 is to the left of the data range, use linear extrapolation.
    if (x0 <= x.front()) {
        double t = (x0 - x.front()) / (x[1] - x.front());
        return y.front() + t * (y[1] - y.front());
    }
    // If x0 is to the right of the data range, use linear extrapolation.
    if (x0 >= x.back()) {
        double t = (x0 - x[n - 2]) / (x.back() - x[n - 2]);
        return y[n - 2] + t * (y.back() - y[n - 2]);
    }

    // Locate the interval [x[i], x[i+1]] where x0 is located.
    size_t i = 0;
    while (i < n - 1 && x0 > x[i + 1])
        i++;

    double h = x[i + 1] - x[i];
    double A = (x[i + 1] - x0) / h;
    double B = (x0 - x[i]) / h;

    // Compute second derivatives using the natural cubic spline method.
    std::vector<double> m = computeSplineSecondDerivatives(x, y);

    // Evaluate the cubic spline interpolation formula.
    double y0 = A * y[i] + B * y[i + 1] +
                ((A * A * A - A) * m[i] + (B * B * B - B) * m[i + 1]) * (h * h) / 6.0;
    return y0;
}

/*========================================================================
 * Function: buildReport
 *------------------------------------------------------------------------
 * Purpose:
 *   Computes the shallow water wave-height distribution parameters and
 *   builds a comprehensive report.
 *
 * Detailed Process:
 *   a) Validates that the user-provided parameters (Hm0 and d) are positive.
 *   b) Computes key wave parameters:
 *         - m0    : Free-surface variance computed as (Hm0/4)^2.
 *         - Hrms  : Mean square wave height computed as 3 * sqrt(m0).
 *         - Htr   : Dimensional transitional wave height computed as
 *                   (0.12 * d / sqrt(m0)) * Hrms.
 *         - H̃_tr : Dimensionless transitional parameter computed as Htr/Hrms.
 *   c) A 25-row table is defined with:
 *         - tableX: Values from 2.3 to 3.5 representing H̃_tr.
 *         - col1 to col7: Predefined ratios representing Hᵢ/Hrms.
 *   d) For each column in the table, natural cubic spline interpolation is
 *      performed at x = H̃_tr to retrieve the dimensionless wave-height ratio (Hᵢ/Hrms).
 *   e) Dimensional wave heights are computed using:
 *         H = (Hᵢ/Hrms) * Hrms.
 *   f) Diagnostic ratios are computed from the set of interpolated dimensionless
 *      values.
 *   g) The report is assembled as a formatted wide string containing:
 *         - Input parameters.
 *         - Computed parameters.
 *         - Dimensionless wave heights (from table, H/Hrms).
 *         - Corresponding dimensional wave heights (in meters).
 *         - Diagnostic ratios.
 *
 * Parameters:
 *   double Hm0 : Local significant spectral wave height in meters.
 *   double d   : Local water depth in meters.
 *
 * Returns:
 *   std::wstring : A formatted report detailing all computations and results.
 *========================================================================*/
static std::wstring buildReport(double Hm0, double d)
{
    std::wstringstream ss;
    ss << std::fixed << std::setprecision(4);

    if (Hm0 <= 0.0 || d <= 0.0) {
        ss << L"ERROR: Both Hm0 and d must be positive.\n";
        return ss.str();
    }

    // Compute fundamental wave parameters.
    double m0 = std::pow(Hm0 / 4.0, 2.0);          // Free-surface variance.
    double Hrms = 3.0 * std::sqrt(m0);               // Mean square wave height.
    // Compute dimensional transitional wave height using the corrected formula.
    double Htr_dim = (0.12 * d / std::sqrt(m0)) * Hrms;
    double Htr_tilde = (Hrms > 0.0) ? (Htr_dim / Hrms) : 0.0; // Dimensionless transitional parameter.

    // Define the 25-row table with tableX values representing H̃_tr.
    std::vector<double> tableX = {
      2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75,
      2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.1, 3.15, 3.2, 3.25,
      3.3, 3.35, 3.4, 3.45, 3.5
    };

    // Table columns: Each column contains predefined ratios representing Hᵢ/Hrms.
    std::vector<double> col1 = {0.988, 0.989, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.997,
                                0.998, 0.998, 0.999, 0.999, 0.999, 0.999, 0.999, 1.000, 1.000, 1.000,
                                1.000, 1.000, 1.000, 1.000, 1.000};  // For H1/Hrms.
    std::vector<double> col2 = {1.419, 1.433, 1.448, 1.462, 1.475, 1.489, 1.502, 1.515, 1.528, 1.540,
                                1.553, 1.565, 1.577, 1.589, 1.601, 1.612, 1.623, 1.635, 1.646, 1.657,
                                1.668, 1.679, 1.689, 1.700, 1.711};  // For H2/Hrms.
    std::vector<double> col3 = {1.397, 1.400, 1.402, 1.404, 1.406, 1.407, 1.409, 1.410, 1.411, 1.412,
                                1.413, 1.413, 1.414, 1.414, 1.414, 1.415, 1.415, 1.415, 1.415, 1.415,
                                1.415, 1.416, 1.416, 1.416, 1.416};  // For H1/3/Hrms.
    std::vector<double> col4 = {1.774, 1.778, 1.781, 1.784, 1.786, 1.789, 1.791, 1.792, 1.794, 1.795,
                                1.796, 1.797, 1.797, 1.798, 1.798, 1.799, 1.799, 1.799, 1.799, 1.799,
                                1.800, 1.800, 1.800, 1.800, 1.800};  // For H1/10/Hrms.
    std::vector<double> col5 = {1.951, 1.954, 1.957, 1.960, 1.962, 1.965, 1.967, 1.969, 1.970, 1.971,
                                1.973, 1.974, 1.974, 1.975, 1.976, 1.976, 1.977, 1.977, 1.977, 1.977,
                                1.977, 1.978, 1.978, 1.978, 1.978};  // For H1/50/Hrms.
    std::vector<double> col6 = {2.116, 2.120, 2.123, 2.126, 2.129, 2.132, 2.134, 2.136, 2.137, 2.139,
                                2.140, 2.141, 2.142, 2.143, 2.144, 2.144, 2.144, 2.145, 2.145, 2.145,
                                2.145, 2.146, 2.146, 2.146, 2.146};  // For H1/100/Hrms.
    std::vector<double> col7 = {2.465, 2.490, 2.515, 2.539, 2.563, 2.586, 2.609, 2.616, 2.618, 2.620,
                                2.621, 2.623, 2.624, 2.625, 2.625, 2.626, 2.626, 2.627, 2.627, 2.627,
                                2.628, 2.628, 2.628, 2.628, 2.628};  // For H1/1000/Hrms.

    // Pack the table columns into a vector-of-vectors for iterative processing.
    std::vector<std::vector<double>> tableY_all = {col1, col2, col3, col4, col5, col6, col7};

    // Interpolate each column at the query point H̃_tr to retrieve Hᵢ/Hrms.
    std::vector<double> interp_Hrms(7, 0.0);
    for (size_t c = 0; c < 7; c++) {
        interp_Hrms[c] = cubicSplineInterpolation(tableX, tableY_all[c], Htr_tilde);
    }
    
    // Compute the dimensional wave heights (in meters) using the relation:
    //      H = (Hᵢ/Hrms) * Hrms.
    std::vector<double> dimensional(7, 0.0);
    for (size_t i = 0; i < 7; i++) {
        dimensional[i] = interp_Hrms[i] * Hrms;
    }
    
    // Compute diagnostic ratios from the interpolated values.
    double ratio_110_13    = (interp_Hrms[2] != 0.0) ? (interp_Hrms[3] / interp_Hrms[2]) : 0.0;
    double ratio_150_13    = (interp_Hrms[2] != 0.0) ? (interp_Hrms[4] / interp_Hrms[2]) : 0.0;
    double ratio_1100_13   = (interp_Hrms[2] != 0.0) ? (interp_Hrms[5] / interp_Hrms[2]) : 0.0;
    double ratio_11000_13  = (interp_Hrms[2] != 0.0) ? (interp_Hrms[6] / interp_Hrms[2]) : 0.0;
    double ratio_11000_110 = (interp_Hrms[3] != 0.0) ? (interp_Hrms[6] / interp_Hrms[3]) : 0.0;
    
    // Assemble the report as a formatted string.
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
    
    // Report the interpolated dimensionless wave heights (H/Hrms).
    ss << L"===============================================\n";
    ss << L"   DIMENSIONLESS WAVE HEIGHTS (H/Hrms) [Table]\n";
    ss << L"===============================================\n";
    ss << L"H1/Hrms        : " << interp_Hrms[0] << L"\n";
    ss << L"H2/Hrms        : " << interp_Hrms[1] << L"\n";
    ss << L"H1/3 / Hrms   : " << interp_Hrms[2] << L"\n";
    ss << L"H1/10 / Hrms  : " << interp_Hrms[3] << L"\n";
    ss << L"H1/50 / Hrms  : " << interp_Hrms[4] << L"\n";
    ss << L"H1/100 / Hrms : " << interp_Hrms[5] << L"\n";
    ss << L"H1/1000 /Hrms : " << interp_Hrms[6] << L"\n\n";
    
    // Report the computed dimensional wave heights (in meters).
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
    
    // Report the diagnostic ratios computed from the wave-height values.
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

/*========================================================================
 * Function: writeReportToFile
 *------------------------------------------------------------------------
 * Purpose:
 *   Writes the generated report to a file named "report.txt" using UTF-8
 *   encoding. This allows for easy sharing and external viewing of the results.
 *
 * Detailed Process:
 *   - Opens an output file stream using std::wofstream.
 *   - Imbues the stream with a UTF-8 codec facet to ensure proper encoding.
 *   - Writes the report string to the file.
 *
 * Parameters:
 *   const std::wstring &report : The complete report to be written.
 *
 * Note:
 *   If the file cannot be opened, the routine simply returns without error.
 *========================================================================*/
static void writeReportToFile(const std::wstring &report)
{
    std::wofstream ofs("report.txt");
    ofs.imbue(std::locale(std::locale(), new std::codecvt_utf8<wchar_t>));
    if (!ofs)
        return;
    ofs << report;
}

/*========================================================================
 * Function: fixNewlinesForEditControl
 *------------------------------------------------------------------------
 * Purpose:
 *   Converts newline characters in the given string from Unix-style ("\n")
 *   to Windows-style CR-LF ("\r\n"). This ensures proper display in a Windows
 *   multiline edit control.
 *
 * Detailed Process:
 *   - Iterates through each character in the input string.
 *   - Whenever a newline character is encountered, it appends the CR-LF sequence.
 *   - Otherwise, it copies the character unchanged.
 *
 * Parameters:
 *   const std::wstring &text : The original string containing Unix-style newlines.
 *
 * Returns:
 *   std::wstring : A new string with Windows-style line breaks.
 *========================================================================*/
static std::wstring fixNewlinesForEditControl(const std::wstring &text)
{
    std::wstring out;
    out.reserve(text.size() + 100);
    for (wchar_t c : text)
    {
        if (c == L'\n')
        {
            out.push_back(L'\r');
            out.push_back(L'\n');
        }
        else {
            out.push_back(c);
        }
    }
    return out;
}

/*========================================================================
 * GUI Definitions:
 *   The following section defines unique IDs for the GUI controls and
 *   declares global variables for the window handles.
 *
 *   Control IDs:
 *     - IDC_EDIT_HM0       : Identifier for the Hm0 (wave height) input edit control.
 *     - IDC_EDIT_D         : Identifier for the water depth (d) input edit control.
 *     - IDC_BUTTON_COMPUTE : Identifier for the "Compute" button.
 *     - IDC_OUTPUT         : Identifier for the output (report) edit control.
 *
 *   Global HWND Variables:
 *     - hEditHm0 : Handle to the Hm0 input edit control.
 *     - hEditD   : Handle to the water depth input edit control.
 *     - hOutput  : Handle to the output edit control.
 *========================================================================*/
#define IDC_EDIT_HM0       101
#define IDC_EDIT_D         102
#define IDC_BUTTON_COMPUTE 103
#define IDC_OUTPUT         104

HWND hEditHm0, hEditD, hOutput;

/*========================================================================
 * Function: WndProc
 *------------------------------------------------------------------------
 * Purpose:
 *   The window procedure processes messages for the main application window.
 *   It handles the creation of controls, command messages from user interactions,
 *   and cleanup during window destruction.
 *
 * Detailed Process:
 *   - WM_CREATE:
 *         * Creates static labels and input edit controls for Hm0 and d.
 *         * Sets default values for Hm0 ("1.5") and d ("10.0") such that H̃_tr
 *           lies within the table range.
 *         * Creates the "Compute" button which initiates the calculation process.
 *         * Creates a multiline, read-only output edit control with vertical
 *           scrolling and assigns it a 20-pt Courier New font.
 *
 *   - WM_COMMAND:
 *         * When the "Compute" button is clicked, the procedure:
 *              1. Retrieves the user-entered Hm0 and d values.
 *              2. Calls buildReport to generate the detailed output report.
 *              3. Writes the report to "report.txt" (UTF-8 encoded).
 *              4. Converts newline characters for proper display in the output control.
 *              5. Displays the report in the output edit control.
 *
 *   - WM_DESTROY:
 *         * Cleans up any resources (such as the font) and posts a quit message.
 *
 * Parameters:
 *   hwnd   : Handle to the current window.
 *   msg    : Message identifier.
 *   wParam : Additional message-specific information.
 *   lParam : Additional message-specific information.
 *
 * Returns:
 *   LRESULT : Result of message processing.
 *========================================================================*/
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    static HFONT hMonoFont = NULL;

    switch(msg)
    {
    case WM_CREATE:
        {
            // Create label and edit control for Hm0 input.
            CreateWindow(L"STATIC", L"Hm0 (m):", WS_CHILD | WS_VISIBLE,
                         10, 10, 80, 20, hwnd, NULL, NULL, NULL);
            // Default value "1.5" m chosen to yield H̃_tr within the table range.
            hEditHm0 = CreateWindow(L"EDIT", L"1.5", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                     100, 10, 100, 20, hwnd, (HMENU)IDC_EDIT_HM0, NULL, NULL);

            // Create label and edit control for water depth (d) input.
            CreateWindow(L"STATIC", L"d (m):", WS_CHILD | WS_VISIBLE,
                         10, 40, 80, 20, hwnd, NULL, NULL, NULL);
            // Default value "10.0" m.
            hEditD = CreateWindow(L"EDIT", L"10.0", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                  100, 40, 100, 20, hwnd, (HMENU)IDC_EDIT_D, NULL, NULL);

            // Create the "Compute" button.
            CreateWindow(L"BUTTON", L"Compute", WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
                         10, 70, 190, 30, hwnd, (HMENU)IDC_BUTTON_COMPUTE, NULL, NULL);

            // Create the multiline, read-only output edit control with vertical scrolling.
            hOutput = CreateWindow(L"EDIT", L"", WS_CHILD | WS_VISIBLE | WS_BORDER |
                                   ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL | ES_READONLY,
                                   10, 110, 760, 440, hwnd, (HMENU)IDC_OUTPUT, NULL, NULL);

            // Create a 20-pt Courier New font for clear, monospaced output.
            hMonoFont = CreateFont(
                20, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                DEFAULT_CHARSET,
                OUT_DEFAULT_PRECIS,
                CLIP_DEFAULT_PRECIS,
                DEFAULT_QUALITY, FIXED_PITCH | FF_DONTCARE, L"Courier New"
            );
            SendMessage(hOutput, WM_SETFONT, (WPARAM)hMonoFont, TRUE);
        }
        break;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDC_BUTTON_COMPUTE)
        {
            // Retrieve Hm0 input from the corresponding edit control.
            wchar_t buffer[64];
            GetWindowText(hEditHm0, buffer, 63);
            double Hm0 = _wtof(buffer);

            // Retrieve water depth (d) input from its edit control.
            GetWindowText(hEditD, buffer, 63);
            double d = _wtof(buffer);

            // Build the detailed report using the input parameters.
            std::wstring report = buildReport(Hm0, d);
            // Save the report to a file ("report.txt") with UTF-8 encoding.
            writeReportToFile(report);
            // Convert newline characters for proper display in the Windows edit control.
            std::wstring guiText = fixNewlinesForEditControl(report);
            // Display the report in the output control.
            SetWindowText(hOutput, guiText.c_str());
        }
        break;

    case WM_DESTROY:
        if (hMonoFont) {
            DeleteObject(hMonoFont);
            hMonoFont = NULL;
        }
        PostQuitMessage(0);
        break;

    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
    return 0;
}

/*========================================================================
 * Function: wWinMain
 *------------------------------------------------------------------------
 * Purpose:
 *   Entry point for the Windows GUI application.
 *
 * Detailed Process:
 *   - Defines and registers a window class named "MyWindowClass" with the
 *     window procedure (WndProc) that handles all messages.
 *   - Creates the main application window using the registered window class.
 *   - The window is created as non-resizable by omitting the resizing border
 *     and maximize box.
 *   - The window is displayed and enters a message loop that processes
 *     incoming events until termination.
 *
 * Parameters:
 *   HINSTANCE hInstance : Handle to the current instance.
 *   HINSTANCE hPrevInstance: Handle to the previous instance (not used).
 *   PWSTR lpCmdLine     : Command-line arguments.
 *   int nCmdShow        : Controls how the window is to be shown.
 *
 * Returns:
 *   int : Exit code of the application.
 *========================================================================*/
int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
                    PWSTR lpCmdLine, int nCmdShow)
{
    const wchar_t CLASS_NAME[] = L"MyWindowClass";

    WNDCLASSEX wc;
    wc.cbSize        = sizeof(WNDCLASSEX);
    wc.style         = 0;
    wc.lpfnWndProc   = WndProc;
    wc.cbClsExtra    = 0;
    wc.cbWndExtra    = 0;
    wc.hInstance     = hInstance;
    wc.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
    wc.lpszMenuName  = NULL;
    wc.lpszClassName = CLASS_NAME;
    wc.hIconSm       = LoadIcon(NULL, IDI_APPLICATION);

    // Register the window class.
    if (!RegisterClassEx(&wc))
    {
        MessageBox(NULL, L"Window Registration Failed!", L"Error!",
                   MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    // Create a non-resizable window by removing the thick frame and maximize box.
    DWORD style = (WS_OVERLAPPEDWINDOW & ~(WS_THICKFRAME | WS_MAXIMIZEBOX));
    HWND hwnd = CreateWindowEx(
        WS_EX_CLIENTEDGE,
        CLASS_NAME,
        L"Shallow Wave Heights GUI",
        style,
        CW_USEDEFAULT, CW_USEDEFAULT, 820, 640,
        NULL, NULL, hInstance, NULL);

    if (!hwnd)
    {
        MessageBox(NULL, L"Window Creation Failed!", L"Error!",
                   MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    // Main message loop: processes all incoming messages.
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return (int)msg.wParam;
}
