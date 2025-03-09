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
 *   native Win32 API. It allows the user to input two key parameters:
 *
 *       1. Hm0 (in meters) - The local significant spectral wave height.
 *       2. d   (in meters) - The local water depth.
 *
 *   Based on these inputs, the program computes:
 *
 *       - Free-surface variance: m0 = (Hm0 / 4)²
 *       - Mean square wave height:
 *             Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0)
 *       - A dimensional transitional wave height:
 *             Htr = (0.12 * d / sqrt(m0)) * Hrms
 *       - The dimensionless transitional parameter H̃_tr = Htr/Hrms.
 *
 *   Using a 70-row table (with columns for H1/Hrms, H2/Hrms, etc.),
 *   a natural cubic spline interpolation is performed at H̃_tr to obtain
 *   dimensionless wave-height ratios. These are then converted to dimensional
 *   quantities (in meters) by multiplying with Hrms.
 *
 *   A detailed report is then generated (and written to "report.txt") with
 *   the input parameters, intermediate values, interpolated ratios and computed
 *   dimensional wave heights, as well as diagnostic ratios.
 *
 *   Compilation Instructions (example using g++ on Windows):
 *
 *         g++ -O2 -Wall -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui \
 *             -mwindows -static -static-libgcc -static-libstdc++
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

// Natural cubic spline functions
std::vector<double> computeSplineSecondDerivatives(const std::vector<double> &x, const std::vector<double> &y)
{
    size_t n = x.size();
    std::vector<double> m(n, 0.0);
    std::vector<double> l(n, 0.0), mu(n, 0.0), z(n, 0.0);

    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    std::vector<double> h(n - 1);
    for (size_t i = 0; i < n - 1; i++)
    {
        h[i] = x[i + 1] - x[i];
    }

    std::vector<double> alpha(n, 0.0);
    for (size_t i = 1; i < n - 1; i++)
    {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    for (size_t i = 1; i < n - 1; i++)
    {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    
    l[n - 1] = 1.0;
    z[n - 1] = 0.0;
    m[n - 1] = 0.0;

    for (int j = static_cast<int>(n) - 2; j >= 0; j--)
    {
        m[j] = z[j] - mu[j] * m[j + 1];
    }
    return m;
}

double cubicSplineInterpolation(const std::vector<double> &x, const std::vector<double> &y, double x0)
{
    size_t n = x.size();
    if (n == 0)
        throw std::runtime_error("No data points provided in x.");
    if (n == 1)
        return y[0];

    if (x0 <= x.front())
    {
        double t = (x0 - x.front()) / (x[1] - x.front());
        return y.front() + t * (y[1] - y.front());
    }

    if (x0 >= x.back())
    {
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

    double y0 = A * y[i] + B * y[i + 1] +
                ((A * A * A - A) * m[i] + (B * B * B - B) * m[i + 1]) * (h * h) / 6.0;
    return y0;
}

// Build the detailed report using the computed values.
static std::wstring buildReport(double Hm0, double d)
{
    std::wstringstream ss;
    ss << std::fixed << std::setprecision(4);

    if (Hm0 <= 0.0 || d <= 0.0)
    {
        ss << L"ERROR: Both Hm0 and d must be positive.\n";
        return ss.str();
    }

    // Compute free-surface variance
    double m0 = std::pow(Hm0 / 4.0, 2.0);
    // Updated calculation of Hrms:
    double Hrms = (2.69 + 3.24 * std::sqrt(m0) / d) * std::sqrt(m0);
    double Htr_dim = (0.12 * d / std::sqrt(m0)) * Hrms;
    double Htr_tilde = (Hrms > 0.0) ? (Htr_dim / Hrms) : 0.0;

    // Updated table: 70 rows from Htr/Hrms = 0.05 to 3.50 in increments of 0.05
    std::vector<double> tableX = {
         0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
         0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
         1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50,
         1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00,
         2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50,
         2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00,
         3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50};

    // New table columns from the provided table:
    // Column order: H1/Hrms, H2/Hrms, H1/3/Hrms, H1/10/Hrms, H1/50/Hrms, H1/100/Hrms, H1/1000/Hrms
    std::vector<double> col1 = {
        12.193, 7.003, 5.063, 4.022, 3.365, 2.908, 2.571, 2.311, 2.104, 1.936,
        1.796, 1.678, 1.578, 1.492, 1.419, 1.356, 1.302, 1.256, 1.216, 1.182,
        1.153, 1.128, 1.108, 1.090, 1.075, 1.063, 1.052, 1.043, 1.036, 1.030,
        1.024, 1.020, 1.016, 1.013, 1.011, 1.009, 1.007, 1.006, 1.004, 1.004,
        1.003, 1.002, 1.002, 1.001, 1.001, 1.001, 1.000, 1.000, 1.000, 1.000,
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000
    };

    std::vector<double> col2 = {
         1.060, 1.060, 1.060, 1.060, 1.060, 1.060, 1.060, 1.060, 1.060, 1.061,
         1.061, 1.062, 1.064, 1.066, 1.069, 1.073, 1.077, 1.083, 1.090, 1.097,
         1.106, 1.116, 1.126, 1.138, 1.150, 1.162, 1.175, 1.189, 1.203, 1.217,
         1.231, 1.246, 1.261, 1.275, 1.290, 1.305, 1.320, 1.334, 1.349, 1.363,
         1.378, 1.392, 1.407, 1.421, 1.435, 1.449, 1.462, 1.476, 1.490, 1.503,
         1.516, 1.529, 1.542, 1.555, 1.568, 1.580, 1.593, 1.605, 1.617, 1.630,
         1.642, 1.653, 1.665, 1.677, 1.689, 1.700, 1.711, 1.723, 1.734, 1.745
    };

    std::vector<double> col3 = {
         1.279, 1.279, 1.279, 1.279, 1.279, 1.279, 1.279, 1.279, 1.279, 1.280,
         1.281, 1.282, 1.284, 1.286, 1.290, 1.294, 1.300, 1.307, 1.315, 1.324,
         1.335, 1.346, 1.359, 1.371, 1.381, 1.389, 1.395, 1.399, 1.403, 1.406,
         1.408, 1.410, 1.411, 1.412, 1.413, 1.413, 1.414, 1.414, 1.415, 1.415,
         1.415, 1.415, 1.415, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416,
         1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416,
         1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416
    };

    std::vector<double> col4 = {
         1.466, 1.466, 1.466, 1.466, 1.466, 1.466, 1.466, 1.466, 1.466, 1.467,
         1.468, 1.469, 1.471, 1.474, 1.478, 1.483, 1.490, 1.498, 1.507, 1.518,
         1.530, 1.543, 1.558, 1.573, 1.590, 1.607, 1.626, 1.644, 1.664, 1.683,
         1.703, 1.721, 1.736, 1.749, 1.759, 1.767, 1.773, 1.779, 1.783, 1.786,
         1.789, 1.791, 1.793, 1.795, 1.796, 1.797, 1.798, 1.798, 1.799, 1.799,
         1.799, 1.800, 1.800, 1.800, 1.800, 1.800, 1.800, 1.800, 1.800, 1.800,
         1.800, 1.800, 1.800, 1.800, 1.800, 1.800, 1.800, 1.800, 1.800, 1.800
    };

    std::vector<double> col5 = {
         1.548, 1.548, 1.548, 1.548, 1.548, 1.548, 1.548, 1.548, 1.549, 1.549,
         1.550, 1.552, 1.554, 1.557, 1.561, 1.567, 1.573, 1.582, 1.591, 1.603,
         1.616, 1.630, 1.645, 1.662, 1.679, 1.698, 1.717, 1.737, 1.757, 1.778,
         1.799, 1.820, 1.841, 1.863, 1.884, 1.906, 1.927, 1.949, 1.970, 1.985,
         1.983, 1.982, 1.981, 1.980, 1.979, 1.978, 1.978, 1.978, 1.978, 1.978,
         1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978,
         1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978, 1.978
    };

    std::vector<double> col6 = {
         1.620, 1.620, 1.620, 1.620, 1.620, 1.620, 1.620, 1.620, 1.620, 1.621,
         1.622, 1.624, 1.626, 1.629, 1.634, 1.639, 1.646, 1.655, 1.665, 1.677,
         1.690, 1.705, 1.721, 1.739, 1.757, 1.776, 1.796, 1.817, 1.838, 1.860,
         1.882, 1.904, 1.927, 1.949, 1.972, 1.994, 2.017, 2.039, 2.062, 2.084,
         2.106, 2.128, 2.150, 2.149, 2.148, 2.148, 2.147, 2.147, 2.147, 2.147,
         2.146, 2.146, 2.146, 2.146, 2.146, 2.146, 2.146, 2.146, 2.146, 2.146,
         2.146, 2.146, 2.146, 2.146, 2.146, 2.146, 2.146, 2.146, 2.146, 2.146
    };

    std::vector<double> col7 = {
         1.813, 1.813, 1.813, 1.813, 1.813, 1.813, 1.813, 1.813, 1.813, 1.814,
         1.815, 1.817, 1.820, 1.823, 1.828, 1.835, 1.843, 1.852, 1.864, 1.877,
         1.892, 1.909, 1.927, 1.946, 1.967, 1.988, 2.011, 2.034, 2.058, 2.082,
         2.106, 2.131, 2.156, 2.182, 2.207, 2.232, 2.257, 2.282, 2.307, 2.332,
         2.357, 2.382, 2.406, 2.430, 2.454, 2.478, 2.502, 2.525, 2.548, 2.571,
         2.593, 2.616, 2.629, 2.629, 2.628, 2.628, 2.628, 2.628, 2.628, 2.628,
         2.628, 2.628, 2.628, 2.628, 2.628, 2.628, 2.628, 2.628, 2.628, 2.628
    };

    // Perform cubic spline interpolation on each column at H̃_tr.
    std::vector<double> interp_Hrms(7, 0.0);
    std::vector<std::vector<double>> tableY_all = {col1, col2, col3, col4, col5, col6, col7};
    for (size_t c = 0; c < 7; c++)
    {
        interp_Hrms[c] = cubicSplineInterpolation(tableX, tableY_all[c], Htr_tilde);
    }

    // Convert dimensionless values to dimensional wave heights.
    std::vector<double> dimensional(7, 0.0);
    for (size_t i = 0; i < 7; i++)
    {
        dimensional[i] = interp_Hrms[i] * Hrms;
    }

    // Calculate diagnostic ratios.
    double ratio_110_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[3] / interp_Hrms[2]) : 0.0;
    double ratio_150_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[4] / interp_Hrms[2]) : 0.0;
    double ratio_1100_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[5] / interp_Hrms[2]) : 0.0;
    double ratio_11000_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[6] / interp_Hrms[2]) : 0.0;
    double ratio_11000_110 = (interp_Hrms[3] != 0.0) ? (interp_Hrms[6] / interp_Hrms[3]) : 0.0;

    // Build the report string.
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
    ss << L"H1/3/Hrms     : " << interp_Hrms[2] << L"\n";
    ss << L"H1/10/Hrms    : " << interp_Hrms[3] << L"\n";
    ss << L"H1/50/Hrms    : " << interp_Hrms[4] << L"\n";
    ss << L"H1/100/Hrms   : " << interp_Hrms[5] << L"\n";
    ss << L"H1/1000/Hrms  : " << interp_Hrms[6] << L"\n\n";
    
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

// Write the report to a file (UTF-8 encoded).
static void writeReportToFile(const std::wstring &report)
{
    std::wofstream ofs("report.txt");
    ofs.imbue(std::locale(std::locale(), new std::codecvt_utf8<wchar_t>));
    if (!ofs)
        return;
    ofs << report;
}

// Convert newline characters for proper display in the edit control.
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
        else
        {
            out.push_back(c);
        }
    }
    return out;
}

#define IDC_EDIT_HM0 101
#define IDC_EDIT_D 102
#define IDC_BUTTON_COMPUTE 103
#define IDC_OUTPUT 104

HWND hEditHm0, hEditD, hOutput;

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    static HFONT hMonoFont = NULL;

    switch (msg)
    {
    case WM_CREATE:
    {
        CreateWindow(L"STATIC", L"Hm0 (m):", WS_CHILD | WS_VISIBLE,
                     10, 10, 80, 20, hwnd, NULL, NULL, NULL);
        hEditHm0 = CreateWindow(L"EDIT", L"1.5", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                100, 10, 100, 20, hwnd, (HMENU)IDC_EDIT_HM0, NULL, NULL);

        CreateWindow(L"STATIC", L"d (m):", WS_CHILD | WS_VISIBLE,
                     10, 40, 80, 20, hwnd, NULL, NULL, NULL);
        hEditD = CreateWindow(L"EDIT", L"10.0", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                              100, 40, 100, 20, hwnd, (HMENU)IDC_EDIT_D, NULL, NULL);

        CreateWindow(L"BUTTON", L"Compute", WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
                     10, 70, 190, 30, hwnd, (HMENU)IDC_BUTTON_COMPUTE, NULL, NULL);

        hOutput = CreateWindow(L"EDIT", L"", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL | ES_READONLY,
                               10, 110, 760, 440, hwnd, (HMENU)IDC_OUTPUT, NULL, NULL);

        hMonoFont = CreateFont(
            20, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
            DEFAULT_CHARSET,
            OUT_DEFAULT_PRECIS,
            CLIP_DEFAULT_PRECIS,
            DEFAULT_QUALITY, FIXED_PITCH | FF_DONTCARE, L"Courier New");
        SendMessage(hOutput, WM_SETFONT, (WPARAM)hMonoFont, TRUE);
    }
    break;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDC_BUTTON_COMPUTE)
        {
            wchar_t buffer[64];
            GetWindowText(hEditHm0, buffer, 63);
            double Hm0 = _wtof(buffer);

            GetWindowText(hEditD, buffer, 63);
            double d = _wtof(buffer);

            std::wstring report = buildReport(Hm0, d);
            writeReportToFile(report);
            std::wstring guiText = fixNewlinesForEditControl(report);
            SetWindowText(hOutput, guiText.c_str());
        }
        break;

    case WM_DESTROY:
        if (hMonoFont)
        {
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

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR lpCmdLine, int nCmdShow)
{
    const wchar_t CLASS_NAME[] = L"MyWindowClass";

    WNDCLASSEX wc;
    wc.cbSize = sizeof(WNDCLASSEX);
    wc.style = 0;
    wc.lpfnWndProc = WndProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = hInstance;
    wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.lpszMenuName = NULL;
    wc.lpszClassName = CLASS_NAME;
    wc.hIconSm = LoadIcon(NULL, IDI_APPLICATION);

    if (!RegisterClassEx(&wc))
    {
        MessageBox(NULL, L"Window Registration Failed!", L"Error!", MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

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
        MessageBox(NULL, L"Window Creation Failed!", L"Error!", MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return (int)msg.wParam;
}
