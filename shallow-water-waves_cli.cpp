/***********************************************************************
 * Program: shallow-water-waves_cli.cpp
 *
 * Detailed Description:
 * This program computes local shallow-foreshore wave-height distribution
 * parameters using a model based on the Composed Weibull distribution.
 *
 * The command-line application performs the following:
 * 1. If three command-line arguments are provided, they are used as
 * Hm0 (local significant spectral wave height), d (local water depth),
 * and slopeM (beach slope 1:m). Otherwise, the program prompts the user
 * for these values.
 * 2. Computes the following intermediate values:
 * - Free-surface variance: m0 = (Hm0 / 4)²
 *
 * - Mean square wave height:
 * Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0)
 *
 * - A dimensional transitional wave height:
 * Htr = (0.35 + 5.8*(1/m)) * d.
 *
 * - The dimensionless transitional parameter:
 * H̃_tr = Htr / Hrms.
 * If H̃_tr is above 3.5 then it is set to 3.5 and Htr is recalculated as
 * Htr = 3.5 * Hrms.
 *
 * 3. Calculates the dimensionless wave-height ratios (Hᵢ/Hrms)
 * by solving a system of non-linear equations derived from the Composite Weibull
 * distribution, ensuring the normalized Hrms of the distribution equals one.
 * This involves using a Newton-Raphson matrix method for simultaneous root-finding
 * and functions for incomplete gamma calculations. These ratios are then
 * converted to dimensional quantities (in meters) by multiplying with Hrms.
 *
 * 4. A detailed report is then generated (and written to "report.txt") with
 * the input parameters, intermediate values, calculated ratios and computed
 * dimensional wave heights, as well as diagnostic ratios.
 *
 * Compilation Instructions (example using g++ on Windows with OpenMP):
 *
 * g++ -O3 -Wall shallow-water-waves_cli.cpp -o shallow-water-waves_cli ^
 * -static -static-libgcc -static-libstdc++
 *
 * To run with command-line arguments (e.g., Hm0=2.5, d=10, slopeM=20):
 * shallow-water-waves_cli 2.5 10 20
 ***********************************************************************/

#define _USE_MATH_DEFINES // Required on some systems (e.g., Windows with MSVC) for M_PI

#include <iostream>     // For input/output operations (e.g., std::cout, std::cerr)
#include <sstream>      // For string stream operations (e.g., std::wstringstream)
#include <vector>       // For dynamic array (std::vector)
#include <cmath>        // For mathematical functions (e.g., std::sqrt, std::log, std::pow, std::sin, std::lgamma, std::nan)
#include <limits>       // For numeric_limits (e.g., std::numeric_limits<double>::min())
#include <iomanip>      // For output formatting (e.g., std::fixed, std::setprecision)
#include <stdexcept>    // For std::invalid_argument
#include <fstream>      // For file output (std::wofstream)
#include <locale>       // For std::locale
#include <codecvt>      // For std::codecvt_utf8 (C++11, for UTF-8 file output)
#include <cstdlib>      // For std::stod (string to double)
#include <string>       // For std::string, std::wstring
#include <algorithm>    // For std::swap
#include <iterator>     // For std::size (C++17)

// Global parameters for the Composite Weibull distribution.
const double k1 = 2.0;    // Exponent for the first part of the Composite Weibull distribution (Rayleigh-shaped)
const double k2 = 3.6;    // Exponent for the second part of the Composite Weibull distribution

// Precision for Numerical Solver Convergence:
const double EPSILON = 1e-12; // A small value (10^-12) indicating the maximum allowable error or difference
const double JACOBIAN_DX = 1e-8; // Small step size for finite difference approximation of Jacobian derivatives.
const double LOCAL_EPS = 1e-16; // A very small tolerance value (10^-16) used for checking convergence
                                // within the gamma function calculations.

/**
 * @brief Computes the normalized lower incomplete gamma function P(a, x) = γ(a, x)/Γ(a).
 *
 * This function is a critical component for calculating moments of Weibull distributions.
 * It employs a hybrid approach for numerical stability and accuracy:
 * - For small values of 'x', it uses a series expansion.
 * - For larger values of 'x', it utilizes a continued fraction expansion.
 *
 * @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
 * @param x The 'x' parameter (upper integration limit) of the incomplete gamma function. Must be non-negative.
 * @return The computed value of P(a, x). Returns `nan("")` (Not-a-Number) if input parameters are invalid.
 */
double incomplete_gamma_p(double a, double x)
{
    const int MAXIT = 500;
    
    if(a <= 0.0 || x < 0.0)
        return std::nan("");
    if(x == 0.0)
        return 0.0;

    double gln = std::lgamma(a);
    
    if(x < a + 1.0) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int n_iter = 1; n_iter <= MAXIT; ++n_iter) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if(std::fabs(del) < std::fabs(sum) * LOCAL_EPS)
                return sum * std::exp(-x + a * std::log(x) - gln);
        }
        return std::nan("");
    } else {
        double b = x + 1.0 - a;
        double c = 1.0 / std::numeric_limits<double>::min();
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= MAXIT; ++i) {
            double an = -1.0 * i * (i - a);
            b += 2.0;
            d = an * d + b;
            if(std::fabs(d) < std::numeric_limits<double>::min())
                d = std::numeric_limits<double>::min();
            c = b + an / c;
            if(std::fabs(c) < std::numeric_limits<double>::min())
                c = std::numeric_limits<double>::min();
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if(std::fabs(del - 1.0) < LOCAL_EPS)
                break;
        }
        double lnPart = -x + a * std::log(x) - gln;
        double Qval = std::exp(lnPart) * h;
        double Pval = 1.0 - Qval;
        return (Pval < 0.0) ? 0.0 : Pval;
    }
}

/**
 * @brief Computes the normalized upper incomplete gamma function Q(a, x) = Γ(a, x)/Γ(a).
 *
 * This function is a simple wrapper that leverages the `incomplete_gamma_p` function.
 * By definition, Q(a, x) is the complement of P(a, x), i.e., 1 - P(a, x).
 *
 * @param a The 'a' parameter of the incomplete gamma function.
 * @param x The 'x' parameter of the incomplete gamma function.
 * @return The computed value of Q(a, x).
 */
inline double incomplete_gamma_q(double a, double x) {
    return 1.0 - incomplete_gamma_p(a, x);
}

/**
 * @brief Computes the unnormalized upper incomplete gamma function Γ(a, x).
 *
 * This function calculates the unnormalized form of the upper incomplete gamma function.
 * It is derived from the normalized upper incomplete gamma function Q(a, x) and the complete gamma function Γ(a).
 *
 * @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
 * @param x The 'x' parameter (lower integration limit) of the incomplete gamma function. Must be non-negative.
 * @return The computed value of the unnormalized upper incomplete gamma function Γ(a, x).
 * @throws `std::invalid_argument` if input parameters ('a' or 'x') are outside their valid ranges.
 */
double incomplete_gamma(double a, double x) {
    if (a <= 0.0) {
        throw std::invalid_argument("incomplete_gamma: 'a' must be positive.");
    }
    if (x < 0.0) {
        throw std::invalid_argument("incomplete_gamma: 'x' must be non-negative.");
    }
    return incomplete_gamma_q(a, x) * std::tgamma(a);
}

/**
 * @brief Calculates HN (wave height with 1/N exceedance probability) for the Composite Weibull distribution.
 *
 * This function determines the specific wave height (H) such that the probability of a wave
 * exceeding this height is 1/N. The calculation depends on whether the target wave height
 * falls into the first or second part of the Composite Weibull distribution, separated by Htr.
 *
 * @param N The N value (e.g., 3 for H1/3). Must be strictly greater than 1.0.
 * @param H1 The scale parameter of the first Weibull distribution. Must be positive.
 * @param H2 The scale parameter of the second Weibull distribution. Must be positive.
 * @param k1 The exponent of the first Weibull distribution. Must be positive.
 * @param k2 The exponent of the second Weibull distribution. Must be positive.
 * @param Htr The transitional wave height.
 * @return The calculated value of HN.
 * @throws `std::invalid_argument` if input parameters are invalid.
 */
double calculate_HN(double N, double H1, double H2, double k1, double k2, double Htr) {
    if (N <= 1.0) {
        throw std::invalid_argument("calculate_HN: N must be greater than 1 for log(N).");
    }
    if (H1 <= 0.0 || H2 <= 0.0) {
        throw std::invalid_argument("calculate_HN: H1 and H2 must be positive.");
    }
    if (k1 <= 0.0 || k2 <= 0.0) {
        throw std::invalid_argument("calculate_HN: k1 and k2 must be positive.");
    }

    double HN_candidate1 = H1 * std::pow(std::log(N), 1.0 / k1);

    if (HN_candidate1 < Htr - EPSILON) {
        return HN_candidate1;
    } else {
        return H2 * std::pow(std::log(N), 1.0 / k2);
    }
}

/**
 * @brief Calculates the mean of the highest 1/N-part of wave heights (H1/N) for the Composite Weibull distribution.
 *
 * This function computes a characteristic wave height that represents the average height of the
 * highest N-th fraction of waves. The calculation depends on whether the relevant wave heights
 * fall into the first or second part of the Composite Weibull distribution.
 *
 * @param N_val The N parameter for H1/N. Must be strictly greater than 1.
 * @param H1 The scale parameter of the first Weibull distribution. Must be positive.
 * @param H2 The scale parameter of the second Weibull distribution. Must be positive.
 * @param k1 The exponent of the first Weibull distribution. Must be positive.
 * @param k2 The exponent of the second Weibull distribution. Must be positive.
 * @param Htr The transitional wave height.
 * @return The calculated value of H1/N.
 * @throws `std::invalid_argument` if input parameters are invalid.
 */
double calculate_H1N(double N_val, double H1, double H2, double k1, double k2, double Htr) {
    if (H1 <= 0.0 || H2 <= 0.0) {
        throw std::invalid_argument("calculate_H1N: H1 and H2 must be positive.");
    }
    if (k1 <= 0.0 || k2 <= 0.0) {
        throw std::invalid_argument("calculate_H1N: k1 and k2 must be positive.");
    }
    if (N_val <= 1.0) {
        throw std::invalid_argument("calculate_H1N: N_val must be greater than 1.");
    }

    double H_N_val = calculate_HN(N_val, H1, H2, k1, k2, Htr);

    if (H_N_val < Htr - EPSILON) {
        double term1_a = 1.0 / k1 + 1.0;
        double term1_x_ln_Nval = std::log(N_val);
        double term1_x_HtrH1 = std::pow(Htr / H1, k1);

        double gamma_term1_part1 = incomplete_gamma(term1_a, term1_x_ln_Nval);
        double gamma_term1_part2 = incomplete_gamma(term1_a, term1_x_HtrH1);
        double gamma_term1 = gamma_term1_part1 - gamma_term1_part2;

        double term2_a = 1.0 / k2 + 1.0;
        double term2_x_HtrH2 = std::pow(Htr / H2, k2);

        double gamma_term2 = incomplete_gamma(term2_a, term2_x_HtrH2);

        return N_val * H1 * gamma_term1 + N_val * H2 * gamma_term2;
    } else {
        double term_a = 1.0 / k2 + 1.0;
        double term_x = std::log(N_val);
        return N_val * H2 * incomplete_gamma(term_a, term_x);
    }
}

// --- Functions for the System of Non-Linear Equations (Newton-Raphson Matrix Method) ---

/**
 * @brief Defines the first non-linear equation F1(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
 *
 * This equation represents the normalized Hrms constraint for the Composite Weibull distribution.
 * It states that the square root of the weighted sum of incomplete gamma functions must equal 1.
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms The normalized scale parameter of the second Weibull distribution.
 * @param Htr_Hrms The normalized transitional wave height (constant for a given solve).
 * @return The value of the first function, which should be driven to zero.
 */
double F1(double H1_Hrms, double H2_Hrms, double Htr_Hrms) {
    if (H1_Hrms <= 0.0 || H2_Hrms <= 0.0) {
        return std::numeric_limits<double>::max();
    }

    double arg1 = std::pow(Htr_Hrms / H1_Hrms, k1);
    double arg2 = std::pow(Htr_Hrms / H2_Hrms, k2);

    double term1 = H1_Hrms * H1_Hrms * incomplete_gamma_p(2.0 / k1 + 1.0, arg1);
    double term2 = H2_Hrms * H2_Hrms * incomplete_gamma_q(2.0 / k2 + 1.0, arg2);

    double sum_terms = term1 + term2;
    if (sum_terms < 0.0) sum_terms = 0.0;

    return std::sqrt(sum_terms) - 1.0;
}

/**
 * @brief Defines the second non-linear equation F2(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
 *
 * This equation represents the continuity condition between the two Weibull distributions
 * at the transitional wave height Htr.
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms The normalized scale parameter of the second Weibull distribution.
 * @param Htr_Hrms The normalized transitional wave height (constant for a given solve).
 * @return The value of the second function, which should be driven to zero.
 */
double F2(double H1_Hrms, double H2_Hrms, double Htr_Hrms) {
    if (H1_Hrms <= 0.0 || H2_Hrms <= 0.0) {
        return std::numeric_limits<double>::max();
    }
    return std::pow(Htr_Hrms / H1_Hrms, k1) - std::pow(Htr_Hrms / H2_Hrms, k2);
}

/**
 * @brief Solves a 2x2 linear system Ax = b for x using Cramer's rule.
 *
 * This function is a helper for the Newton-Raphson method for systems.
 *
 * @param J11 Element (1,1) of the Jacobian matrix (dF1/dH1).
 * @param J12 Element (1,2) of the Jacobian matrix (dF1/dH2).
 * @param J21 Element (2,1) of the Jacobian matrix (dF2/dH1).
 * @param J22 Element (2,2) of the Jacobian matrix (dF2/dH2).
 * @param b1 Right-hand side for the first equation (-F1).
 * @param b2 Right-hand side for the second equation (-F2).
 * @param dx1 Output: The calculated change for H1_Hrms.
 * @param dx2 Output: The calculated change for H2_Hrms.
 */
void solve_linear_system_2x2(double J11, double J12, double J21, double J22,
                             double b1, double b2, double &dx1, double &dx2) {
    double determinant = J11 * J22 - J12 * J21;

    if (std::abs(determinant) < std::numeric_limits<double>::epsilon() * 100) {
        dx1 = 0.0;
        dx2 = 0.0;
        return;
    }

    dx1 = (b1 * J22 - b2 * J12) / determinant;
    dx2 = (J11 * b2 - J21 * b1) / determinant;
}

/**
 * @brief Provides initial guesses for H1_Hrms and H2_Hrms based on Htr_Hrms.
 *
 * A good initial guess is crucial for the efficiency and robustness of the Newton-Raphson method.
 *
 * @param Htr_Hrms The normalized transitional wave height.
 * @param H1_initial Output: The initial guess for H1_Hrms.
 * @param H2_initial Output: The initial guess for H2_Hrms.
 */
void get_initial_guesses(double Htr_Hrms, double &H1_initial, double &H2_initial) {
    double logHtr = std::log(Htr_Hrms);
    H1_initial = 1.0 + std::exp(-1.705 - 2.329 * logHtr - 0.313 * logHtr * logHtr);

    H2_initial = Htr_Hrms * std::pow(H1_initial / Htr_Hrms, k1 / k2);

    if (H1_initial <= 0.0) H1_initial = std::numeric_limits<double>::min();
    if (H2_initial <= 0.0) H2_initial = std::numeric_limits<double>::min();
}

/**
 * @brief Solves for H1_Hrms and H2_Hrms simultaneously using the Newton-Raphson method for systems.
 *
 * This function implements the multi-dimensional Newton-Raphson algorithm to find the roots
 * of the system of non-linear equations F1 and F2.
 *
 * @param Htr_Hrms The normalized transitional wave height (constant for this solve).
 * @param H1_Hrms Output: The converged normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms Output: The converged normalized scale parameter of the second Weibull distribution.
 * @param tol The desired tolerance for convergence (maximum absolute value of F1 and F2).
 * @param maxit The maximum number of iterations allowed.
 * @return `true` if the solver successfully converges, `false` otherwise.
 */
bool newtonRaphsonSystemSolver(double Htr_Hrms, double &H1_Hrms, double &H2_Hrms,
                               double tol, int maxit) {
    get_initial_guesses(Htr_Hrms, H1_Hrms, H2_Hrms);

    for (int iter = 0; iter < maxit; ++iter) {
        double f1_val = F1(H1_Hrms, H2_Hrms, Htr_Hrms);
        double f2_val = F2(H1_Hrms, H2_Hrms, Htr_Hrms);

        if (std::abs(f1_val) < tol && std::abs(f2_val) < tol) {
            return true;
        }

        double J11 = (F1(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F1(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        double J12 = (F1(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F1(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        double J21 = (F2(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F2(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        double J22 = (F2(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F2(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0 * JACOBIAN_DX);

        double dH1, dH2;
        solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val, dH1, dH2);

        H1_Hrms += dH1;
        H2_Hrms += dH2;

        if (H1_Hrms <= 0.0) H1_Hrms = std::numeric_limits<double>::min();
        if (H2_Hrms <= 0.0) H2_Hrms = std::numeric_limits<double>::min();
    }

    std::wcerr << L"Newton-Raphson system solver failed to converge for Htr_Hrms = " << Htr_Hrms << L" after " << maxit << L" iterations.\n";
    return false;
}


//---------------------------------------------------------------------
// Function: buildReport
// Purpose:
//   Computes all wave parameters, performs calculations based on the
//   Composed Weibull model, and builds a detailed report as a
//   formatted wide string.
//---------------------------------------------------------------------
static std::wstring buildReport(double Hm0, double d, double slopeM)
{
    std::wstringstream ss;
    ss << std::fixed << std::setprecision(4);

    // Input validation for physical parameters
    if (Hm0 <= 0.0 || d <= 0.0)
    {
        ss << L"ERROR: Hm0 and d must be positive.\n";
        return ss.str();
    }
    if (slopeM <= 0.0)
    {
        ss << L"ERROR: Beach slope (m) must be positive.\n";
        return ss.str();
    }

    // Compute free-surface variance: m0 = (Hm0 / 4)^2
    double m0 = std::pow(Hm0 / 4.0, 2.0);

    // Mean square wave height: Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0)
    double Hrms = (2.69 + 3.24 * std::sqrt(m0) / d) * std::sqrt(m0);

    // Compute the actual tangent of the beach slope: tan(alpha) = 1.0 / slopeM.
    double tanAlpha = 1.0 / slopeM;
    // Dimensional transitional wave height: Htr = (0.35 + 5.8*(1/m)) * d.
    double Htr_dim = (0.35 + 5.8 * tanAlpha) * d;
    // Dimensionless transitional parameter: H̃_tr = Htr / Hrms.
    double Htr_tilde = (Hrms > 0.0) ? (Htr_dim / Hrms) : 0.0;

    // If dimensionless Htr exceeds 3.5, it is capped at 3.5 and Htr_dim is recalculated
    // to maintain consistency with the capped dimensionless value.
    if (Htr_tilde > 3.5)
    {
        Htr_tilde = 3.5;
        Htr_dim = 3.5 * Hrms;
    }

    // Solve for H1_Hrms and H2_Hrms simultaneously using the Newton-Raphson matrix solver.
    double H1_Hrms;
    double H2_Hrms;
    if (!newtonRaphsonSystemSolver(Htr_tilde, H1_Hrms, H2_Hrms, EPSILON, 100)) {
        ss << L"ERROR: Solver failed to converge for Htr_tilde = " << Htr_tilde << L". Cannot generate full report.\n";
        return ss.str();
    }

    // Define the N values for which H(1/N)/Hrms quantiles are calculated.
    const int N_values_for_H1N[] = {3, 10, 50, 100, 250, 1000};
    const size_t num_H1N_values = std::size(N_values_for_H1N);
    // Vector to store normalized wave heights (H1/Hrms, H2/Hrms, and H1/N values)
    std::vector<double> interp_Hrms(2 + num_H1N_values); 

    // Store H1/Hrms and H2/Hrms directly at the beginning of the vector.
    interp_Hrms[0] = H1_Hrms;
    interp_Hrms[1] = H2_Hrms;

    // Calculate H1/N values using the `calculate_H1N` function.
    for (size_t j = 0; j < num_H1N_values; ++j) {
        int N = N_values_for_H1N[j];
        interp_Hrms[j + 2] = calculate_H1N(static_cast<double>(N), H1_Hrms, H2_Hrms, k1, k2, Htr_tilde);
    }

    // Convert dimensionless values to dimensional wave heights (in meters)
    std::vector<double> dimensional(interp_Hrms.size(), 0.0);
    for (size_t i = 0; i < interp_Hrms.size(); i++)
    {
        dimensional[i] = interp_Hrms[i] * Hrms;
    }

    // Calculate diagnostic ratios.
    double ratio_110_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[3] / interp_Hrms[2]) : 0.0;
    double ratio_150_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[4] / interp_Hrms[2]) : 0.0;
    double ratio_1100_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[5] / interp_Hrms[2]) : 0.0;
    double ratio_1250_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[6] / interp_Hrms[2]) : 0.0;
    double ratio_11000_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[7] / interp_Hrms[2]) : 0.0;


    // Build the report string with formatted output.
    ss << L"======================\n";
    ss << L"   INPUT PARAMETERS\n";
    ss << L"======================\n";
    ss << L"Hm0 (m)         : " << Hm0 << L"\n";
    ss << L"d (m)           : " << d << L"\n";
    ss << L"Beach slope (1:m): " << slopeM << L"   (tan(alpha) = " << tanAlpha << L")\n\n";

    ss << L"===========================\n";
    ss << L"   CALCULATED PARAMETERS\n";
    ss << L"===========================\n";
    ss << L"Free surface variance m0         : " << m0 << L"\n";
    ss << L"Mean square wave height Hrms (m) : " << Hrms << L"\n";
    ss << L"Dimensionless H~_tr (Htr/Hrms)   : " << Htr_tilde << L"\n";
    ss << L"Transitional wave height Htr (m) : " << Htr_dim << L"\n\n";

    ss << L"=========================================\n";
    ss << L"   DIMENSIONLESS WAVE HEIGHTS (H/Hrms)\n";
    ss << L"=========================================\n";
    ss << L"H1/Hrms       : " << interp_Hrms[0] << L"\n";
    ss << L"H2/Hrms       : " << interp_Hrms[1] << L"\n";
    ss << L"H1/3 / Hrms   : " << interp_Hrms[2] << L"\n";
    ss << L"H1/10 / Hrms  : " << interp_Hrms[3] << L"\n";
    ss << L"H1/50 / Hrms  : " << interp_Hrms[4] << L"\n";
    ss << L"H1/100 / Hrms : " << interp_Hrms[5] << L"\n";
    ss << L"H1/250 / Hrms : " << interp_Hrms[6] << L"\n";
    ss << L"H1/1000 /Hrms : " << interp_Hrms[7] << L"\n\n";

    ss << L"==================================\n";
    ss << L"   DIMENSIONAL WAVE HEIGHTS (m)\n";
    ss << L"==================================\n";
    ss << L"H1 (m)        : " << dimensional[0] << L"\n";
    ss << L"H2 (m)        : " << dimensional[1] << L"\n";
    ss << L"H1/3 (m)      : " << dimensional[2] << L"\n";
    ss << L"H1/10 (m)     : " << dimensional[3] << L"\n";
    ss << L"H1/50 (m)     : " << dimensional[4] << L"\n";
    ss << L"H1/100 (m)    : " << dimensional[5] << L"\n";
    ss << L"H1/250 (m)    : " << dimensional[6] << L"\n";
    ss << L"H1/1000 (m)   : " << dimensional[7] << L"\n\n";

    ss << L"=======================\n";
    ss << L"   DIAGNOSTIC RATIOS\n";
    ss << L"=======================\n";
    ss << L"(H1/10)/(H1/3)   : " << ratio_110_13 << L"\n";
    ss << L"(H1/50)/(H1/3)   : " << ratio_150_13 << L"\n";
    ss << L"(H1/100)/(H1/3)  : " << ratio_1100_13 << L"\n";
    ss << L"(H1/250)/(H1/3)  : " << ratio_1250_13 << L"\n";
    ss << L"(H1/1000)/(H1/3) : " << ratio_11000_13 << L"\n";

    ss << L"End of Report\n";
    return ss.str();
}

//---------------------------------------------------------------------
// Function: writeReportToFile
// Purpose:
//   Writes the generated report string to a file named "report.txt".
//   The file is encoded in UTF-8 to support wide characters.
//---------------------------------------------------------------------
static void writeReportToFile(const std::wstring &report)
{
    std::wofstream ofs("report.txt");
    ofs.imbue(std::locale(std::locale(), new std::codecvt_utf8<wchar_t>));
    if (!ofs) {
        std::wcerr << L"Error: Could not open report.txt for writing.\n";
        return;
    }
    ofs << report;
    ofs.close();
}

//---------------------------------------------------------------------
// Main: Command-line interface entry point.
//---------------------------------------------------------------------
int main(int argc, char* argv[]) {
    std::wcout.imbue(std::locale(""));
    
    double Hm0 = 0.0;
    double d = 0.0;
    double slopeM = 0.0;

    if (argc >= 4) {
        try {
            Hm0 = std::stod(argv[1]);
            d = std::stod(argv[2]);
            slopeM = std::stod(argv[3]);
        } catch (const std::invalid_argument& ia) {
            std::wcerr << L"Invalid argument: " << ia.what() << L". Please enter numeric values.\n";
            return 1;
        } catch (const std::out_of_range& oor) {
            std::wcerr << L"Argument out of range: " << oor.what() << L". Value too large or too small.\n";
            return 1;
        }
    } else {
        std::wcout << L"Enter Hm0 (m): ";
        std::wcin >> Hm0;
        
        std::wcout << L"Enter water depth d (m): ";
        std::wcin >> d;
        
        std::wcout << L"Enter beach slope (1:m): ";
        std::wcin >> slopeM;
    }

    std::wstring report = buildReport(Hm0, d, slopeM);
    
    writeReportToFile(report);

    std::wcout << L"\n" << report << std::endl;

    return 0;
}
