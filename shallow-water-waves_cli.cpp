/***********************************************************************
 * Program: shallow-water-waves_cli.cpp
 *
 * Detailed Description:
 * This program computes local shallow-foreshore wave-height distribution
 * parameters using a model based on the Composed Weibull distribution.
 *
 * The command-line application performs the following:
 * 1. If four command-line arguments are provided, they are used as
 * Hm0 (local significant spectral wave height), m0 (free-surface variance),
 * d (local water depth), and slopeM (beach slope 1:m). Otherwise, the program
 * prompts the user for these values.
 * 2. Computes the following intermediate values:
 * - Mean square wave height:
 * Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0)
 *
 * - A dimensional transitional wave height:
 * Htr = (0.35 + 5.8*(1/m)) * d.
 *
 * - The dimensionless transitional parameter:
 * H̃_tr = Htr / Hrms.
 *
 * 3. Calculates the dimensionless wave-height ratios (Hᵢ/Hrms)
 * by solving a system of non-linear equations derived from the Composite Weibull
 * distribution, ensuring the normalized Hrms of the distribution equals one.
 * This involves using a Newton-Raphson matrix method for simultaneous root-finding
 * and functions for unnormalized incomplete gamma calculations. These ratios are then
 * converted to dimensional quantities (in meters) by multiplying with Hrms.
 *
 * 4. A detailed report is then generated (and written to "report.txt") with
 * the input parameters, intermediate values, calculated ratios and computed
 * dimensional wave heights, as well as diagnostic ratios.
 *
 * Compilation Instructions (example using g++ on Windows with OpenMP):
 *
 * g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic \
 * -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
 * -o shallow-water-waves_cli shallow-water-waves_cli.cpp
 *
 * To run with command-line arguments (e.g., Hm0=2.5, m0=0.3906, d=5, slopeM=100):
 * shallow-water-waves_cli 2.5 0.3906 5 100
 ***********************************************************************/

#define _USE_MATH_DEFINES // Required on some systems (e.g., Windows with MSVC) for M_PI

#include <iostream>     // For input/output operations (e.g., std::cout, std::cerr)
#include <sstream>      // For string stream operations (e.g., std::wstringstream)
#include <vector>       // For dynamic array (std::vector)
#include <cmath>        // For mathematical functions (e.g., std::sqrt, std::log, std::pow, std::sin, std::lgamma, std::nan, std::tgamma)
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
 * @brief Computes the unnormalized lower incomplete gamma function γ(a, x) = ∫[0,x] t^(a-1)e^(-t) dt.
 *
 * This function is a critical component for calculating moments of Weibull distributions,
 * which are fundamental to the Composite Weibull model. It employs a hybrid approach
 * for numerical stability and accuracy:
 * - For small values of 'x' (specifically, x < a + 1.0), it uses a series expansion.
 * - For larger values of 'x', it utilizes a continued fraction expansion.
 * This adaptive strategy ensures robust and precise computation across different input ranges.
 *
 * @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
 * @param x The 'x' parameter (upper integration limit) of the incomplete gamma function. Must be non-negative.
 * @return The computed value of γ(a, x). Returns `nan("")` (Not-a-Number) if input parameters are invalid
 * or if the series/continued fraction fails to converge within the maximum iterations.
 */
double incomplete_gamma_lower(double a, double x)
{
    const int MAXIT = 500;      // Maximum number of iterations allowed for either the series expansion or the continued fraction.

    // Input validation: Essential for preventing mathematical errors and ensuring function robustness.
    if(a <= 0.0 || x < 0.0)
        return std::nan(""); // If 'a' is not positive or 'x' is negative, return NaN as the result is undefined.
    if(x == 0.0)
        return 0.0; // By definition, γ(a, 0) is 0 for any positive 'a'.

    // Compute the natural logarithm of the complete gamma function, ln(Γ(a)).
    // `std::lgamma` is used instead of `std::tgamma` followed by `log` for improved numerical stability,
    // especially when 'a' is large, as `lgamma` avoids potential overflow/underflow issues with very large/small gamma values.
    double gln = std::lgamma(a);

    // Conditional logic to select the appropriate numerical method (series or continued fraction).
    // This choice is based on the relationship between 'x' and 'a', a common heuristic for optimal performance and accuracy.
    if(x < a + 1.0) {  // Use series expansion for x < a+1.
        double ap = a;          // `ap` is a mutable copy of 'a', used to increment in the series terms.
        double sum = 1.0 / a;   // Initialize the sum with the first term of the series expansion.
        double del = sum;       // `del` stores the value of the current term being added to the sum.
        for (int n_iter = 1; n_iter <= MAXIT; ++n_iter) { // Iterate up to MAXIT to compute series terms.
            ap += 1.0;          // Increment 'a' for the next term in the series (a+1, a+2, ...).
            del *= x / ap;      // Calculate the next term: `del_n = del_{n-1} * (x / (a + n))`.
            sum += del;         // Add the newly calculated term to the cumulative sum.
            // Check for convergence: If the absolute value of the current term (`del`) is extremely small
            // relative to the absolute value of the cumulative sum (`sum`), then the series has converged.
            if(std::abs(del) < std::abs(sum) * LOCAL_EPS)
                // Return the final result for γ(a,x).
                return sum * std::exp(-x + a * std::log(x) - gln) * std::tgamma(a);
        }
        return std::nan(""); // If the loop finishes without `del` becoming sufficiently small, it means convergence failed.
    } else {  // Use continued fraction for x >= a+1.
        double b = x + 1.0 - a; // Initialize the 'b' term (B_0) for the continued fraction expansion.
        // Initialize 'c' and 'd' to extremely large/small values. This is a common technique
        // to prevent division by zero or underflow during the initial steps of the continued fraction
        // algorithm, where denominators might otherwise be zero or very close to zero.
        double c = 1.0 / std::numeric_limits<double>::min();
        double d = 1.0 / b;
        double h = d; // `h` accumulates the value of the continued fraction.
        for (int i = 1; i <= MAXIT; ++i) { // Iterate up to MAXIT to compute continued fraction terms.
            double an = -1.0 * i * (i - a); // Calculate the 'a_n' term (numerator) for the current iteration.
            b += 2.0; // Update the 'b_n' term (denominator) for the current iteration.
            d = an * d + b; // Apply the recurrence relation for D_n (denominator part of the fraction).
            // Safeguard against division by zero or extremely small numbers for 'd'.
            if(std::abs(d) < std::numeric_limits<double>::min())
                d = std::numeric_limits<double>::min();
            c = b + an / c; // Apply the recurrence relation for C_n (numerator part of the fraction).
            // Safeguard against division by zero or extremely small numbers for 'c'.
            if(std::abs(c) < std::numeric_limits<double>::min())
                c = std::numeric_limits<double>::min();
            d = 1.0 / d;    // Invert 'd' for the next step.
            double del = d * c; // Calculate the current term (delta_n) of the continued fraction.
            h *= del;       // Multiply `h` by the current term to update the accumulated fraction value.
            // Check for convergence: If the current term `del` is very close to 1.0, the continued fraction has converged.
            if(std::abs(del - 1.0) < LOCAL_EPS)
                break; // Exit the loop if converged.
        }
        // This part of the calculation results in Q(a,x) before multiplication by Gamma(a).
        double lnPart = -x + a * std::log(x) - gln;
        double Qval_normalized = std::exp(lnPart) * h;
        // P(a,x) = 1 - Q(a,x).
        double Pval_normalized = 1.0 - Qval_normalized;
        // It's clamped to 0.0 to prevent very small negative floating-point results due to precision issues.
        Pval_normalized = (Pval_normalized < 0.0) ? 0.0 : Pval_normalized;
        // Return unnormalized lower incomplete gamma γ(a,x) = P(a,x) * Γ(a).
        return Pval_normalized * std::tgamma(a);
    }
}


/**
 * @brief Computes the unnormalized upper incomplete gamma function Γ(a, x) = ∫[x,∞] t^(a-1)e^(-t) dt.
 *
 * This function calculates the unnormalized form of the upper incomplete gamma function.
 * It is derived from the complete gamma function Γ(a) and the unnormalized lower incomplete gamma function γ(a, x),
 * using the identity Γ(a,x) = Γ(a) - γ(a,x).
 *
 * @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
 * @param x The 'x' parameter (lower integration limit) of the incomplete gamma function. Must be non-negative.
 * @return The computed value of the unnormalized upper incomplete gamma function Γ(a, x).
 * @throws `std::invalid_argument` if input parameters ('a' or 'x') are outside their valid ranges.
 */
double incomplete_gamma_upper(double a, double x) {
    // Input validation: Ensures 'a' is positive and 'x' is non-negative, as required by the gamma function definitions.
    if (a <= 0.0) {
        throw std::invalid_argument("incomplete_gamma_upper: 'a' must be positive.");
    }
    if (x < 0.0) {
        throw std::invalid_argument("incomplete_gamma_upper: 'x' must be non-negative.");
    }
    // Calculation: Γ(a,x) = Γ(a) - γ(a,x).
    // `std::tgamma` computes the complete gamma function Γ(a).
    return std::tgamma(a) - incomplete_gamma_lower(a, x);
}


/**
 * @brief Calculates HN (wave height with 1/N exceedance probability) for the Composite Weibull distribution.
 *
 * This function determines the specific wave height (H) such that the probability of a wave
 * exceeding this height is 1/N. This is a key statistical measure. The calculation depends
 * on whether the target wave height falls into the first or second part of the Composite Weibull
 * distribution, which is separated by the transitional wave height (Htr).
 *
 * @param N The N value (e.g., 3 for H1/3, 100 for H1%). Must be strictly greater than 1.0
 * because `log(N)` is used, and `N=1` would result in `log(1)=0`.
 * @param H1 The scale parameter of the first Weibull distribution. Must be positive.
 * @param H2 The scale parameter of the second Weibull distribution. Must be positive.
 * @param k1 The exponent (shape parameter) of the first Weibull distribution. Must be positive.
 * @param k2 The exponent (shape parameter) of the second Weibull distribution. Must be positive.
 * @param Htr The transitional wave height, which defines the boundary between the two Weibull parts.
 * @return The calculated value of HN.
 * @throws `std::invalid_argument` if input parameters are invalid (e.g., N <= 1, H1/H2 <= 0, k1/k2 <= 0).
 */
double calculate_HN(double N, double H1, double H2, double k1, double k2, double Htr) {
    // Input validation: Ensures all parameters are within their mathematically valid ranges.
    if (N <= 1.0) {
        throw std::invalid_argument("calculate_HN: N must be greater than 1 for log(N).");
    }
    if (H1 <= 0.0 || H2 <= 0.0) {
        throw std::invalid_argument("calculate_HN: H1 and H2 must be positive.");
    }
    if (k1 <= 0.0 || k2 <= 0.0) {
        throw std::invalid_argument("calculate_HN: k1 and k2 must be positive.");
    }

    // Calculate a candidate HN assuming it falls within the *first* part of the distribution (H <= Htr).
    double HN_candidate1 = H1 * std::pow(std::log(N), 1.0 / k1);

    // Decision point: Check if the `HN_candidate1` is indeed less than the transitional wave height (Htr).
    // A small `EPSILON` is subtracted from `Htr` to account for floating-point inaccuracies when comparing.
    if (HN_candidate1 < Htr - EPSILON) {
        // If true, it means the wave height for the given exceedance probability is governed by the first Weibull part.
        return HN_candidate1;
    } else {
        // If `HN_candidate1` is not less than `Htr`, it implies that the wave height for the given
        // exceedance probability is determined by the *second* part of the distribution (H > Htr).
        return H2 * std::pow(std::log(N), 1.0 / k2);
    }
}

/**
 * @brief Calculates the mean of the highest 1/N-part of wave heights (H1/N) for the Composite Weibull distribution.
 *
 * This function computes a characteristic wave height that represents the average height of the
 * highest N-th fraction of waves in a given wave field. This is a commonly used metric in
 * oceanography and coastal engineering (e.g., H1/3 for significant wave height).
 * The calculation involves integrals of the probability density function and depends on
 * whether the relevant wave heights fall into the first or second part of the Composite Weibull distribution.
 *
 * @param N_val The N parameter for H1/N (e.g., 3 for H1/3, 10 for H1/10). Must be strictly greater than 1.
 * @param H1 The scale parameter of the first Weibull distribution. Must be positive.
 * @param H2 The scale parameter of the second Weibull distribution. Must be positive.
 * @param k1 The exponent (shape parameter) of the first Weibull distribution. Must be positive.
 * @param k2 The exponent (shape parameter) of the second Weibull distribution. Must be positive.
 * @param Htr The transitional wave height.
 * @return The calculated value of H1/N.
 * @throws `std::invalid_argument` if input parameters are invalid (e.g., H1/H2 <= 0, k1/k2 <= 0, N_val <= 1).
 */
double calculate_H1N(double N_val, double H1, double H2, double k1, double k2, double Htr) {
    // Input validation: Ensures all parameters are within their mathematically valid ranges.
    if (H1 <= 0.0 || H2 <= 0.0) {
        throw std::invalid_argument("calculate_H1N: H1 and H2 must be positive.");
    }
    if (k1 <= 0.0 || k2 <= 0.0) {
        throw std::invalid_argument("calculate_H1N: k1 and k2 must be positive.");
    }
    if (N_val <= 1.0) {
        throw std::invalid_argument("calculate_H1N: N_val must be greater than 1.");
    }

    // First, determine HN (the wave height with 1/N exceedance probability).
    // This is crucial because the integration limits and the specific formula for H1/N
    // depend on where HN falls relative to Htr.
    double H_N_val = calculate_HN(N_val, H1, H2, k1, k2, Htr);

    // Case 1: H_N_val is smaller than Htr.
    // This implies that the integration for H1/N spans both parts of the Composite Weibull distribution.
    double term1_a = 1.0 / k1 + 1.0;
    double term1_x_ln_Nval = std::log(N_val);
    double term1_x_HtrH1 = std::pow(Htr / H1, k1);

    double term2_a = 1.0 / k2 + 1.0;
    double term2_x_HtrH2 = std::pow(Htr / H2, k2);

    if (H_N_val < Htr - EPSILON) {
        // Contribution from the first Weibull distribution.
        double gamma_term1_part1 = incomplete_gamma_upper(term1_a, term1_x_ln_Nval);
        double gamma_term1_part2 = incomplete_gamma_upper(term1_a, term1_x_HtrH1);
        double gamma_term1 = gamma_term1_part1 - gamma_term1_part2;

        // Contribution from the second Weibull distribution.
        double gamma_term2 = incomplete_gamma_upper(term2_a, term2_x_HtrH2);

        return N_val * H1 * gamma_term1 + N_val * H2 * gamma_term2;
    } else {
        // Case 2: H_N_val is greater than or equal to Htr.
        // This means the integration for H1/N only involves the second part of the Composite Weibull distribution.
        return N_val * H2 * incomplete_gamma_upper(term2_a, term1_x_ln_Nval);
    }
}

// --- Functions for the System of Non-Linear Equations (Newton-Raphson Matrix Method) ---

/**
 * @brief Defines the first non-linear equation F1(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
 *
 * This equation represents the normalized Hrms constraint for the Composite Weibull distribution.
 * It states that the overall normalized Hrms of the Composite Weibull distribution must precisely equal 1.
 * The equation directly uses the unnormalized lower incomplete gamma function, γ(a,x),
 * and the unnormalized upper incomplete gamma function, Γ(a,x).
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms The normalized scale parameter of the second Weibull distribution.
 * @param Htr_Hrms The normalized transitional wave height (constant for a given solve).
 * @return The value of the first function, which should be driven to zero.
 */
double F1(double H1_Hrms, double H2_Hrms, double Htr_Hrms) {
    // Input validation for H1_Hrms and H2_Hrms to prevent issues like log(0) or sqrt(negative).
    // While the solver should ideally keep values positive, these checks add robustness.
    if (H1_Hrms <= 0.0 || H2_Hrms <= 0.0) {
        // Return a large value to push the solver away from invalid regions.
        return std::numeric_limits<double>::max();
    }

    double arg1 = std::pow(Htr_Hrms / H1_Hrms, k1);
    double arg2 = std::pow(Htr_Hrms / H2_Hrms, k2);

    // Calculate terms using unnormalized incomplete gamma functions.
    double term1 = H1_Hrms * H1_Hrms * incomplete_gamma_lower(2.0 / k1 + 1.0, arg1);
    double term2 = H2_Hrms * H2_Hrms * incomplete_gamma_upper(2.0 / k2 + 1.0, arg2);

    // Ensure the argument to sqrt is non-negative, though theoretically it should be.
    double sum_terms = term1 + term2;
    if (sum_terms < 0.0) sum_terms = 0.0; // Clamp to zero to avoid NaN from sqrt of negative.

    return std::sqrt(sum_terms) - 1.0;
}

/**
 * @brief Defines the second non-linear equation F2(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
 *
 * This equation represents the continuity condition between the two Weibull distributions
 * at the transitional wave height Htr: `(Htr/H1)^k1 = (Htr/H2)^k2`.
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution.
 * @param H2_Hrms The normalized scale parameter of the second Weibull distribution.
 * @param Htr_Hrms The normalized transitional wave height (constant for a given solve).
 * @return The value of the second function, which should be driven to zero.
 */
double F2(double H1_Hrms, double H2_Hrms, double Htr_Hrms) {
    // Input validation for H1_Hrms and H2_Hrms to prevent division by zero or log(0).
    if (H1_Hrms <= 0.0 || H2_Hrms <= 0.0) {
        // Return a large value to push the solver away from invalid regions.
        return std::numeric_limits<double>::max();
    }
    return std::pow(Htr_Hrms / H1_Hrms, k1) - std::pow(Htr_Hrms / H2_Hrms, k2);
}

/**
 * @brief Solves a 2x2 linear system Ax = b for x using Cramer's rule.
 *
 * This function is a helper for the Newton-Raphson method for systems.
 * It takes the Jacobian matrix elements (J11, J12, J21, J22) and the negative
 * function values (-F1, -F2) as the right-hand side, and computes the updates
 * (dx1, dx2) for H1_Hrms and H2_Hrms.
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

    // Check for a singular or nearly singular Jacobian matrix.
    if (std::abs(determinant) < std::numeric_limits<double>::epsilon() * 100) {
        // If determinant is too small, the matrix is singular or ill-conditioned.
        // In this case, we cannot reliably solve the system.
        // Set dx1 and dx2 to zero to prevent large, unstable steps.
        dx1 = 0.0;
        dx2 = 0.0;
        return;
    }

    // Apply Cramer's rule:
    dx1 = (b1 * J22 - b2 * J12) / determinant;
    dx2 = (J11 * b2 - J21 * b1) / determinant;
}

/**
 * @brief Provides initial guesses for H1_Hrms and H2_Hrms based on Htr_Hrms.
 *
 * A good initial guess is crucial for the efficiency and robustness of the Newton-Raphson method.
 * This function uses an empirical regression for H1_Hrms, and then derives H2_Hrms from the
 * continuity condition, assuming an initial relationship.
 *
 * @param Htr_Hrms The normalized transitional wave height.
 * @param H1_initial Output: The initial guess for H1_Hrms.
 * @param H2_initial Output: The initial guess for H2_Hrms.
 */
void get_initial_guesses(double Htr_Hrms, double &H1_initial, double &H2_initial) {
    // Empirical regression for H1/Hrms.
    H1_initial = 0.9552427998926 / (1.0 - 0.992405988921401 * exp(-1.42537392576977 * Htr_Hrms));

    // Empirical regression for H2_initial
    H2_initial = 1.054085273232950 + 0.9369023639428842 * pow(Htr_Hrms, 2.980718327103574) / (pow(2.549022900471753, 2.980718327103574) + pow(Htr_Hrms, 2.980718327103574));

    // Ensure initial guesses are positive, as physical parameters cannot be zero or negative.
    if (H1_initial <= 0.0) H1_initial = std::numeric_limits<double>::min(); // Small positive value
    if (H2_initial <= 0.0) H2_initial = std::numeric_limits<double>::min(); // Small positive value
}

/**
 * @brief Solves for H1_Hrms and H2_Hrms simultaneously using the Newton-Raphson method for systems.
 *
 * This function implements the multi-dimensional Newton-Raphson algorithm to find the roots
 * of the system of non-linear equations F1 and F2. It iteratively refines the guesses
 * for H1_Hrms and H2_Hrms until the functions F1 and F2 are sufficiently close to zero.
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
    // Get initial guesses for H1_Hrms and H2_Hrms.
    get_initial_guesses(Htr_Hrms, H1_Hrms, H2_Hrms);

    for (int iter = 0; iter < maxit; ++iter) {
        // Evaluate the functions at the current guesses.
        double f1_val = F1(H1_Hrms, H2_Hrms, Htr_Hrms);
        double f2_val = F2(H1_Hrms, H2_Hrms, Htr_Hrms);

        // Check for convergence. If both function values are close to zero, we've converged.
        if (std::abs(f1_val) < tol && std::abs(f2_val) < tol) {
            return true;
        }

        // Calculate the Jacobian matrix elements using central finite differences.
        // J11 = dF1/dH1
        double J11 = (F1(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F1(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        // J12 = dF1/dH2
        double J12 = (F1(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F1(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        // J21 = dF2/dH1
        double J21 = (F2(H1_Hrms + JACOBIAN_DX, H2_Hrms, Htr_Hrms) - F2(H1_Hrms - JACOBIAN_DX, H2_Hrms, Htr_Hrms)) / (2.0 * JACOBIAN_DX);
        // J22 = dF2/dH2
        double J22 = (F2(H1_Hrms, H2_Hrms + JACOBIAN_DX, Htr_Hrms) - F2(H1_Hrms, H2_Hrms - JACOBIAN_DX, Htr_Hrms)) / (2.0 * JACOBIAN_DX);

        // Solve the linear system J * dx = -F for dx.
        // Here, dx = [dH1, dH2]^T and F = [f1_val, f2_val]^T.
        double dH1, dH2;
        solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val, dH1, dH2);

        // Update the guesses.
        H1_Hrms += dH1;
        H2_Hrms += dH2;

        // Ensure H1_Hrms and H2_Hrms remain positive. If they become non-positive,
        // clamp them to a small positive value to prevent mathematical errors in subsequent iterations.
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
static std::wstring buildReport(double Hm0, double m0, double d, double slopeM)
{
    std::wstringstream ss;
    ss << std::fixed << std::setprecision(4);

    // Input validation for physical parameters
    if (Hm0 <= 0.0 || m0 <= 0.0 || d <= 0.0)
    {
        ss << L"ERROR: Hm0, m0, and d must be positive.\n";
        return ss.str();
    }
    if (slopeM <= 0.0)
    {
        ss << L"ERROR: Beach slope (m) must be positive.\n";
        return ss.str();
    }

    // Mean square wave height: Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0)
    double Hrms = (2.69 + 3.24 * std::sqrt(m0) / d) * std::sqrt(m0);

    // Compute the actual tangent of the beach slope: tan(alpha) = 1.0 / slopeM.
    double tanAlpha = 1.0 / slopeM;
    // Dimensional transitional wave height: Htr = (0.35 + 5.8*(1/m)) * d.
    double Htr_dim = (0.35 + 5.8 * tanAlpha) * d;
    // Dimensionless transitional parameter: H̃_tr = Htr / Hrms.
    double Htr_tilde = (Hrms > 0.0) ? (Htr_dim / Hrms) : 0.0;

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
    ss << L"m0 (m^2)        : " << m0 << L"\n";
    ss << L"d (m)           : " << d << L"\n";
    ss << L"Beach slope (m) : " << slopeM << L"   (tan(alpha) = " << tanAlpha << L")\n\n";

    ss << L"===========================\n";
    ss << L"   CALCULATED PARAMETERS\n";
    ss << L"===========================\n";
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
    double m0 = 0.0;
    double d = 0.0;
    double slopeM = 0.0;

    if (argc >= 5) {
        try {
            Hm0 = std::stod(argv[1]);
            m0 = std::stod(argv[2]);
            d = std::stod(argv[3]);
            slopeM = std::stod(argv[4]);
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

        std::wcout << L"Enter free-surface variance m0 (m^2): ";
        std::wcin >> m0;
        
        std::wcout << L"Enter water depth d (m): ";
        std::wcin >> d;
        
        std::wcout << L"Enter beach slope (m): ";
        std::wcin >> slopeM;
    }

    std::wstring report = buildReport(Hm0, m0, d, slopeM);
    
    writeReportToFile(report);

    std::wcout << L"\n" << report << std::endl;

    return 0;
}