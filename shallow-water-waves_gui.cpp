/***********************************************************************
 * Program: shallow-water-waves_gui.cpp
 *
 * Detailed Description:
 * This program computes local shallow-foreshore wave-height distribution
 * parameters using a model based on the Composed Weibull distribution.
 *
 * The program is implemented as a Windows GUI application using only the
 * native Win32 API and standard C++ (with OpenMP directives for parallelism).
 * It allows the user to input three key parameters:
 *
 * 1. Hm0 (in meters)          - The local significant spectral wave height.
 * 2. d (in meters)            - The local water depth.
 * 3. Beach slope (1:m)        - The beach slope expressed as “1:m”.
 * (For example, enter 20 for a slope of 1/20 = 0.05.)
 *
 * Based on these inputs, the program computes:
 *
 * - Free-surface variance: m0 = (Hm0 / 4)²
 *
 * - Mean square wave height:
 * Hrms = (3.00 + 3.50*sqrt(m0)/d)*sqrt(m0)
 * (Empirical coefficients have been slightly increased relative to the
 * original deep-water formula to better capture the shallow-water
 * distribution of extreme waves.)
 *
 * - A dimensional transitional wave height:
 * Htr = (0.35 + 5.8*(1/m)) * d.
 * (For example, if m = 20 then tan(alpha)=1/20=0.05 and
 * Htr = (0.35 + 5.8*0.05)*d = 0.64*d.)
 *
 * - The dimensionless transitional parameter:
 * H̃_tr = Htr / Hrms.
 * If H̃_tr is above 3.5 then it is set to 3.5 and Htr is recalculated as
 * Htr = 3.5 * Hrms.
 *
 * - The dimensionless wave-height ratios (Hᵢ/Hrms) are calculated
 * by solving a nonlinear equation derived from the Composite Weibull
 * distribution, ensuring the normalized Hrms of the distribution equals one.
 * This involves using numerical methods such as Newton-Raphson with a bisection fallback for root-finding
 * and functions for incomplete gamma calculations. These ratios are then
 * converted to dimensional quantities (in meters) by multiplying with Hrms.
 *
 * - A detailed report is then generated (and written to "report.txt") with
 * the input parameters, intermediate values, calculated ratios and computed
 * dimensional wave heights, as well as diagnostic ratios.
 *
 * Compilation Instructions (example using g++ on Windows with OpenMP):
 *
 * g++ -O3 -Wall -municode shallow-water-waves_gui.cpp -o shallow-water-waves_gui ^
 * -mwindows -static -static-libgcc -static-libstdc++ -fopenmp
 *
 ***********************************************************************/

#ifdef _OPENMP
#include <omp.h>
#endif

#define _USE_MATH_DEFINES // Required on some systems (e.g., Windows with MSVC) for M_PI

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
#include <limits>       // For numeric_limits (e.g., std::numeric_limits<double>::min())
#include <algorithm>    // For std::swap
#include <iterator>     // For std::size (C++17)

// Global parameters for the Composite Weibull distribution.
const double k1 = 2.0;    // Exponent for the first part of the Composite Weibull distribution (Rayleigh-shaped)
const double k2 = 3.6;    // Exponent for the second part of the Composite Weibull distribution

// Precision for Numerical Solver Convergence:
const double EPSILON = 1e-9; // A small value (10^-9) indicating the maximum allowable error or difference

/**
 * @brief Computes the normalized lower incomplete gamma function P(a, x) = γ(a, x)/Γ(a).
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
 * @return The computed value of P(a, x). Returns `nan("")` (Not-a-Number) if input parameters are invalid
 * or if the series/continued fraction fails to converge within the maximum iterations.
 */
double incomplete_gamma_p(double a, double x)
{
    const int MAXIT = 500;      // Maximum number of iterations allowed for either the series expansion or the continued fraction.
    const double LOCAL_EPS = 1e-14;   // A very small tolerance value (10^-14) used for checking convergence
                                      // within the gamma function calculations, ensuring high numerical precision.

    // Input validation: Essential for preventing mathematical errors and ensuring function robustness.
    if(a <= 0.0 || x < 0.0)
        return std::nan(""); // If 'a' is not positive or 'x' is negative, return NaN as the result is undefined.
    if(x == 0.0)
        return 0.0; // By definition, P(a, 0) is 0 for any positive 'a'.

    // Compute the natural logarithm of the complete gamma function, ln(Γ(a)).
    // `std::lgamma` is used instead of `std::tgamma` followed by `log` for improved numerical stability,
    // especially when 'a' is large, as `lgamma` avoids potential overflow/underflow issues with very large/small gamma values.
    double gln = std::lgamma(a); 
    
    // Conditional logic to select the appropriate numerical method (series or continued fraction).
    // This choice is based on the relationship between 'x' and 'a', a common heuristic for optimal performance and accuracy.
    if(x < a + 1.0) {  // Use series expansion for x < a+1 (often referred to as the "gamma series" method).
        double ap = a;          // `ap` is a mutable copy of 'a', used to increment in the series terms.
        double sum = 1.0 / a;   // Initialize the sum with the first term of the series expansion.
        double del = sum;       // `del` stores the value of the current term being added to the sum.
        for (int n_iter = 1; n_iter <= MAXIT; ++n_iter) { // Iterate up to MAXIT to compute series terms.
            ap += 1.0;          // Increment 'a' for the next term in the series (a+1, a+2, ...).
            del *= x / ap;      // Calculate the next term: `del_n = del_{n-1} * (x / (a + n))`.
            sum += del;         // Add the newly calculated term to the cumulative sum.
            // Check for convergence: If the absolute value of the current term (`del`) is extremely small
            // relative to the absolute value of the cumulative sum (`sum`), then the series has converged.
            if(std::fabs(del) < std::fabs(sum) * LOCAL_EPS)
                // Return the final result. The series sum is multiplied by `exp(-x + a*log(x) - gln)`
                // to account for the common factors factored out for numerical stability in the series derivation.
                return sum * std::exp(-x + a * std::log(x) - gln);
        }
        return std::nan(""); // If the loop finishes without `del` becoming sufficiently small, it means convergence failed.
    } else {  // Use continued fraction for x >= a+1 (often referred to as the "gamma fraction" method).
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
            if(std::fabs(d) < std::numeric_limits<double>::min())
                d = std::numeric_limits<double>::min();
            c = b + an / c; // Apply the recurrence relation for C_n (numerator part of the fraction).
            // Safeguard against division by zero or extremely small numbers for 'c'.
            if(std::fabs(c) < std::numeric_limits<double>::min())
                c = std::numeric_limits<double>::min();
            d = 1.0 / d;    // Invert 'd' for the next step.
            double del = d * c; // Calculate the current term (delta_n) of the continued fraction.
            h *= del;       // Multiply `h` by the current term to update the accumulated fraction value.
            // Check for convergence: If the current term `del` is very close to 1.0, the continued fraction has converged.
            if(std::fabs(del - 1.0) < LOCAL_EPS)
                break; // Exit the loop if converged.
        }
        // Calculate the exponential and power part of the incomplete gamma function: `exp(-x + a*log(x) - gln)`.
        double lnPart = -x + a * std::log(x) - gln;
        // `Qval` represents Γ(a,x)/Γ(a) (the normalized upper incomplete gamma function)
        // as computed by the continued fraction method.
        double Qval = std::exp(lnPart) * h;
        // The normalized lower incomplete gamma function P(a,x) is 1 - Q(a,x).
        // It's clamped to 0.0 to prevent very small negative floating-point results due to precision issues.
        double Pval = 1.0 - Qval;
        return (Pval < 0.0) ? 0.0 : Pval; 
    }
}

/**
 * @brief Computes the normalized upper incomplete gamma function Q(a, x) = Γ(a, x)/Γ(a).
 *
 * This function is a simple wrapper that leverages the `incomplete_gamma_p` function.
 * By definition, Q(a, x) is the complement of P(a, x), i.e., 1 - P(a, x).
 * It is declared `inline` to suggest to the compiler that its code should be inserted directly
 * at the call site, potentially optimizing performance by avoiding function call overhead.
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
 * This function calculates the unnormalized form of the upper incomplete gamma function,
 * which is used in the integral formulations for wave height statistics. It is derived
 * from the normalized upper incomplete gamma function Q(a, x) and the complete gamma function Γ(a).
 *
 * @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
 * @param x The 'x' parameter (lower integration limit) of the incomplete gamma function. Must be non-negative.
 * @return The computed value of the unnormalized upper incomplete gamma function Γ(a, x).
 * @throws `std::invalid_argument` if input parameters ('a' or 'x') are outside their valid ranges.
 */
double incomplete_gamma(double a, double x) {
    // Input validation: Ensures 'a' is positive and 'x' is non-negative, as required by the gamma function definitions.
    if (a <= 0.0) {
        throw std::invalid_argument("incomplete_gamma: 'a' must be positive.");
    }
    if (x < 0.0) {
        throw std::invalid_argument("incomplete_gamma: 'x' must be non-negative.");
    }
    // Calculation: Γ(a,x) = Q(a,x) * Γ(a).
    // `std::tgamma` computes the complete gamma function Γ(a).
    return incomplete_gamma_q(a, x) * std::tgamma(a); 
}

/**
 * @brief Calculates HN (wave height with 1/N exceedance probability) for the Composite Weibull distribution.
 *
 * This function determines the specific wave height (H) such that the probability of a wave
 * exceeding this height is 1/N. This is a key statistical measure. The calculation depends
 * on whether the target wave height falls into the first or second part of the Composite Weibull
 * distribution, which is separated by the transitional wave height (Htr).
 * The logic directly implements the formulas from Appendix A.2.1 of Groenendijk (1998).
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
    // This uses Equation A.5 from Groenendijk (1998): `H_N,1 = H1 * (ln(N))^(1/k1)`.
    double HN_candidate1 = H1 * std::pow(std::log(N), 1.0 / k1);

    // Decision point: Check if the `HN_candidate1` is indeed less than the transitional wave height (Htr).
    // A small `EPSILON` is subtracted from `Htr` to account for floating-point inaccuracies when comparing.
    if (HN_candidate1 < Htr - EPSILON) {
        // If true, it means the wave height for the given exceedance probability is governed by the first Weibull part.
        return HN_candidate1;
    } else {
        // If `HN_candidate1` is not less than `Htr`, it implies that the wave height for the given
        // exceedance probability is determined by the *second* part of the distribution (H > Htr).
        // This uses Equation A.7 from Groenendijk (1998): `H_N = H2 * (ln(N))^(1/k2)`.
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
 * The implementation follows the detailed derivations in Appendix A.2.2 of Groenendijk (1998).
 *
 * @param N_val The N parameter for H1/N (e.g., 3.0 for H1/3, 10.0 for H1/10). Must be strictly greater than 1.
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
    // The formula used here corresponds to Equation A.15 in Groenendijk (1998).
    if (H_N_val < Htr - EPSILON) {
        // Term 1: Contribution from the first Weibull distribution (F1(H)).
        // The 'a' parameter for the incomplete gamma function in this context is (1/k1 + 1).
        double term1_a = 1.0 / k1 + 1.0;
        // The 'x' parameter for the first incomplete gamma term is `ln(N_val)`.
        double term1_x_ln_Nval = std::log(N_val);
        // The 'x' parameter for the second incomplete gamma term (related to Htr) is `(Htr/H1)^k1`.
        double term1_x_HtrH1 = std::pow(Htr / H1, k1);

        // Calculate the two incomplete gamma terms for the first part of the integral.
        // This represents `Γ[a, ln(N)] - Γ[a, (Htr/H1)^k1]`.
        double gamma_term1_part1 = incomplete_gamma(term1_a, term1_x_ln_Nval);
        double gamma_term1_part2 = incomplete_gamma(term1_a, term1_x_HtrH1);
        double gamma_term1 = gamma_term1_part1 - gamma_term1_part2;

        // Term 2: Contribution from the second Weibull distribution (F2(H)).
        // The 'a' parameter for the incomplete gamma function is (1/k2 + 1).
        double term2_a = 1.0 / k2 + 1.0;
        // The 'x' parameter for this term is `(Htr/H2)^k2`.
        double term2_x_HtrH2 = std::pow(Htr / H2, k2);

        // Calculate the incomplete gamma term for the second part of the integral.
        // This represents `Γ[a, (Htr/H2)^k2]`.
        double gamma_term2 = incomplete_gamma(term2_a, term2_x_HtrH2);

        // Combine the terms as per Equation A.15: `N_val * H1 * (gamma_term1) + N_val * H2 * (gamma_term2)`.
        return N_val * H1 * gamma_term1 + N_val * H2 * gamma_term2;
    } else {
        // Case 2: H_N_val is greater than or equal to Htr.
        // This means the integration for H1/N only involves the second part of the Composite Weibull distribution.
        // The formula used here corresponds to Equation A.20 in Groenendijk (1998).
        double term_a = 1.0 / k2 + 1.0; // The 'a' parameter for the incomplete gamma function.
        double term_x = std::log(N_val);     // The 'x' parameter is `ln(N_val)`.
        // Calculate `N_val * H2 * Γ[a, ln(N_val)]`.
        return N_val * H2 * incomplete_gamma(term_a, term_x);
    }
}

/**
 * @brief Computes the residual function f(H1_Hrms) for the root-finding problem.
 *
 * The objective of the numerical solver is to find the value of H1_Hrms that makes this
 * residual function equal to zero. This function embodies the core constraint of the
 * Composite Weibull distribution: that its overall root-mean-square wave height (Hrms),
 * when normalized by itself, must equal 1.0.
 *
 * The function is formulated as `sqrt(sum_of_terms) - 1.0`. Finding the root means
 * `sqrt(sum_of_terms) = 1.0`, which implies `sum_of_terms = 1.0`.
 * This corresponds to Equation 1.7 in Groenendijk & Van Gent (1998) or Equation A.26
 * in Appendix A.2.3 of Groenendijk (1998) (the normalized version of Hrms).
 *
 * @param H1_Hrms The normalized scale parameter of the first Weibull distribution (H1/Hrms).
 * This is the variable for which the root is being sought.
 * @param Htr_Hrms The normalized transitional wave height (Htr/Hrms). This value is constant
 * for a given iteration of the main program loop.
 * @return The calculated value of the residual function.
 */
double residual(double H1_Hrms, double Htr_Hrms)
{
    // Calculate H2_Hrms (normalized) based on the continuity condition between the two Weibull distributions.
    // The continuity condition states that F1(Htr) = F2(Htr), which simplifies to (Htr/H1)^k1 = (Htr/H2)^k2.
    // Rearranging this equation to solve for H2_Hrms: H2_Hrms = Htr_Hrms * (H1_Hrms / Htr_Hrms)^(k1 / k2).
    double H2_Hrms = Htr_Hrms * std::pow(H1_Hrms / Htr_Hrms, k1 / k2);
    
    // Calculate the arguments for the incomplete gamma functions. These arguments are derived
    // from the integration limits (0 to Htr_Hrms for the first part, Htr_Hrms to infinity for the second)
    // and the specific form of the Weibull probability density function when calculating moments.
    double arg1 = std::pow(Htr_Hrms / H1_Hrms, k1); // Argument for the incomplete gamma function related to the first Weibull part.
    double arg2 = std::pow(Htr_Hrms / H2_Hrms, k2); // Argument for the incomplete gamma function related to the second Weibull part.
    
    // Compute the normalized lower incomplete gamma function P(a,x) for the first term.
    // The 'a' parameter for the Hrms calculation integral (integral of H^2 * f(H)) is (2/k + 1).
    // P(a,x) is used here because the integral for the first part of the Hrms calculation goes from 0 to Htr_Hrms.
    double g1_val = incomplete_gamma_p(2.0 / k1 + 1.0, arg1); 
    
    // Compute the normalized upper incomplete gamma function Q(a,x) for the second term.
    // Q(a,x) is used here because the integral for the second part of the Hrms calculation goes from Htr_Hrms to infinity.
    double g2_val = incomplete_gamma_q(2.0 / k2 + 1.0, arg2); 
    
    // Calculate the sum of terms inside the square root from the Hrms equation.
    // This effectively represents (Hrms_calculated / Hrms_actual)^2, where Hrms_actual is the true Hrms
    // of the distribution, which is normalized to 1.0 in this context.
    double sum_terms = H1_Hrms * H1_Hrms * g1_val + H2_Hrms * H2_Hrms * g2_val;
    
    // The residual is the difference between the calculated normalized Hrms and 1.0.
    // The goal of the solver is to make this value as close to zero as possible.
    return std::sqrt(sum_terms) - 1.0;
}

/**
 * @brief Computes the derivative of the residual function using a central finite difference approximation.
 *
 * The derivative of the residual function is a crucial input for the Newton-Raphson method.
 * Newton-Raphson uses the function's value and its slope (derivative) at the current guess
 * to estimate a better next guess for the root. A central finite difference is preferred
 * over a forward or backward difference for its improved accuracy.
 *
 * @param H1_Hrms The normalized H1 parameter at which the derivative is to be calculated.
 * @param Htr_Hrms The normalized transitional wave height (constant for this calculation).
 * @return The approximate value of the derivative of the residual function at H1_Hrms.
 */
double derivative_of_residual(double H1_Hrms, double Htr_Hrms) {
    const double dx = 1e-6; // A small step size used for the finite difference approximation.
                            // This value balances accuracy (smaller dx is more accurate) with
                            // numerical stability (too small dx can lead to precision issues).
    // Apply the central finite difference formula: f'(x) ≈ (f(x + dx) - f(x - dx)) / (2 * dx).
    return (residual(H1_Hrms + dx, Htr_Hrms) - residual(H1_Hrms - dx, Htr_Hrms)) / (2.0 * dx);
}

/**
 * @brief Solves for the root of the function f(H1_Hrms)=0 using the Newton-Raphson method.
 *
 * The Newton-Raphson method is an efficient iterative numerical technique for finding
 * the roots (or zeros) of a real-valued function. It starts with an initial guess and
 * iteratively refines it by moving along the tangent line of the function at the current guess.
 *
 * @param Htr_Hrms The normalized transitional wave height, which remains constant during a single solve.
 * @param initial_guess The starting point for the iterative process. A good initial guess can
 * significantly improve convergence speed and robustness.
 * @param tol The desired tolerance for convergence. The iteration stops when the absolute
 * difference between successive guesses is less than this tolerance.
 * @param maxit The maximum number of iterations allowed. This prevents infinite loops
 * in cases where the solver might not converge or converges very slowly.
 * @return The estimated root of the function (H1_Hrms) if convergence is achieved within `maxit`.
 * If the derivative becomes too small (indicating a flat region or local extremum),
 * it returns the current best guess to prevent division by zero.
 */
double newtonRaphsonSolver(double Htr_Hrms, double initial_guess, double tol, int maxit) {
    double x_old = initial_guess; // Stores the current approximation of the root.
    double x_new = initial_guess; // Stores the next, refined approximation of the root.

    for (int iter = 0; iter < maxit; ++iter) { // Loop for a maximum of `maxit` iterations.
        double fx = residual(x_old, Htr_Hrms); // Calculate the value of the function (residual) at the current guess.
        double f_prime_x = derivative_of_residual(x_old, Htr_Hrms); // Calculate the derivative of the function at the current guess.

        // Critical check: If the absolute value of the derivative is extremely small, it indicates
        // that the function is very flat around `x_old`, or `x_old` is near a local extremum.
        // In such cases, the Newton-Raphson formula (division by `f_prime_x`) becomes unstable.
        if (std::fabs(f_prime_x) < std::numeric_limits<double>::min()) { 
            return x_old; // Return the current best guess as the solver cannot proceed reliably.
        }

        // Newton-Raphson iteration formula: x_{n+1} = x_n - f(x_n) / f'(x_n).
        x_new = x_old - fx / f_prime_x;

        // Check for convergence: If the absolute difference between the new and old guesses
        // is less than the specified tolerance, the solution has converged.
        if (std::fabs(x_new - x_old) < tol) {
            return x_new; // Return the converged root.
        }
        x_old = x_new; // Update the old guess to the new guess for the next iteration.
    }
    return x_new; // If the loop completes without converging (i.e., `maxit` iterations are reached),
                  // return the last computed `x_new` as the best approximation found.
}

/**
 * @brief Provides an initial guess for H1_Hrms based on Htr_Hrms using an empirical regression formula.
 *
 * Supplying a good initial guess is crucial for the efficiency and reliability of iterative
 * numerical solvers like Newton-Raphson. This empirical formula is derived from a pre-analysis
 * or regression fit of the relationship between Htr_Hrms and H1_Hrms, allowing the solver
 * to start very close to the actual root.
 *
 * @param Htr_Hrms The normalized transitional wave height, used as input for the regression.
 * @return An empirically derived initial guess for H1_Hrms.
 */
double regression_initial_guess(double Htr_Hrms)
{
    double logHtr = std::log(Htr_Hrms); // Calculate the natural logarithm of Htr_Hrms.
                                   // Logarithmic transformations are common in empirical fits.
    // This is the empirical regression formula. The coefficients (-1.705, -2.329, -0.313)
    // are specific to this model and were likely determined through curve fitting to data.
    return 1.0 + std::exp(-1.705 - 2.329 * logHtr - 0.313 * logHtr * logHtr);
}

/**
 * @brief Computes H1_Hrms with high accuracy, employing a robust numerical strategy.
 *
 * This function orchestrates the root-finding process for H1_Hrms. It first attempts to use
 * the fast and efficient Newton-Raphson method. If Newton-Raphson fails to converge
 * (e.g., due to a poor initial guess or problematic function behavior), it falls back to
 * the bisection method, which is slower but guaranteed to converge if a root is successfully
 * bracketed. This hybrid approach significantly enhances the overall robustness of the solver.
 *
 * @param Htr_Hrms The normalized transitional wave height, which is a fixed input for this solve.
 * @param H1_Hrms A reference to a double variable where the computed H1_Hrms value will be stored.
 * @return `true` if the solver successfully finds a root within the specified tolerance, `false` otherwise.
 */
bool highAccuracySolver(double Htr_Hrms, double &H1_Hrms)
{
    // Obtain an initial guess for H1_Hrms using the empirical regression model.
    double initial_guess = regression_initial_guess(Htr_Hrms);
    const double tol = EPSILON; // The very high desired tolerance for the root, indicating extreme precision.
    const int maxit = 200;    // The maximum number of iterations allowed for both Newton-Raphson and bisection.

    // Attempt to solve for H1_Hrms using the Newton-Raphson method.
    H1_Hrms = newtonRaphsonSolver(Htr_Hrms, initial_guess, tol, maxit);

    // After attempting Newton-Raphson, check if the residual at the found H1_Hrms is sufficiently close to zero.
    // A slightly relaxed tolerance (tol * 100) is used for this initial check to avoid overly strict failure.
    if (std::fabs(residual(H1_Hrms, Htr_Hrms)) > tol * 100) { 
        // If Newton-Raphson did not converge satisfactorily, initiate the bisection fallback.
        // Bisection is a robust but slower method that guarantees convergence if a root is bracketed.
        double a = 0.01, b = 50.0; // Define a wide initial bracket [a, b] for the possible range of H1_Hrms.
                                   // This bracket must contain the root.
        double fa = residual(a, Htr_Hrms); // Calculate the residual at the lower bound 'a'.
        double fb = residual(b, Htr_Hrms); // Calculate the residual at the upper bound 'b'.

        // Check if a root is successfully bracketed (i.e., the function values at the bounds have opposite signs).
        if (fa * fb < 0) { 
            for (int i = 0; i < maxit; ++i) { // Iterate for bisection up to `maxit`.
                double c = (a + b) / 2.0; // Calculate the midpoint of the current interval.
                double fc = residual(c, Htr_Hrms); // Calculate the residual at the midpoint.

                // Check for convergence: if the residual at 'c' is very close to zero, or if the interval
                // [a, b] has become extremely small, the root has been found.
                if (std::fabs(fc) < tol || std::fabs(b - a) < tol) {
                    H1_Hrms = c; // Assign the midpoint as the found root.
                    return true; // Indicate successful convergence.
                }

                // Narrow the bracket: If the product `fa * fc` is negative, the root is in [a, c].
                // Otherwise, the root is in [c, b].
                if (fa * fc < 0) { 
                    b = c; // The new upper bound is 'c'.
                } else { 
                    a = c;   // The new lower bound is 'c'.
                    fa = fc; // Update `fa` to `fc` for the next iteration (important for bisection logic).
                }
            }
        }
        // If neither Newton-Raphson nor the bisection fallback converges, print an error message
        // to the standard error stream, indicating the failure for the specific Htr_Hrms value.
        // In a GUI, this would ideally be a message box.
        // MessageBox(NULL, L"Solver: Failed to converge for Htr_Hrms after Newton-Raphson and fallback.", L"Error", MB_ICONERROR | MB_OK);
        return false; // Indicate that the solver failed to converge.
    }
    return true; // Indicate that Newton-Raphson successfully converged.
}


// Build the detailed report using the computed values.
static std::wstring buildReport(double Hm0, double d, double slopeM)
{
    std::wstringstream ss;
    ss << std::fixed << std::setprecision(4);

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

    // Compute free-surface variance
    double m0 = std::pow(Hm0 / 4.0, 2.0);

    // Mean square wave height (using formula that incorporates 'd' and 'm0')
    double Hrms = (3.00 + 3.50 * std::sqrt(m0) / d) * std::sqrt(m0);

    // Compute the actual tangent of the beach slope: tan(alpha) = 1.0 / slopeM.
    double tanAlpha = 1.0 / slopeM;
    // Dimensional transitional wave height (using formula that incorporates 'tanAlpha' and 'd')
    double Htr_dim = (0.35 + 5.8 * tanAlpha) * d;
    double Htr_tilde = (Hrms > 0.0) ? (Htr_dim / Hrms) : 0.0;

    // If dimensionless Htr exceeds 3.5, then adopt Htr_tilde = 3.5 and recalc Htr_dim.
    if (Htr_tilde > 3.5)
    {
        Htr_tilde = 3.5;
        Htr_dim = 3.5 * Hrms;
    }

    // Solve for H1_Hrms using the high accuracy solver
    double H1_Hrms;
    if (!highAccuracySolver(Htr_tilde, H1_Hrms)) {
        ss << L"ERROR: Solver failed to converge for Htr_tilde = " << Htr_tilde << L". Cannot generate full report.\n";
        return ss.str();
    }

    // Calculate H2_Hrms based on the derived H1_Hrms and the relationship
    double H2_Hrms = Htr_tilde * std::pow(H1_Hrms / Htr_tilde, k1 / k2);

    // Define the N values for which H(1/N)/Hrms quantiles are calculated.
    // H1/Hrms and H2/Hrms are directly computed, not H1/N quantiles.
    const int N_values_for_H1N[] = {3, 10, 50, 100, 250, 1000};
    const size_t num_H1N_values = std::size(N_values_for_H1N);
    std::vector<double> interp_Hrms(2 + num_H1N_values); // For H1/Hrms, H2/Hrms, and H1/N values

    // Store H1/Hrms and H2/Hrms directly
    interp_Hrms[0] = H1_Hrms;
    interp_Hrms[1] = H2_Hrms;

    // Calculate H1/N values using the new function
    for (size_t j = 0; j < num_H1N_values; ++j) {
        int N = N_values_for_H1N[j];
        // Note: H1_Hrms and H2_Hrms are already normalized by Hrms.
        // The calculate_H1N function expects H1 and H2 as scale parameters, so we pass the normalized ones.
        interp_Hrms[j + 2] = calculate_H1N(static_cast<double>(N), H1_Hrms, H2_Hrms, k1, k2, Htr_tilde);
    }

    // Convert dimensionless values to dimensional wave heights.
    std::vector<double> dimensional(interp_Hrms.size(), 0.0);
    for (size_t i = 0; i < interp_Hrms.size(); i++)
    {
        dimensional[i] = interp_Hrms[i] * Hrms;
    }

    // Calculate diagnostic ratios.
    // Indices for the specific ratios:
    // H1/3 is at index 2 (N=3)
    // H1/10 is at index 3 (N=10)
    // H1/50 is at index 4 (N=50)
    // H1/100 is at index 5 (N=100)
    // H1/250 is at index 6 (N=250)
    // H1/1000 is at index 7 (N=1000)
    double ratio_110_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[3] / interp_Hrms[2]) : 0.0;
    double ratio_150_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[4] / interp_Hrms[2]) : 0.0;
    double ratio_1100_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[5] / interp_Hrms[2]) : 0.0;
    double ratio_1250_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[6] / interp_Hrms[2]) : 0.0;
    double ratio_11000_13 = (interp_Hrms[2] != 0.0) ? (interp_Hrms[7] / interp_Hrms[2]) : 0.0;


    // Build the report string.
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
#define IDC_EDIT_SLOPE 105
#define IDC_BUTTON_COMPUTE 103
#define IDC_OUTPUT 104

HWND hEditHm0, hEditD, hEditSlope, hOutput;

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    static HFONT hMonoFont = NULL;

    switch (msg)
    {
    case WM_CREATE:
    {
        CreateWindow(L"STATIC", L"Hm0 (m):", WS_CHILD | WS_VISIBLE,
                     10, 10, 80, 20, hwnd, NULL, NULL, NULL);
        hEditHm0 = CreateWindow(L"EDIT", L"2.5", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                100, 10, 100, 20, hwnd, (HMENU)IDC_EDIT_HM0, NULL, NULL);

        CreateWindow(L"STATIC", L"d (m):", WS_CHILD | WS_VISIBLE,
                     10, 40, 80, 20, hwnd, NULL, NULL, NULL);
        hEditD = CreateWindow(L"EDIT", L"10.0", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                              100, 40, 100, 20, hwnd, (HMENU)IDC_EDIT_D, NULL, NULL);

        CreateWindow(L"STATIC", L"Beach slope m:", WS_CHILD | WS_VISIBLE,
                     10, 70, 120, 20, hwnd, NULL, NULL, NULL);
        // Default value is "20" (i.e. slope m = 20 gives tan(alpha)=1/20=0.05)
        hEditSlope = CreateWindow(L"EDIT", L"20", WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                                  140, 70, 60, 20, hwnd, (HMENU)IDC_EDIT_SLOPE, NULL, NULL);

        CreateWindow(L"BUTTON", L"Compute", WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
                     10, 100, 190, 30, hwnd, (HMENU)IDC_BUTTON_COMPUTE, NULL, NULL);

        hOutput = CreateWindow(L"EDIT", L"", WS_CHILD | WS_VISIBLE | WS_BORDER |
                               ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL | ES_READONLY,
                               10, 140, 760, 440, hwnd, (HMENU)IDC_OUTPUT, NULL, NULL);

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

            GetWindowText(hEditSlope, buffer, 63);
            double slopeM = _wtof(buffer);

            std::wstring report = buildReport(Hm0, d, slopeM);
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
