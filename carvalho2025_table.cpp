/***********************************************************************
 * Program: carvalho2025_table.cpp
 *
 * Detailed Description:
 * This program is designed to compute and tabulate normalized wave height parameters
 * for shallow-water environments, utilizing the Composite Weibull distribution model.
 * All wave heights are normalized by the root-mean-square wave height (Hrms) to provide
 * dimensionless ratios, which are fundamental for understanding wave statistics in complex
 * shallow foreshore conditions where wave breaking significantly alters traditional
 * wave height distributions (e.g., Rayleigh distribution).
 *
 * The core functionality of this program involves an iterative numerical solution:
 *
 * 1.  **Iterating through Normalized Transitional Wave Height (Htr/Hrms)**:
 * The program systematically varies the normalized transitional wave height (Htr_Hrms)
 * from 0.01 to 3.50 in increments of 0.01. The transitional wave height (Htr) is a
 * critical parameter in the Composite Weibull distribution, marking the point where
 * the wave height distribution transitions from a Rayleigh-like behavior (for smaller waves)
 * to a different, more complex behavior influenced by wave breaking (for larger waves).
 *
 * 2.  **Solving for Normalized Scale Parameter (H1/Hrms)**:
 * For each Htr_Hrms value, the program solves a non-linear equation to determine
 * H1_Hrms, which is the normalized scale parameter of the first part of the
 * Composite Weibull distribution. This non-linear equation is derived from the
 * fundamental constraint that the overall normalized Hrms of the Composite Weibull
 * distribution must precisely equal one. This ensures consistency and proper normalization
 * of the distribution. The solution employs a robust numerical root-finding method,
 * specifically Newton-Raphson with a bisection fallback for enhanced accuracy and stability.
 *
 * 3.  **Calculating Normalized Second Scale Parameter (H2/Hrms)**:
 * Once H1_Hrms is successfully determined, the program calculates H2_Hrms, the
 * normalized scale parameter for the second part of the Composite Weibull distribution.
 * This calculation is based on the continuity condition of the Composite Weibull
 * distribution at the transitional wave height (Htr). This condition ensures a smooth
 * transition between the two Weibull components of the distribution.
 *
 * 4.  **Computing Normalized Quantile Wave Heights (H(1/N)/Hrms)**:
 * Finally, for a predefined set of exceedance probabilities (represented by N values,
 * e.g., N=3 for H1/3, N=10 for H1/10, N=50 for H1/50, N=100 for H1/100, N=250 for H1/250,
 * and N=1000 for H1/1000), the program computes the corresponding normalized
 * quantile wave heights. These calculations utilize specific formulas (Equations A.15 or A.20)
 * from Groenendijk (1998), depending on the relationship between Htr_Hrms and the
 * specific quantile being calculated. These quantiles represent the mean of the highest
 * 1/N-part of the wave heights, providing crucial insights into the statistical properties
 * of extreme waves.
 *
 * Output:
 * The program generates a formatted text file named "carvalho2025_table.txt". This file
 * contains a comprehensive table of Htr_Hrms, H1_Hrms, H2_Hrms, and the various
 * H(1/N)/Hrms values, providing a valuable reference for wave height analysis in
 * shallow-water coastal engineering applications.
 *
 * Compilation Instructions:
 * To compile this program with high optimization levels and OpenMP support (though OpenMP
 * is not explicitly used in this single-threaded version, the flags are retained for consistency
 * with typical scientific computing setups), you can use a g++ command similar to this:
 *
 * g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic \
 * -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
 * -o carvalho2025_table.exe carvalho2025_table.cpp
 *
 ***********************************************************************/

// This preprocessor directive is crucial for ensuring that the M_PI constant (representing the mathematical constant Pi)
// is defined and available for use within the <cmath> header. This is particularly relevant when compiling
// on certain systems, such as Windows with Microsoft Visual C++ (MSVC), where it might not be defined by default.
#define _USE_MATH_DEFINES 

// Standard C++ Library Includes:
// These lines import necessary functionalities from the C++ Standard Library, providing tools
// for input/output, file handling, mathematical operations, data structures, and formatting.
#include <iostream>     // Provides standard input/output streams (e.g., `std::cout` for console output, `std::cerr` for error messages).
#include <fstream>      // Enables file stream operations, allowing the program to read from and write to files (e.g., `std::ofstream` for output files).
#include <sstream>      // Offers string stream capabilities (e.g., `std::ostringstream`), useful for in-memory string manipulation, though not directly used in the main logic of this specific version.
#include <vector>       // Provides the `std::vector` container, a dynamic array that can resize itself, used here to store calculated H1/N values efficiently.
#include <cmath>        // Contains a wide range of mathematical functions (e.g., `std::sqrt` for square root, `std::log` for natural logarithm, `std::pow` for exponentiation, `std::sin` for sine, `std::lgamma` for the natural logarithm of the gamma function, `std::tgamma` for the true gamma function, and `M_PI` if `_USE_MATH_DEFINES` is effective).
#include <limits>       // Provides `std::numeric_limits`, which allows querying properties of numeric types, such as the smallest representable double (`std::numeric_limits<double>::min()`), crucial for handling very small numbers and avoiding division by zero in numerical algorithms.
#include <iomanip>      // Offers manipulators for output formatting (e.g., `std::fixed` for fixed-point notation, `std::setprecision` to control decimal places, `std::setw` to set field width for formatted output).
#include <algorithm>    // Includes general-purpose algorithms (e.g., `std::max`, `std::swap`, `std::round`), although not explicitly called in the main execution path of this version, they are commonly useful.

// Using Namespace:
// This directive brings all identifiers (like `cout`, `endl`, `vector`, `sqrt`, etc.) from the
// `std` (standard) namespace directly into the current scope. This avoids the need to prefix
// these identifiers with `std::`, making the code more concise.
using namespace std;

// Global Parameters for the Composite Weibull Distribution:
// These constants define the fixed exponents for the two parts of the Composite Weibull distribution.
// They are declared globally as `const double` because their values are fundamental to the model
// and remain unchanged throughout the program's execution.
const double k1 = 2.0;    // Exponent for the first part of the Composite Weibull distribution.
                          // As per the theoretical foundation (Groenendijk, 1998, Section 2.1 and 3.3.2),
                          // the initial part of the wave height distribution in shallow water is assumed
                          // to be Rayleigh-shaped. A Rayleigh distribution is a special case of the Weibull
                          // distribution with an exponent (shape parameter) of 2.0.
const double k2 = 3.6;    // Exponent for the second part of the Composite Weibull distribution.
                          // This value was empirically determined through calibration and optimization
                          // processes described in Groenendijk (1998, Section 2.1) and Groenendijk & Van Gent (1998).
                          // It reflects the altered shape of the wave height distribution for larger waves
                          // due to depth-induced breaking.

// Precision for Numerical Solver Convergence:
// This constant defines the tolerance level used in numerical methods (like Newton-Raphson)
// to determine when an iterative solution has converged to an acceptable accuracy.
const double EPSILON = 1e-9; // A small value (10^-9) indicating the maximum allowable error or difference
                             // between successive iterations for a solution to be considered converged.

// Forward Declarations of Functions:
// These declarations inform the C++ compiler about the existence and signature (return type, name, and parameters)
// of functions that are defined later in the source code. This is necessary because some functions might call
// other functions that have not yet been fully defined at the point of the call. This prevents compilation errors.

// Declares a function to compute the normalized lower incomplete gamma function P(a, x) = γ(a, x)/Γ(a).
double incomplete_gamma_p(double a, double x);
// Declares a function to compute the normalized upper incomplete gamma function Q(a, x) = Γ(a, x)/Γ(a).
double incomplete_gamma_q(double a, double x);
// Declares a function to compute the unnormalized upper incomplete gamma function Γ(a, x).
double incomplete_gamma(double a, double x); 
// Declares a function to calculate HN (the wave height with a 1/N exceedance probability).
double calculate_HN(double N, double H1, double H2, double k1, double k2, double Htr);
// Declares a function to calculate H1/N (the mean of the highest 1/N-th part of wave heights).
double calculate_H1N(double N_val, double H1, double H2, double k1, double k2, double Htr);
// Declares a function that represents the residual (error) in the non-linear equation
// whose root needs to be found by the solver.
double residual(double H1_Hrms, double Htr_Hrms);
// Declares a function to compute the derivative of the residual function, essential for
// the Newton-Raphson numerical method.
double derivative_of_residual(double H1_Hrms, double Htr_Hrms); 
// Declares the Newton-Raphson solver function, an iterative method for finding roots of functions.
double newtonRaphsonSolver(double Htr_Hrms, double initial_guess, double tol, int maxit);
// Declares a function to provide an initial guess for the numerical solver,
// typically based on an empirical regression to speed up convergence.
double regression_initial_guess(double Htr_Hrms);
// Declares a wrapper function that orchestrates the high-accuracy root-finding process for H1/Hrms.
bool highAccuracySolver(double Htr_Hrms, double &H1_Hrms);

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
        return nan(""); // If 'a' is not positive or 'x' is negative, return NaN as the result is undefined.
    if(x == 0.0)
        return 0.0; // By definition, P(a, 0) is 0 for any positive 'a'.

    // Compute the natural logarithm of the complete gamma function, ln(Γ(a)).
    // `std::lgamma` is used instead of `std::tgamma` followed by `log` for improved numerical stability,
    // especially when 'a' is large, as `lgamma` avoids potential overflow/underflow issues with very large/small gamma values.
    double gln = lgamma(a); 
    
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
            if(abs(del) < abs(sum) * LOCAL_EPS)
                // Return the final result. The series sum is multiplied by `exp(-x + a*log(x) - gln)`
                // to account for the common factors factored out for numerical stability in the series derivation.
                return sum * exp(-x + a * log(x) - gln);
        }
        return nan(""); // If the loop finishes without `del` becoming sufficiently small, it means convergence failed.
    } else {  // Use continued fraction for x >= a+1 (often referred to as the "gamma fraction" method).
        double b = x + 1.0 - a; // Initialize the 'b' term (B_0) for the continued fraction expansion.
        // Initialize 'c' and 'd' to extremely large/small values. This is a common technique
        // to prevent division by zero or underflow during the initial steps of the continued fraction
        // algorithm, where denominators might otherwise be zero or very close to zero.
        double c = 1.0 / numeric_limits<double>::min(); 
        double d = 1.0 / b;
        double h = d; // `h` accumulates the value of the continued fraction.
        for (int i = 1; i <= MAXIT; ++i) { // Iterate up to MAXIT to compute continued fraction terms.
            double an = -1.0 * i * (i - a); // Calculate the 'a_n' term (numerator) for the current iteration.
            b += 2.0; // Update the 'b_n' term (denominator) for the current iteration.
            d = an * d + b; // Apply the recurrence relation for D_n (denominator part of the fraction).
            // Safeguard against division by zero or extremely small numbers for 'd'.
            if(abs(d) < numeric_limits<double>::min())
                d = numeric_limits<double>::min();
            c = b + an / c; // Apply the recurrence relation for C_n (numerator part of the fraction).
            // Safeguard against division by zero or extremely small numbers for 'c'.
            if(abs(c) < numeric_limits<double>::min())
                c = numeric_limits<double>::min();
            d = 1.0 / d;    // Invert 'd' for the next step.
            double del = d * c; // Calculate the current term (delta_n) of the continued fraction.
            h *= del;       // Multiply `h` by the current term to update the accumulated fraction value.
            // Check for convergence: If the current term `del` is very close to 1.0, the continued fraction has converged.
            if(abs(del - 1.0) < LOCAL_EPS)
                break; // Exit the loop if converged.
        }
        // Calculate the exponential and power part of the incomplete gamma function: `exp(-x + a*log(x) - gln)`.
        double lnPart = -x + a * log(x) - gln;
        // `Qval` represents Γ(a,x)/Γ(a) (the normalized upper incomplete gamma function)
        // as computed by the continued fraction method.
        double Qval = exp(lnPart) * h;
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
        throw invalid_argument("incomplete_gamma: 'a' must be positive.");
    }
    if (x < 0.0) {
        throw invalid_argument("incomplete_gamma: 'x' must be non-negative.");
    }
    // Calculation: Γ(a,x) = Q(a,x) * Γ(a).
    // `std::tgamma` computes the complete gamma function Γ(a).
    return incomplete_gamma_q(a, x) * tgamma(a); 
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
        throw invalid_argument("calculate_HN: N must be greater than 1 for log(N).");
    }
    if (H1 <= 0.0 || H2 <= 0.0) {
        throw invalid_argument("calculate_HN: H1 and H2 must be positive.");
    }
    if (k1 <= 0.0 || k2 <= 0.0) {
        throw invalid_argument("calculate_HN: k1 and k2 must be positive.");
    }

    // Calculate a candidate HN assuming it falls within the *first* part of the distribution (H <= Htr).
    // This uses Equation A.5 from Groenendijk (1998): `H_N,1 = H1 * (ln(N))^(1/k1)`.
    double HN_candidate1 = H1 * pow(log(N), 1.0 / k1);

    // Decision point: Check if the `HN_candidate1` is indeed less than the transitional wave height (Htr).
    // A small `EPSILON` is subtracted from `Htr` to account for floating-point inaccuracies when comparing.
    if (HN_candidate1 < Htr - EPSILON) {
        // If true, it means the wave height for the given exceedance probability is governed by the first Weibull part.
        return HN_candidate1;
    } else {
        // If `HN_candidate1` is not less than `Htr`, it implies that the wave height for the given
        // exceedance probability is determined by the *second* part of the distribution (H > Htr).
        // This uses Equation A.7 from Groenendijk (1998): `H_N = H2 * (ln(N))^(1/k2)`.
        return H2 * pow(log(N), 1.0 / k2);
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
        throw invalid_argument("calculate_H1N: H1 and H2 must be positive.");
    }
    if (k1 <= 0.0 || k2 <= 0.0) {
        throw invalid_argument("calculate_H1N: k1 and k2 must be positive.");
    }
    if (N_val <= 1.0) {
        throw invalid_argument("calculate_H1N: N_val must be greater than 1.");
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
        double term1_x_ln_Nval = log(N_val);
        // The 'x' parameter for the second incomplete gamma term (related to Htr) is `(Htr/H1)^k1`.
        double term1_x_HtrH1 = pow(Htr / H1, k1);

        // Calculate the two incomplete gamma terms for the first part of the integral.
        // This represents `Γ[a, ln(N)] - Γ[a, (Htr/H1)^k1]`.
        double gamma_term1_part1 = incomplete_gamma(term1_a, term1_x_ln_Nval);
        double gamma_term1_part2 = incomplete_gamma(term1_a, term1_x_HtrH1);
        double gamma_term1 = gamma_term1_part1 - gamma_term1_part2;

        // Term 2: Contribution from the second Weibull distribution (F2(H)).
        // The 'a' parameter for the incomplete gamma function is (1/k2 + 1).
        double term2_a = 1.0 / k2 + 1.0;
        // The 'x' parameter for this term is `(Htr/H2)^k2`.
        double term2_x_HtrH2 = pow(Htr / H2, k2);

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
        double term_x = log(N_val);     // The 'x' parameter is `ln(N_val)`.
        // Calculate `N_val * H2 * Γ[a, ln(N_val)]`.
        return N_val * H2 * incomplete_gamma(term_a, term_x);
    }
}

// --- Numerical Solver for H1/Hrms and H2/Hrms ---
// The following functions implement the numerical methods required to solve the non-linear
// equation that ensures the normalized Hrms of the Composite Weibull distribution equals one.

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
    double H2_Hrms = Htr_Hrms * pow(H1_Hrms / Htr_Hrms, k1 / k2);
    
    // Calculate the arguments for the incomplete gamma functions. These arguments are derived
    // from the integration limits (0 to Htr_Hrms for the first part, Htr_Hrms to infinity for the second)
    // and the specific form of the Weibull probability density function when calculating moments.
    double arg1 = pow(Htr_Hrms / H1_Hrms, k1); // Argument for the incomplete gamma function related to the first Weibull part.
    double arg2 = pow(Htr_Hrms / H2_Hrms, k2); // Argument for the incomplete gamma function related to the second Weibull part.
    
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
    return sqrt(sum_terms) - 1.0;
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
        if (abs(f_prime_x) < numeric_limits<double>::min()) { 
            return x_old; // Return the current best guess as the solver cannot proceed reliably.
        }

        // Newton-Raphson iteration formula: x_{n+1} = x_n - f(x_n) / f'(x_n).
        x_new = x_old - fx / f_prime_x;

        // Check for convergence: If the absolute difference between the new and old guesses
        // is less than the specified tolerance, the solution has converged.
        if (abs(x_new - x_old) < tol) {
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
    double logHtr = log(Htr_Hrms); // Calculate the natural logarithm of Htr_Hrms.
                                   // Logarithmic transformations are common in empirical fits.
    // This is the empirical regression formula. The coefficients (-1.705, -2.329, -0.313)
    // are specific to this model and were likely determined through curve fitting to data.
    return 1.0 + exp(-1.705 - 2.329 * logHtr - 0.313 * logHtr * logHtr);
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
    const double tol = 1e-14; // The very high desired tolerance for the root, indicating extreme precision.
    const int maxit = 200;    // The maximum number of iterations allowed for both Newton-Raphson and bisection.

    // Attempt to solve for H1_Hrms using the Newton-Raphson method.
    H1_Hrms = newtonRaphsonSolver(Htr_Hrms, initial_guess, tol, maxit);

    // After attempting Newton-Raphson, check if the residual at the found H1_Hrms is sufficiently close to zero.
    // A slightly relaxed tolerance (tol * 100) is used for this initial check to avoid overly strict failure.
    if (abs(residual(H1_Hrms, Htr_Hrms)) > tol * 100) { 
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
                if (abs(fc) < tol || abs(b - a) < tol) {
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
        cerr << "Solver: Failed to converge for Htr_Hrms = " << Htr_Hrms << " after Newton-Raphson and fallback.\n";
        return false; // Indicate that the solver failed to converge.
    }
    return true; // Indicate that Newton-Raphson successfully converged.
}

/**
 * @brief Main function of the program.
 *
 * This is the entry point of the C++ program. It orchestrates the entire computation
 * and output process. The function iterates through a predefined range of normalized
 * transitional wave heights (Htr_Hrms), for each value performing the necessary
 * numerical solves and calculations to determine the corresponding normalized
 * wave parameters of the Composite Weibull distribution. Finally, it writes all
 * computed results to a formatted text file.
 *
 * @return 0 if the program executes successfully, 1 if there is an error (e.g., file opening failure).
 */
int main()
{
    // Open an output file named "carvalho2025_table.txt" for writing.
    // `std::ofstream` is used to create an output file stream.
    ofstream fout("carvalho2025_table.txt");
    // Check if the file was successfully opened. If `fout` is in a bad state (e.g., due to permissions or disk full),
    // it evaluates to `false`.
    if (!fout) {
        cerr << "Error opening output file: carvalho2025_table.txt" << endl; // Print an error message to console.
        return 1; // Return a non-zero exit code to indicate an error.
    }

    // Define the range and step size for the normalized transitional wave height (Htr/Hrms).
    // These constants control the granularity and extent of the generated table.
    const double Htr_Hrms_start = 0.01; // The starting value for Htr/Hrms.
    const double Htr_Hrms_end = 3.50;   // The ending value for Htr/Hrms.
    const double Htr_Hrms_step = 0.01;  // The increment step size for Htr/Hrms.
    
    // Define an array of integer N values for which the H(1/N)/Hrms quantiles will be calculated.
    // These N values correspond to specific exceedance probabilities or mean-of-highest-N-part values.
    // For example, N=3 corresponds to H1/3 (significant wave height).
    const int N_values[] = {3, 10, 50, 100, 250, 1000}; 

    // Write the header row to the output file.
    // `std::setw(X)` is a manipulator that sets the field width for the *next* output item to X characters.
    // This ensures that the columns in the output file are neatly aligned.
    fout << setw(9) << "Htr/Hrms"
         << setw(9) << "H1/Hrms"
         << setw(9) << "H2/Hrms"
         << setw(11) << "H1/3/Hrms"    // Normalized significant wave height.
         << setw(11) << "H1/10/Hrms"   // Normalized mean of the highest 1/10th waves.
         << setw(11) << "H1/50/Hrms"   // Normalized mean of the highest 1/50th waves.
         << setw(12) << "H1/100/Hrms"  // Normalized mean of the highest 1/100th waves (often H1%).
         << setw(12) << "H1/250/Hrms"  // Normalized mean of the highest 1/250th waves.
         << setw(13) << "H1/1000/Hrms" // Normalized mean of the highest 1/1000th waves.
         << endl; // Insert a newline character to move to the next line after the header.

    // Main loop: Iterate through the defined range of Htr_Hrms values.
    // The `+ EPSILON` in the loop condition helps to ensure that the `Htr_Hrms_end` value
    // is included in the iteration, compensating for potential floating-point arithmetic inaccuracies.
    for (double Htr_Hrms = Htr_Hrms_start; Htr_Hrms <= Htr_Hrms_end + EPSILON; Htr_Hrms += Htr_Hrms_step) {
        try { // Begin a try-catch block to gracefully handle potential runtime errors during calculations.
            double H1_normalized; // Declare a variable to store the computed normalized H1 (H1/Hrms).
            double H2_normalized; // Declare a variable to store the computed normalized H2 (H2/Hrms).
            
            // Call the `highAccuracySolver` function to find the value of H1_normalized.
            // This function attempts to solve the non-linear equation for H1/Hrms.
            if (!highAccuracySolver(Htr_Hrms, H1_normalized)) {
                // If the solver returns `false`, it indicates a failure to converge, so throw a runtime_error.
                throw runtime_error("Solver failed to find H1_Hrms.");
            }

            // Calculate H2_normalized based on the computed H1_normalized and the continuity condition.
            // The formula used is derived from `(Htr/H1)^k1 = (Htr/H2)^k2`, rearranged for H2_Hrms.
            H2_normalized = Htr_Hrms * pow(H1_normalized / Htr_Hrms, k1 / k2);

            // Create a dynamic array (vector) to store the calculated H1/N values for the current Htr_Hrms.
            vector<double> h1n_values; 
            // Iterate through each predefined N value.
            for (int N_val : N_values) { 
                // Calculate the H1/N value using the `calculate_H1N` function and add it to the vector.
                // `static_cast<double>(N_val)` ensures that N_val is treated as a double in the function call.
                h1n_values.push_back(calculate_H1N(static_cast<double>(N_val), H1_normalized, H2_normalized, k1, k2, Htr_Hrms));
            }

            // Set output formatting for the current row in the file.
            // `std::fixed` ensures floating-point numbers are printed in fixed-point notation (not scientific).
            // `std::setprecision(5)` sets the number of digits after the decimal point to 5.
            fout << fixed << setprecision(5); 
            // Write the primary normalized parameters (Htr/Hrms, H1/Hrms, H2/Hrms) to the file,
            // formatted with specified widths.
            fout << setw(9) << Htr_Hrms
                 << setw(9) << H1_normalized
                 << setw(9) << H2_normalized;
            
            // Output each calculated H1/N value from the vector to the file,
            // also formatted with specified widths.
            fout << setw(11) << h1n_values[0]  // H1/3/Hrms
                 << setw(11) << h1n_values[1]  // H1/10/Hrms
                 << setw(11) << h1n_values[2]  // H1/50/Hrms
                 << setw(12) << h1n_values[3]  // H1/100/Hrms
                 << setw(12) << h1n_values[4]  // H1/250/Hrms
                 << setw(13) << h1n_values[5]  // H1/1000/Hrms
                 << endl; // Move to the next line in the output file.

        } catch (const invalid_argument& e) { // Catch `std::invalid_argument` exceptions.
            // These typically occur due to invalid input parameters passed to functions (e.g., log of non-positive number).
            cerr << "Input error for Htr_Hrms = " << Htr_Hrms << ": " << e.what() << endl;
            // If an error occurs, write "ERROR" for all calculated values in the output file for this row.
            fout << setw(9) << Htr_Hrms
                 << setw(9) << "ERROR"
                 << setw(9) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(12) << "ERROR" 
                 << setw(13) << "ERROR"
                 << endl;
        } catch (const runtime_error& e) { // Catch `std::runtime_error` exceptions.
            // These typically indicate issues during the execution, such as solver failure to converge.
            cerr << "Calculation error for Htr_Hrms = " << Htr_Hrms << ": " << e.what() << endl;
            // If an error occurs, write "ERROR" for all calculated values in the output file for this row.
            fout << setw(9) << Htr_Hrms
                 << setw(9) << "ERROR"
                 << setw(9) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(12) << "ERROR" 
                 << setw(13) << "ERROR"
                 << endl;
        } catch (const exception& e) { // Catch any other standard C++ exceptions.
            // This is a general catch-all for unexpected errors.
            cerr << "An unexpected error occurred for Htr_Hrms = " << Htr_Hrms << ": " << e.what() << endl;
            // If an error occurs, write "ERROR" for all calculated values in the output file for this row.
            fout << setw(9) << Htr_Hrms
                 << setw(9) << "ERROR"
                 << setw(9) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(11) << "ERROR"
                 << setw(12) << "ERROR"
                 << setw(12) << "ERROR" 
                 << setw(13) << "ERROR"
                 << endl;
        }
    }

    fout.close(); // Close the output file stream, ensuring all buffered data is written to the file.
    cout << "Calculation complete. Results saved to carvalho2025_table.txt" << endl; // Inform the user via console.

    return 0; // Return 0 to indicate that the program executed successfully without critical errors.
}
