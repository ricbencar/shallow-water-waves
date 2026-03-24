/***********************************************************************
 * Program: rayleigh_battjes.cpp
 * Author: Gemini (Updated)
 *
 * Description:
 * This program calculates and compares an extensive set of theoretical
 * wave height ratios for two different deep-water wave models, for N
 * ranging from 1 to 10,000.
 *
 * The header of the output file contains an extremely detailed set of
 * constant diagnostic ratios for the Rayleigh distribution and its
 * underlying JONSWAP and Pierson-Moskowitz (PM) spectral shapes.
 *
 * Compilation (using g++):
 * g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic \
 * -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
 * -o rayleigh_battjes rayleigh_battjes.cpp
 *
 * Usage:
 * ./rayleigh_battjes
 ***********************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

// Use double for high-precision floating-point calculations.
using Real = double;

// Define mathematical constants with high precision.
const Real PI = acos(-1.0);
// Euler-Mascheroni constant (gamma), value from C++20 <numbers> header
const Real EULER_MASCHERONI = 0.57721566490153286;


// --- FORWARD DECLARATIONS ---

Real get_H1N_over_Hrms(int N);
Real get_HN_over_Hrms(int N);
Real get_Hmax_over_Hrms(int N);

// Structure to hold pre-computed constant ratios to avoid recalculation.
struct PrecomputedRatios {
    Real Hm0_over_Hrms_Rayleigh; // For pure Rayleigh: 4/sqrt(8) = sqrt(2)
    Real Hm0_over_Hrms_BG;       // For B&G convention: 4/2.69
};


/**
 * @brief Helper function to perform all calculations and write a formatted row.
 * @param N The current N value.
 * @param outFile The output file stream.
 * @param ratios A struct containing all pre-computed constant ratios.
 */
void calculateAndWrite(int N, std::ofstream& outFile, const PrecomputedRatios& ratios) {
    // --- N-dependent calculations ---
    // Base ratio for the current N (valid for any Rayleigh-shaped distribution).
    Real H1_N_over_Hrms = get_H1N_over_Hrms(N);
    // Theoretical maximum wave height for a sea state of N waves.
    Real Hmax_over_Hrms = get_Hmax_over_Hrms(N);

    // --- Final Ratio Calculations ---
    Real ratio_Rayleigh_vs_Hm0 = H1_N_over_Hrms / ratios.Hm0_over_Hrms_Rayleigh;
    Real ratio_BG_vs_Hm0 = H1_N_over_Hrms / ratios.Hm0_over_Hrms_BG;
    Real Hmax_over_Hm0 = Hmax_over_Hrms / ratios.Hm0_over_Hrms_Rayleigh;

    // --- Write Row to File ---
    outFile << std::left
            << std::setw(8) << N
            << std::setw(25) << H1_N_over_Hrms
            << std::setw(25) << ratio_Rayleigh_vs_Hm0
            << std::setw(25) << ratio_BG_vs_Hm0
            << std::setw(25) << Hmax_over_Hm0 << "\n";
}

/**
 * @brief Helper function to evaluate the polynomial approximations from ITTC Table A.4
 * for JONSWAP spectral properties.
 * @param gamma The peak enhancement factor.
 * @param c0 Constant term.
 * @param c1 Coefficient for gamma.
 * @param c2 Coefficient for gamma^2.
 * @param c3 Coefficient for gamma^3.
 * @return The result of the polynomial evaluation.
 */
Real eval_jonswap_poly(Real gamma, Real c0, Real c1, Real c2, Real c3) {
    return c0 + c1 * gamma + c2 * pow(gamma, 2) + c3 * pow(gamma, 3);
}


int main() {
    // Attempt to open the output file for writing.
    std::ofstream outFile("output.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output.txt for writing." << std::endl;
        return 1;
    }

    // Set output stream formatting for the entire file.
    outFile << std::fixed << std::setprecision(15);

    // --- Pre-compute all constant ratios ---
    PrecomputedRatios ratios;
    ratios.Hm0_over_Hrms_Rayleigh = sqrt(2.0); // Theoretical ratio: 4.0 / sqrt(8.0)
    ratios.Hm0_over_Hrms_BG = 4.0 / 2.69;      // B&G empirical ratio

    // --- Calculate Diagnostic Ratios for Header (Based on Pure Rayleigh) ---
    Real H_mean_over_Hrms = get_H1N_over_Hrms(1);
    Real H1_3_over_Hrms = get_H1N_over_Hrms(3);
    Real H1_10_over_Hrms = get_H1N_over_Hrms(10);
    Real H1_100_over_Hrms = get_H1N_over_Hrms(100);
    Real H1_1000_over_Hrms = get_H1N_over_Hrms(1000);
    Real Hmax10000_over_Hrms = get_Hmax_over_Hrms(10000);
    Real Hmax3000_over_Hrms = get_Hmax_over_Hrms(3000); // For rule-of-thumb check
    Real H_mode_over_Hrms = 1.0 / sqrt(2.0);
    Real H_median_over_Hrms = get_HN_over_Hrms(2);
    Real H_10_percent_over_Hrms = get_HN_over_Hrms(10);
    
    // Fundamental statistical properties
    Real coeff_of_variation = sqrt(4.0 / PI - 1.0);
    Real skewness = (2.0 * sqrt(PI) * (PI - 3.0)) / pow(4.0 - PI, 1.5);
    Real kurtosis = (-6.0 * pow(PI, 2) + 24.0 * PI - 16.0) / pow(4.0 - PI, 2);
    Real IQR_over_Hrms = sqrt(log(4.0)) - sqrt(log(4.0/3.0));
    Real IQR_div_Hmean = IQR_over_Hrms / H_mean_over_Hrms;

    // Central tendency ratios
    Real Hrms_div_Hmean = 1.0 / H_mean_over_Hrms;
    Real Hs_div_Hmean = H1_3_over_Hrms / H_mean_over_Hrms;
    Real Hmode_div_Hmean = H_mode_over_Hrms / H_mean_over_Hrms;
    Real Hmedian_div_Hmean = H_median_over_Hrms / H_mean_over_Hrms;
    Real Hrms_div_Hmode = 1.0 / H_mode_over_Hrms;

    // Significant & characteristic wave ratios
    Real H1_3_over_Hm0 = H1_3_over_Hrms / ratios.Hm0_over_Hrms_Rayleigh;
    Real H1_10_over_Hs = H1_10_over_Hrms / H1_3_over_Hrms;
    Real H1_100_over_Hs = H1_100_over_Hrms / H1_3_over_Hrms;
    Real H_10_percent_over_Hs = H_10_percent_over_Hrms / H1_3_over_Hrms;

    // --- Underlying spectral properties ---
    // JONSWAP Calculations (gamma = 3.3)
    Real gamma_j = 3.3;
    Real m0_j_prop = eval_jonswap_poly(gamma_j, 0.1475, 0.05617, -0.003077, 0.0001618);
    Real m1_j_prop = eval_jonswap_poly(gamma_j, 0.2059, 0.05705, -0.003154, 0.0001661);
    Real m2_j_prop = eval_jonswap_poly(gamma_j, 0.3420, 0.05827, -0.003269, 0.0001723);

    Real fE_div_fp_j = 1.0 / eval_jonswap_poly(gamma_j, 0.8255, 0.03852, -0.005537, 0.0003154);
    Real fmean_div_fp_j = 1.0 / eval_jonswap_poly(gamma_j, 0.7303, 0.04936, -0.006556, 0.0003610);
    Real fz_div_fp_j = 1.0 / eval_jonswap_poly(gamma_j, 0.6673, 0.05037, -0.006230, 0.0003341);
    
    Real Tp_div_TE_j = fE_div_fp_j;
    Real Tp_div_Tmean_j = fmean_div_fp_j;
    Real Tp_div_Tz_j = fz_div_fp_j;
    Real moment_ratio_j = m1_j_prop / sqrt(m0_j_prop * m2_j_prop);
    Real nu_j = sqrt((m0_j_prop * m2_j_prop) / pow(m1_j_prop, 2) - 1.0);

    // PM Calculations (gamma = 1.0), direct from ITTC PDF Table A.2
    Real m0_pm = 1.0; // Proportionality constant
    Real m1_pm = tgamma(0.75);
    Real m2_pm = sqrt(PI);
    Real moment_ratio_pm = m1_pm / sqrt(m0_pm * m2_pm);
    Real nu_pm = sqrt((m0_pm * m2_pm) / pow(m1_pm, 2) - 1.0);
    Real Tp_div_Tz_pm = 1.331 / 0.946;
    Real Tp_div_Tmean_pm = 1.225 / 0.946;
    Real Tp_div_TE_pm = 1.103 / 0.946;


    // Extreme wave ratios
    Real H1000_div_H100_ratio = H1_1000_over_Hrms / H1_100_over_Hrms;
    Real E_max_div_MP_max_10000 = 1.0 + (0.5 * EULER_MASCHERONI) / log(10000.0);
    Real Hmax10000_over_Hm0_ratio = Hmax10000_over_Hrms / ratios.Hm0_over_Hrms_Rayleigh;
    Real Hmax3000_over_Hs_ratio = Hmax3000_over_Hrms / H1_3_over_Hrms;
    Real Hmax10000_over_H1_100_ratio = Hmax10000_over_Hrms / H1_100_over_Hrms;
    Real prob_rogue_2_5 = exp(-12.5);

    // --- Write Introductory Note ---
    outFile << "# Wave Height Ratio Comparison: Theoretical Rayleigh vs. Battjes & Groenendijk (2000) Convention\n\n"
            << "# This file calculates and compares wave height statistics for two fundamental deep-water models.\n"
            << "# Both models are based on the Rayleigh probability distribution, but they differ in how they relate wave statistics to the total energy of the sea state.\n\n"
            << "# --- Model Descriptions ---\n\n"
            << "# 1. Theoretical Rayleigh Model:\n"
            << "#    - Context: This is the pure, classical model for wave heights, first derived for sea waves by Longuet-Higgins (1952).\n"
            << "#    - Assumption: It assumes the sea surface is a linear Gaussian process with a narrow energy spectrum (like a clean, long-distance swell).\n"
            << "#    - Key Relationship: It uses the theoretical link: Hrms = sqrt(8 * m0) approx 2.828 * sqrt(m0).\n\n"
            << "# 2. Battjes & Groenendijk (B&G) Deep-Water Convention:\n"
            << "#    - Context: This is an empirical adjustment to the Rayleigh model's parameters to better fit real, wind-driven seas with broad energy spectra.\n"
            << "#      It represents the DEEP-WATER LIMIT of the full B&G shallow-water model.\n"
            << "#    - Key Relationship: It uses an empirical link based on field data (Goda, 1979): Hrms approx. 2.69 * sqrt(m0).\n\n"
            << "# --- Constant Ratios for a Pure Rayleigh Distribution ---\n"
            << "# These diagnostic ratios are constant and describe the fundamental shape and properties of the distribution.\n\n"
            << "# --- Fundamental Statistical Properties ---\n"
            << "# Coeff. of Variat. (StDev/Mean): " << std::left << std::setw(20) << coeff_of_variation << " (Dimensionless measure of spread)\n"
            << "# Skewness                      : " << std::left << std::setw(20) << skewness << " (Measure of asymmetry; positive indicates a tail to the right)\n"
            << "# Kurtosis (non-excess)         : " << std::left << std::setw(20) << kurtosis << " (Measure of 'tailedness'; Normal dist = 3)\n"
            << "# IQR / Hmean                   : " << std::left << std::setw(20) << IQR_div_Hmean << " (Normalized Interquartile Range)\n\n"
            << "# --- Central Tendency & Common Wave Ratios ---\n"
            << "# Hrms / Hmean                  : " << std::left << std::setw(20) << Hrms_div_Hmean << " (Ratio of root-mean-square to the mean wave height)\n"
            << "# Hs / Hmean                    : " << std::left << std::setw(20) << Hs_div_Hmean << " (Ratio of significant to the mean wave height)\n"
            << "# Hrms / Hmode                  : " << std::left << std::setw(20) << Hrms_div_Hmode << " (Ratio of RMS to most probable height; should be sqrt(2))\n"
            << "# Hmode / Hmean                 : " << std::left << std::setw(20) << Hmode_div_Hmean << " (Ratio of the most probable wave height to the mean)\n"
            << "# Hmedian / Hmean               : " << std::left << std::setw(20) << Hmedian_div_Hmean << " (Ratio of the median wave height to the mean)\n\n"
            << "# --- Significant & Characteristic Wave Ratios ---\n"
            << "# Hs / Hrms                     : " << std::left << std::setw(20) << H1_3_over_Hrms << " (Also H(1/3)/Hrms; should be approx. sqrt(2) = 1.414)\n"
            << "# H(1/3) / Hm0                  : " << std::left << std::setw(20) << H1_3_over_Hm0 << " (Statistical H(1/3)=4.004*sqrt(m0) vs spectral Hm0=4*sqrt(m0))\n"
            << "# H(1/10) / Hs                  : " << std::left << std::setw(20) << H1_10_over_Hs << " (Mean of highest 10% vs significant; approx. 1.27)\n"
            << "# H(1/100) / Hs                 : " << std::left << std::setw(20) << H1_100_over_Hs << " (Mean of highest 1% vs significant; approx. 1.67)\n"
            << "# H(10%) / Hs                   : " << std::left << std::setw(20) << H_10_percent_over_Hs << " (Wave height exceeded by 10% of waves vs significant height)\n\n"
            << "# --- Underlying Spectral Properties & Shape (for JONSWAP Spectrum, gamma=3.3) ---\n"
            << "# Theory                        : The JONSWAP spectrum describes developing seas with a higher peak than a fully-developed sea.\n"
            << "# Spectral Width (nu)           : " << std::left << std::setw(20) << nu_j << " (Bandwidth parameter based on moments m0, m1, m2)\n"
            << "# Moment Ratio m1/sqrt(m0*m2)   : " << std::left << std::setw(20) << moment_ratio_j << " (Dimensionless spectral moment ratio describing shape)\n"
            << "# Tp / Tz (Peak/Zero-Cross)     : " << std::left << std::setw(20) << Tp_div_Tz_j << " (Ratio of peak period to zero-crossing period)\n"
            << "# Tp / Tmean (Peak/Mean)        : " << std::left << std::setw(20) << Tp_div_Tmean_j << " (Ratio of peak period to mean period)\n"
            << "# Tp / TE (Peak/Energy)         : " << std::left << std::setw(20) << Tp_div_TE_j << " (Ratio of peak period to energy period)\n\n"
            << "# --- Underlying Spectral Properties & Shape (for Pierson-Moskowitz Spectrum) ---\n"
            << "# Theory                        : The PM spectrum (JONSWAP with gamma=1.0) describes a 'fully developed' sea in equilibrium.\n"
            << "# High-Frequency Tail           : The PM spectrum assumes an f^-5 decay in the high-frequency range, based on Phillips' equilibrium range theory.\n"
            << "# Spectral Width (nu)           : " << std::left << std::setw(20) << nu_pm << " (Bandwidth parameter based on moments m0, m1, m2)\n"
            << "# Moment Ratio m1/sqrt(m0*m2)   : " << std::left << std::setw(20) << moment_ratio_pm << " (Dimensionless spectral moment ratio describing shape)\n"
            << "# Tp / Tz (Peak/Zero-Cross)     : " << std::left << std::setw(20) << Tp_div_Tz_pm << " (Ratio of peak period to zero-crossing period)\n"
            << "# Tp / Tmean (Peak/Mean)        : " << std::left << std::setw(20) << Tp_div_Tmean_pm << " (Ratio of peak period to mean period)\n"
            << "# Tp / TE (Peak/Energy)         : " << std::left << std::setw(20) << Tp_div_TE_pm << " (Ratio of peak period to energy period)\n\n"
            << "# --- Extreme & Rogue Wave Statistics ---\n"
            << "# H(1/1000) / H(1/100)          : " << std::left << std::setw(20) << H1000_div_H100_ratio << " (Growth ratio of extreme wave averages)\n"
            << "# Hmax(N=10000) / H(1/100)      : " << std::left << std::setw(20) << Hmax10000_over_H1_100_ratio << " (Expected max wave vs. the mean of the top 1%)\n"
            << "# E[Hmax] / MP[Hmax] (N=10k)    : " << std::left << std::setw(20) << E_max_div_MP_max_10000 << " (Ratio of Expected Maximum to Most Probable Maximum)\n"
            << "# Hmax(N=10000) / Hm0           : " << std::left << std::setw(20) << Hmax10000_over_Hm0_ratio << " (Expected max wave in 10000 waves vs. spectral significant wave)\n"
            << "# Hmax(N=3000) / Hs             : " << std::left << std::setw(20) << Hmax3000_over_Hs_ratio << " (Rule-of-Thumb check; should be approx. 2.0 for a typical storm)\n"
            << "# Prob(H > 2.5*Hs)              : " << std::scientific << std::setprecision(4) << prob_rogue_2_5 << std::fixed << std::setprecision(15) << " (Prob. of a 'rogue wave'; implies 1 in " << int(1.0/prob_rogue_2_5) << " waves)\n\n"
            << "# --- Column Explanations ---\n"
            << "# N               : Number of waves in the sea state.\n"
            << "# H(1/N)/Hrms     : The mean of the highest 1/N waves, normalized by Hrms. This is the base ratio derived from the Rayleigh shape.\n"
            << "# Rayleigh/Hm0    : Theoretical model ratio, normalized by the spectral significant wave height (Hm0 = 4*sqrt(m0) using Hrms=sqrt(8*m0)).\n"
            << "# B&G/Hm0         : B&G convention ratio, normalized by spectral significant wave height (Hm0 = 4*sqrt(m0) using Hrms=2.69*sqrt(m0)).\n"
            << "# Hmax/Hm0        : The ratio of the expected maximum wave height in N waves to the spectral significant wave height.\n"
            << "-----------------------------------------------------------------------------------------------------------------\n";

    // --- Write Header ---
    outFile << std::left
            << std::setw(8) << "N"
            << std::setw(25) << "H(1/N)/Hrms"
            << std::setw(25) << "Rayleigh/Hm0"
            << std::setw(25) << "B&G/Hm0"
            << std::setw(25) << "Hmax/Hm0" << "\n";


    // --- Main Calculation Loop with Custom Steps ---
    for (int N = 1; N < 10; ++N) {
        calculateAndWrite(N, outFile, ratios);
    }
    for (int N = 10; N <= 100; N += 10) {
        calculateAndWrite(N, outFile, ratios);
    }
    for (int N = 150; N < 1000; N += 50) {
        calculateAndWrite(N, outFile, ratios);
    }
    for (int N = 1000; N <= 10000; N += 1000) {
        calculateAndWrite(N, outFile, ratios);
    }

    // Close the file stream.
    outFile.close();

    // Print a success message to the console.
    std::cout << "Successfully calculated and compared wave height ratios for N=1 to 10000." << std::endl;
    std::cout << "Results have been written to output.txt" << std::endl;

    return 0; // Indicate successful execution.
}


// --- FUNCTION DEFINITIONS ---

/**
 * @brief Calculates the ratio H(1/N)/Hrms for a pure Rayleigh distribution.
 * H(1/N) is the mean of the highest 1/N waves.
 * The formula is: H(1/N)/Hrms = sqrt(ln(N)) + N * sqrt(PI)/2 * erfc(sqrt(ln(N)))
 */
Real get_H1N_over_Hrms(int N) {
    if (N <= 0) return 0.0;
    if (N == 1) return sqrt(PI) / 2.0;
    Real n_real = static_cast<Real>(N);
    Real log_N = log(n_real);
    Real sqrt_log_N = sqrt(log_N);
    Real term1 = sqrt_log_N;
    Real term2 = n_real * (sqrt(PI) / 2.0) * erfc(sqrt_log_N);
    return term1 + term2;
}

/**
 * @brief Calculates the ratio H_N/Hrms for a pure Rayleigh distribution.
 * H_N is the wave height with an exceedance probability of 1/N.
 * The formula is: H_N/Hrms = sqrt(ln(N))
 */
Real get_HN_over_Hrms(int N) {
    if (N < 1) return 0.0; // Probability must be <= 1.
    if (N == 1) return 0.0; // H exceeded by all waves is 0.
    return sqrt(log(static_cast<Real>(N)));
}

/**
 * @brief Calculates the approximate ratio H_max/Hrms for a given N.
 * H_max is the expected maximum wave height in a sequence of N waves.
 * This formula represents the MOST PROBABLE maximum.
 * The formula is: H_max/Hrms = sqrt(2 * ln(N))
 */
Real get_Hmax_over_Hrms(int N) {
    if (N <= 1) return get_H1N_over_Hrms(1); // For N=1, H_max is just the mean.
    return sqrt(2.0 * log(static_cast<Real>(N)));
}