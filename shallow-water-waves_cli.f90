!***********************************************************************
! Program: shallow-water-waves_cli.f90
!
! Detailed Description:
! This program computes local shallow-foreshore wave-height distribution
! parameters using a model based on the Composed Weibull distribution.
!
! The command-line application performs the following:
! 1. If three command-line arguments are provided, they are used as
! Hm0 (local significant spectral spectral wave height), d (local water depth),
! and slopeM (beach slope 1:m). Otherwise, the program prompts the user
! for these values.
! 2. Computes the following intermediate values:
! - Free-surface variance: m0 = (Hm0 / 4)²
!
! - Mean square wave height:
! Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0)
!
! - A dimensional transitional wave height:
! Htr = (0.35 + 5.8*(1/m)) * d.
!
! - The dimensionless transitional parameter:
! H̃_tr = Htr / Hrms.
!
! 3. Calculates the dimensionless wave-height ratios (Hᵢ/Hrms)
! by solving a system of non-linear equations derived from the Composite Weibull
! distribution, ensuring the normalized Hrms of the distribution equals one.
! This involves using a Newton-Raphson matrix method for simultaneous root-finding
! and functions for unnormalized incomplete gamma calculations. These ratios are then
! converted to dimensional quantities (in meters) by multiplying with Hrms.
!
! 4. A detailed report is then generated (and written to "report.txt") with
! the input parameters, intermediate values, calculated ratios and computed
! dimensional wave heights, as well as diagnostic ratios.
!
! Compilation Instructions (example using gfortran):
!
! gfortran -O3 -march=native -std=f2008 -Wall -Wextra -pedantic \
! -fno-underscoring -o shallow-water-waves_cli.exe shallow-water-waves_cli.f90
!
! To run with command-line arguments (e.g., Hm0=2.5, d=10, slopeM=20):
! shallow-water-waves_cli 2.5 10 20
!***********************************************************************

PROGRAM ShallowWaterWaves
    IMPLICIT NONE

! Define a portable double precision kind
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)

! Global parameters for the Composite Weibull distribution.
    REAL(KIND=dp), PARAMETER :: k1 = 2.0_dp! Exponent for the first part of the Composite Weibull distribution (Rayleigh-shaped)
    REAL(KIND=dp), PARAMETER :: k2 = 3.6_dp! Exponent for the second part of the Composite Weibull distribution

! Precision for Numerical Solver Convergence:
    REAL(KIND=dp), PARAMETER :: EPSILON = 1.0E-12_dp! A small value (10^-12) indicating the maximum allowable error or difference
                                        ! within numerical computations.
    REAL(KIND=dp), PARAMETER :: JACOBIAN_DX = 1.0E-8_dp! Small step size for finite difference approximation of Jacobian derivatives.

    INTEGER, PARAMETER :: REPORT_UNIT = 10
    REAL(KIND=dp) :: Hm0, d, slopeM
    CHARACTER(LEN=256) :: arg_str
    INTEGER :: arg_count, iostat_val

! Main execution starts here

! Check command-line arguments
    arg_count = COMMAND_ARGUMENT_COUNT()

    IF (arg_count >= 3) THEN
    ! Read from command-line arguments
        CALL GET_COMMAND_ARGUMENT(1, arg_str)
        READ(arg_str, *, IOSTAT=iostat_val) Hm0
        IF (iostat_val /= 0) THEN
            WRITE(*,*) "Invalid argument: Hm0. Please enter a numeric value."
            ERROR STOP
        END IF

        CALL GET_COMMAND_ARGUMENT(2, arg_str)
        READ(arg_str, *, IOSTAT=iostat_val) d
        IF (iostat_val /= 0) THEN
            WRITE(*,*) "Invalid argument: water depth d. Please enter a numeric value."
            ERROR STOP
        END IF

        CALL GET_COMMAND_ARGUMENT(3, arg_str)
        READ(arg_str, *, IOSTAT=iostat_val) slopeM
        IF (iostat_val /= 0) THEN
            WRITE(*,*) "Invalid argument: beach slope m. Please enter a numeric value."
            ERROR STOP
        END IF
    ELSE
    ! Prompt user for input
        WRITE(*,*) "Enter Hm0 (m): "
        READ(*,*,IOSTAT=iostat_val) Hm0
        IF (iostat_val /= 0) THEN
            WRITE(*,*) "Input error for Hm0. Exiting."
            ERROR STOP
        END IF

        WRITE(*,*) "Enter water depth d (m): "
        READ(*,*,IOSTAT=iostat_val) d
        IF (iostat_val /= 0) THEN
            WRITE(*,*) "Input error for water depth d. Exiting."
            ERROR STOP
        END IF

        WRITE(*,*) "Enter beach slope (m): "
        READ(*,*,IOSTAT=iostat_val) slopeM
        IF (iostat_val /= 0) THEN
            WRITE(*,*) "Input error for beach slope (m). Exiting."
            ERROR STOP
        END IF
    END IF

! Validate input parameters
    IF (Hm0 <= 0.0_dp .OR. d <= 0.0_dp) THEN
        WRITE(*,*) "ERROR: Hm0 and d must be positive."
        CALL writeReportToFile("ERROR: Hm0 and d must be positive.")
        ERROR STOP
    END IF
    IF (slopeM <= 0.0_dp) THEN
        WRITE(*,*) "ERROR: Beach slope (m) must be positive."
        CALL writeReportToFile("ERROR: Beach slope (m) must be positive.")
        ERROR STOP
    END IF

    CALL buildReportAndDisplay(Hm0, d, slopeM)

CONTAINS

! All gamma functions applied, either complete or incomplete, are non-normalized.

    REAL(KIND=dp) FUNCTION incomplete_gamma_lower(a, x)
        IMPLICIT NONE
    ! @brief Computes the unnormalized lower incomplete gamma function γ(a, x) = ∫[0,x] t^(a-1)e^(-t) dt.
    !
    ! This function is a critical component for calculating moments of Weibull distributions,
    ! which are fundamental to the Composite Weibull model. It employs a hybrid approach
    ! for numerical stability and accuracy:
    ! - For small values of 'x' (specifically, x < a + 1.0), it uses a series expansion.
    ! - For larger values of 'x', it utilizes a continued fraction expansion.
    ! This adaptive strategy ensures robust and precise computation across different input ranges.
    ! This is a standard numerical implementation for the incomplete gamma function.
    !
    ! @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
    ! @param x The 'x' parameter (upper integration limit) of the incomplete gamma function. Must be non-negative.
    ! @return The computed value of γ(a, x). Returns `Huge(1.0_dp)` (representing an, error) if input parameters are invalid
    ! or if the series/continued fraction fails to converge within the maximum iterations.
        REAL(KIND=dp), INTENT(IN) :: a, x
        INTEGER, PARAMETER :: MAXIT = 500 ! Maximum number of iterations allowed for either the series expansion or the continued fraction.
        REAL(KIND=dp), PARAMETER :: LOCAL_EPS_GAMMA = 1.0E-16_dp! A very small tolerance value (10^-16) used for checking convergence
                                        ! within the gamma function calculations.
        REAL(KIND=dp) :: gln, ap, sum_val, del_val, b, c, d_val, h, an, del, exp_term
        INTEGER :: n_iter ! Explicitly declare n_iter

    ! Compute the natural logarithm of the complete gamma function, ln(Γ(a)).
    ! `LOG_GAMMA` is used for improved numerical stability.
        gln = LOG_GAMMA(a)

    ! Conditional logic to select the appropriate numerical method (series or continued fraction).
        IF (x < a + 1.0_dp) THEN! Use series expansion for x < a+1.
            ap = a
            sum_val = 1.0_dp / a
            del_val = sum_val
            DO n_iter = 1, MAXIT
                ap = ap + 1.0_dp
                del_val = del_val * x / ap
                sum_val = sum_val + del_val
                IF (ABS(del_val) < ABS(sum_val) * LOCAL_EPS_GAMMA) THEN
                    exp_term = EXP(-x + a * LOG(x) - gln)
                    incomplete_gamma_lower = sum_val * exp_term * GAMMA(a)
                    RETURN
                END IF
            END DO
            incomplete_gamma_lower = HUGE(1.0_dp)! Convergence failed
        ELSE! Use continued fraction for x >= a+1.
            b = x + 1.0_dp - a
            c = 1.0_dp / TINY(1.0_dp)! Equivalent to C++ numeric_limits<double>::min()
            d_val = 1.0_dp / b
            h = d_val
            DO n_iter = 1, MAXIT! Using n_iter for consistency
                an = -1.0_dp * REAL(n_iter, KIND=dp) * (REAL(n_iter, KIND=dp) - a)
                b = b + 2.0_dp
                d_val = an * d_val + b
                IF (ABS(d_val) < TINY(1.0_dp)) d_val = TINY(1.0_dp)
                c = b + an / c
                IF (ABS(c) < TINY(1.0_dp)) c = TINY(1.0_dp)
                d_val = 1.0_dp / d_val
                del = d_val * c
                h = h * del
                IF (ABS(del - 1.0_dp) < LOCAL_EPS_GAMMA) EXIT
            END DO
            exp_term = EXP(-x + a * LOG(x) - gln)
            incomplete_gamma_lower = (1.0_dp - exp_term * h) * GAMMA(a)
        END IF
    END FUNCTION incomplete_gamma_lower


    REAL(KIND=dp) FUNCTION incomplete_gamma_upper(a, x)
        IMPLICIT NONE
    ! @brief Computes the unnormalized upper incomplete gamma function Γ(a, x) = ∫[x,∞] t^(a-1)e^(-t) dt.
    !
    ! This function calculates the unnormalized form of the upper incomplete gamma function.
    ! It is derived from the complete gamma function Γ(a) and the unnormalized lower incomplete gamma function γ(a, x),
    ! using the identity Γ(a,x) = Γ(a) - γ(a,x).
    ! This is a standard numerical implementation for the incomplete gamma function.
    !
    ! @param a The 'a' parameter (shape parameter) of the incomplete gamma function. Must be positive.
    ! @param x The 'x' parameter (lower integration limit) of the incomplete gamma function. Must be non-negative.
    ! @return The computed value of the unnormalized upper incomplete gamma function Γ(a, x).
    ! Returns `Huge(1.0_dp)` if input parameters are invalid.
        REAL(KIND=dp), INTENT(IN) :: a, x

        IF (a <= 0.0_dp) THEN
            WRITE(*,*) "Error: incomplete_gamma_upper: 'a' must be positive."
            incomplete_gamma_upper = HUGE(1.0_dp)
            RETURN
        END IF
        IF (x < 0.0_dp) THEN
            WRITE(*,*) "Error: incomplete_gamma_upper: 'x' must be non-negative."
            incomplete_gamma_upper = HUGE(1.0_dp)
            RETURN
        END IF
        incomplete_gamma_upper = GAMMA(a) - incomplete_gamma_lower(a, x)

    END FUNCTION incomplete_gamma_upper

    REAL(KIND=dp) FUNCTION calculate_HN(N, H1, H2, k1_val, k2_val, Htr)
        IMPLICIT NONE
    ! @brief Calculates HN (wave height with 1/N exceedance probability) for the Composite Weibull distribution.
    !
    ! This function determines the specific wave height (H) such that the probability of a wave
    ! exceeding this height is 1/N. This is a key statistical measure. The calculation depends
    ! on whether the target wave height falls into the first or second part of the Composite Weibull
    ! distribution, which is separated by the transitional wave height (Htr).
    !
    ! @param N The N value (e.g., 3 for H1/3, 100 for H1%). Must be strictly greater than 1.0
    ! because `log(N)` is used, and `N=1` would result in `log(1)=0`.
    ! @param H1 The scale parameter of the first Weibull distribution. Must be positive.
    ! @param @param H2 The scale parameter of the second Weibull distribution. Must be positive.
    ! @param k1_val The exponent (shape parameter) of the first Weibull distribution. Must be positive.
    ! @param k2_val The exponent (shape parameter) of the second Weibull distribution. Must be positive.
    ! @param Htr The transitional wave height, which defines the boundary between the two Weibull parts.
    ! @return The calculated value of HN.
    ! Returns `Huge(1.0_dp)` if input parameters are invalid.
        REAL(KIND=dp), INTENT(IN) :: N, H1, H2, k1_val, k2_val, Htr
        REAL(KIND=dp) :: HN_candidate1, log_N_val, term_power1, term_power2

        log_N_val = LOG(N)

        term_power1 = (1.0_dp / k1_val)
        HN_candidate1 = H1 * (log_N_val)**term_power1

        IF (HN_candidate1 < Htr - EPSILON) THEN
            calculate_HN = HN_candidate1
        ELSE
            term_power2 = (1.0_dp / k2_val)
            calculate_HN = H2 * (log_N_val)**term_power2
        END IF
    END FUNCTION calculate_HN

    REAL(KIND=dp) FUNCTION calculate_H1N(N_val, H1, H2, k1_val, k2_val, Htr)
        IMPLICIT NONE
    ! @brief Calculates the mean of the highest 1/N-part of wave heights (H1/N) for the Composite Weibull distribution.
    !
    ! This function computes a characteristic wave height that represents the average height of the
    ! highest N-th fraction of waves in a given wave field. This is a commonly used metric in
    ! oceanography and coastal engineering (e.g., H1/3 for significant wave height).
    ! The calculation involves integrals of the probability density function and depends on
    ! whether the relevant wave heights fall into the first or second part of the Composite Weibull distribution.
    !
    ! @param N_val The N parameter for H1/N (e.g., 3 for H1/3, 10 for H1/10). Must be strictly greater than 1.
    ! @param H1 The scale parameter of the first Weibull distribution. Must be positive.
    ! @param H2 The scale parameter of the second Weibull distribution. Must be positive.
    ! @param k1_val The exponent (shape parameter) of the first Weibull distribution. Must be positive.
    ! @param k2_val The exponent (shape parameter) of the second Weibull distribution. Must be positive.
    ! @param Htr The transitional wave height.
    ! @return The calculated value of H1/N.
    ! Returns `Huge(1.0_dp)` if input parameters are invalid.
        REAL(KIND=dp), INTENT(IN) :: N_val, H1, H2, k1_val, k2_val, Htr
        REAL(KIND=dp) :: H_N_val, term1_a, term1_x_ln_Nval, term1_x_HtrH1, &
                             gamma_term1_part1, gamma_term1_part2, gamma_term1, &
                             term2_a, term2_x_HtrH2, gamma_term2, log_N_val, Htr_H1_ratio, Htr_H2_ratio

        IF (H1 <= 0.0_dp .OR. H2 <= 0.0_dp) THEN
            WRITE(*,*) "Error: calculate_H1N: H1 and H2 must be positive."
            calculate_H1N = HUGE(1.0_dp)
            RETURN
        END IF
        IF (k1_val <= 0.0_dp .OR. k2_val <= 0.0_dp) THEN
            WRITE(*,*) "Error: calculate_H1N: k1 and k2 must be positive."
            calculate_H1N = HUGE(1.0_dp)
            RETURN
        END IF
        IF (N_val <= 1.0_dp) THEN
            WRITE(*,*) "Error: calculate_H1N: N_val must be greater than 1."
            calculate_H1N = HUGE(1.0_dp)
            RETURN
        END IF

        H_N_val = calculate_HN(N_val, H1, H2, k1_val, k2_val, Htr)

        term1_a = 1.0_dp / k1_val + 1.0_dp
        log_N_val = LOG(N_val)
        term1_x_ln_Nval = log_N_val

        Htr_H1_ratio = Htr / H1
        term1_x_HtrH1 = (Htr_H1_ratio)**k1_val

        term2_a = 1.0_dp / k2_val + 1.0_dp
        Htr_H2_ratio = Htr / H2
        term2_x_HtrH2 = (Htr_H2_ratio)**k2_val

        IF (H_N_val < Htr - EPSILON) THEN
        ! Contribution from the first Weibull distribution.
            gamma_term1_part1 = incomplete_gamma_upper(term1_a, term1_x_ln_Nval)
            gamma_term1_part2 = incomplete_gamma_upper(term1_a, term1_x_HtrH1)
            gamma_term1 = gamma_term1_part1 - gamma_term1_part2

        ! Contribution from the second Weibull distribution.
            gamma_term2 = incomplete_gamma_upper(term2_a, term2_x_HtrH2)

            calculate_H1N = N_val * H1 * gamma_term1 + N_val * H2 * gamma_term2
        ELSE
        ! Case 2: H_N_val is greater than or equal to Htr.
        ! This means the integration for H1/N only involves the second part of the Composite Weibull distribution.
            calculate_H1N = N_val * H2 * incomplete_gamma_upper(term2_a, term1_x_ln_Nval)
        END IF
    END FUNCTION calculate_H1N

    REAL(KIND=dp) FUNCTION F1(H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val)
        IMPLICIT NONE
    ! @brief Defines the first non-linear equation F1(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
    !
    ! This equation represents the normalized Hrms constraint for the Composite Weibull distribution.
    ! It states that the overall normalized Hrms of the Composite Weibull distribution must precisely equal 1.
    ! The equation directly uses the unnormalized lower incomplete gamma function, γ(a,x),
    ! and the unnormalized upper incomplete gamma function, Γ(a,x).
    !
    ! @param H1_Hrms_val The normalized scale parameter of the first Weibull distribution.
    ! @param H2_Hrms_val The normalized scale parameter of the second Weibull distribution.
    ! @param Htr_Hrms_val The normalized transitional wave height (constant for a given solve).
    ! @return The value of the first function, which should be driven to zero.
        REAL(KIND=dp), INTENT(IN) :: H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val
        REAL(KIND=dp) :: u, v, term1, term2, sum_terms

        IF (H1_Hrms_val <= 0.0_dp .OR. H2_Hrms_val <= 0.0_dp) THEN
            F1 = HUGE(1.0_dp)
            RETURN
        END IF

        u = (Htr_Hrms_val / H1_Hrms_val)**k1
        v = (Htr_Hrms_val / H2_Hrms_val)**k2

    ! Calculate terms using unnormalized incomplete gamma functions.
        term1 = H1_Hrms_val * H1_Hrms_val * incomplete_gamma_lower(2.0_dp / k1 + 1.0_dp, u)
        term2 = H2_Hrms_val * H2_Hrms_val * incomplete_gamma_upper(2.0_dp / k2 + 1.0_dp, v)
        sum_terms = term1 + term2

        F1 = SQRT(sum_terms) - 1.0_dp
    END FUNCTION F1


    REAL(KIND=dp) FUNCTION F2(H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val)
        IMPLICIT NONE
    ! @brief Defines the second non-linear equation F2(H1_Hrms, H2_Hrms, Htr_Hrms) = 0.
    !
    ! This equation represents the continuity condition between the two Weibull distributions
    ! at the transitional wave height Htr: `(Htr/H1)^k1 = (Htr/H2)^k2`.
    !
    ! @param H1_Hrms_val The normalized scale parameter of the first Weibull distribution.
    ! @param H2_Hrms_val The normalized scale parameter of the second Weibull distribution.
    ! @param Htr_Hrms_val The normalized transitional wave height (constant for a given solve).
    ! @return The value of the second function, which should be driven to zero.
        REAL(KIND=dp), INTENT(IN) :: H1_Hrms_val, H2_Hrms_val, Htr_Hrms_val
        REAL(KIND=dp) :: term1_F2, term2_F2

        IF (H1_Hrms_val <= 0.0_dp .OR. H2_Hrms_val <= 0.0_dp) THEN
        ! Return a large value to push the solver away from invalid regions.
            F2 = HUGE(1.0_dp)
            RETURN
        END IF

        term1_F2 = (Htr_Hrms_val / H1_Hrms_val)**k1
        term2_F2 = (Htr_Hrms_val / H2_Hrms_val)**k2

        F2 = term1_F2 - term2_F2
    END FUNCTION F2

    SUBROUTINE solve_linear_system_2x2(J11, J12, J21, J22, b1, b2, dx1, dx2)
        IMPLICIT NONE
    ! @brief Solves a 2x2 linear system Ax = b for x using Cramer's rule.
    !
    ! This function is a helper for the Newton-Raphson method for systems.
    ! It takes the Jacobian matrix elements (J11, J12, J21, J22) and the negative
    ! function values (-F1, -F2) as the right-hand side, and computes the updates
    ! (dx1, dx2) for H1_Hrms and H2_Hrms.
    !
    ! @param J11 Element (1,1) of the Jacobian matrix (dF1/dH1).
    ! @param J12 Element (1,2) of the Jacobian matrix (dF1/dH2).
    ! @param J21 Element (2,1) of the Jacobian matrix (dF2/dH1).
    ! @param J22 Element (2,2) of the Jacobian matrix (dF2/dH2).
    ! @param b1 Right-hand side for the first equation (-F1).
    ! @param b2 Right-hand side for the second equation (-F2).
    ! @param dx1 Output: The calculated change for H1_Hrms.
    ! @param dx2 Output: The calculated change for H2_Hrms.
        REAL(KIND=dp), INTENT(IN) :: J11, J12, J21, J22, b1, b2
        REAL(KIND=dp), INTENT(OUT) :: dx1, dx2
        REAL(KIND=dp) :: determinant

        determinant = J11 * J22 - J12 * J21

    ! Check for a singular or nearly singular Jacobian matrix.
        IF (ABS(determinant) < TINY(1.0_dp) * 100.0_dp) THEN
            dx1 = 0.0_dp
            dx2 = 0.0_dp
            RETURN
        END IF

    ! Apply Cramer's rule:
        dx1 = (b1 * J22 - b2 * J12) / determinant
        dx2 = (J11 * b2 - J21 * b1) / determinant
    END SUBROUTINE solve_linear_system_2x2

    SUBROUTINE calculate_initial_H_guesses_from_Htr(Htr_Hrms_val, H1_initial, H2_initial)
        IMPLICIT NONE
    ! @brief Provides initial guesses for H1_Hrms and H2_Hrms based on a known Htr_Hrms value.
    !
    ! A good initial guess is crucial for the efficiency and robustness of the Newton-Raphson method.
    ! This function uses an empirical regression for H1_Hrms and H2_Hrms based on the pre-calculated,
    ! deterministic value of Htr_Hrms.
    !
    ! @param Htr_Hrms_val The known, normalized transitional wave height.
    ! @param H1_initial Output: The initial guess for H1_Hrms.
    ! @param H2_initial Output: The initial guess for H2_Hrms.
        REAL(KIND=dp), INTENT(IN) :: Htr_Hrms_val
        REAL(KIND=dp), INTENT(OUT) :: H1_initial, H2_initial
        REAL(KIND=dp) :: exp_term, Htr_Hrms_val_power

    ! Empirical regression for H1_initial
        exp_term = EXP(-1.42537392576977_dp * Htr_Hrms_val)
        H1_initial = 0.9552427998926_dp / (1.0_dp - 0.992405988921401_dp * exp_term)

    ! Empirical regression for H2_initial
        Htr_Hrms_val_power = (Htr_Hrms_val**2.980718327103574_dp)
        H2_initial = 1.054085273232950_dp + 0.9369023639428842_dp * Htr_Hrms_val_power / &
                     ((2.549022900471753_dp**2.980718327103574_dp) + Htr_Hrms_val_power)

    END SUBROUTINE calculate_initial_H_guesses_from_Htr

    LOGICAL FUNCTION newtonRaphsonSystemSolver(Htr_Hrms_val, H1_Hrms_out, H2_Hrms_out, tol, maxit)
        IMPLICIT NONE
    ! @brief Solves for H1_Hrms and H2_Hrms simultaneously using the Newton-Raphson method for systems.
    !
    ! This function implements the multi-dimensional Newton-Raphson algorithm to find the roots
    ! of the system of non-linear equations F1 and F2. It iteratively refines the guesses
    ! for H1_Hrms and H2_Hrms until the functions F1 and F2 are sufficiently close to zero.
    !
    ! @param Htr_Hrms_val The normalized transitional wave height (constant for this solve).
    ! @param H1_Hrms_out Output: The converged normalized scale parameter of the first Weibull distribution.
    ! @param H2_Hrms_out Output: The converged normalized scale parameter of the second Weibull distribution.
    ! @param tol The desired tolerance for convergence (maximum absolute value of F1 and F2).
    ! @param maxit The maximum number of iterations allowed.
    ! @return `.TRUE.` if the solver successfully converges, `.FALSE.` otherwise.
        REAL(KIND=dp), INTENT(IN) :: Htr_Hrms_val, tol
        REAL(KIND=dp), INTENT(OUT) :: H1_Hrms_out, H2_Hrms_out
        INTEGER, INTENT(IN) :: maxit
        REAL(KIND=dp) :: f1_val, f2_val
        REAL(KIND=dp) :: J11, J12, J21, J22
        REAL(KIND=dp) :: dH1, dH2
        INTEGER :: iter

        CALL calculate_initial_H_guesses_from_Htr(Htr_Hrms_val, H1_Hrms_out, H2_Hrms_out)

        DO iter = 0, maxit - 1
            f1_val = F1(H1_Hrms_out, H2_Hrms_out, Htr_Hrms_val)
            f2_val = F2(H1_Hrms_out, H2_Hrms_out, Htr_Hrms_val)

            IF (ABS(f1_val) < tol .AND. ABS(f2_val) < tol) THEN
                newtonRaphsonSystemSolver = .TRUE.
                RETURN
            END IF

        ! Calculate the Jacobian matrix elements using central finite differences.
        ! J11 = dF1/dH1_Hrms_out
            J11 = (F1(H1_Hrms_out + JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val) - &
                   F1(H1_Hrms_out - JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val)) / (2.0_dp * JACOBIAN_DX)

        ! J12 = dF1/dH2_Hrms_out
            J12 = (F1(H1_Hrms_out, H2_Hrms_out + JACOBIAN_DX, Htr_Hrms_val) - &
                   F1(H1_Hrms_out, H2_Hrms_out - JACOBIAN_DX, Htr_Hrms_val)) / (2.0_dp * JACOBIAN_DX)

        ! J21 = dF2/dH1_Hrms_out
            J21 = (F2(H1_Hrms_out + JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val) - &
                   F2(H1_Hrms_out - JACOBIAN_DX, H2_Hrms_out, Htr_Hrms_val)) / (2.0_dp * JACOBIAN_DX)

        ! J22 = dF2/dH2_Hrms_out
            J22 = (F2(H1_Hrms_out, H2_Hrms_out + JACOBIAN_DX, Htr_Hrms_val) - &
                   F2(H1_Hrms_out, H2_Hrms_out - JACOBIAN_DX, Htr_Hrms_val)) / (2.0_dp * JACOBIAN_DX)

            CALL solve_linear_system_2x2(J11, J12, J21, J22, -f1_val, -f2_val, dH1, dH2)

            H1_Hrms_out = H1_Hrms_out + dH1
            H2_Hrms_out = H2_Hrms_out + dH2

        ! Ensure H1_Hrms_out and H2_Hrms_out remain positive
            IF (H1_Hrms_out <= 0.0_dp) H1_Hrms_out = TINY(1.0_dp)
            IF (H2_Hrms_out <= 0.0_dp) H2_Hrms_out = TINY(1.0_dp)

        END DO

        WRITE(*,*) "Newton-Raphson system solver failed to converge for Htr_Hrms =", Htr_Hrms_val, &
             " after ", maxit, " iterations."
        newtonRaphsonSystemSolver = .FALSE.
    END FUNCTION newtonRaphsonSystemSolver

!---------------------------------------------------------------------
! Subroutine: buildReportAndDisplay
! Purpose:
!   Computes all wave parameters, performs calculations based on the
!   Composed Weibull model, and builds a detailed report. It then
!   writes this report to "report.txt".
!---------------------------------------------------------------------
    SUBROUTINE buildReportAndDisplay(Hm0, d, slopeM)
        IMPLICIT NONE
        REAL(KIND=dp), INTENT(IN) :: Hm0, d, slopeM
        CHARACTER(LEN=4000) :: report_str! A large enough string to hold the report
        CHARACTER(LEN=100) :: temp_line! Temporary string for building each line
        REAL(KIND=dp) :: m0, Hrms, tanAlpha, Htr_dim, Htr_tilde
        REAL(KIND=dp) :: H1_Hrms, H2_Hrms
        INTEGER, PARAMETER :: N_values_for_H1N(6) = (/3, 10, 50, 100, 250, 1000/)
        REAL(KIND=dp) :: interp_Hrms(2 + SIZE(N_values_for_H1N))! For H1/Hrms, H2/Hrms, and H1/N values
        REAL(KIND=dp) :: dimensional(2 + SIZE(N_values_for_H1N))
        REAL(KIND=dp) :: ratio_110_13, ratio_150_13, ratio_1100_13, ratio_1250_13, ratio_11000_13
        INTEGER :: j, i_temp! Renamed i to i_temp to avoid conflict with implicit type if i was used elsewhere

    ! Initialize report_str to empty
        report_str = ''

    ! Compute free-surface variance: m0 = (Hm0 / 4)^2
        m0 = (Hm0 / 4.0_dp)**2

    ! Mean square wave height: Hrms = (2.69 + 3.24*sqrt(m0)/d)*sqrt(m0)
        Hrms = (2.69_dp + 3.24_dp * SQRT(m0) / d) * SQRT(m0)

    ! Compute the actual tangent of the beach slope: tan(alpha) = 1.0 / slopeM.
        tanAlpha = 1.0_dp / slopeM

    ! Dimensional transitional wave height: Htr = (0.35 + 5.8*(1/m)) * d.
        Htr_dim = (0.35_dp + 5.8_dp * tanAlpha) * d

    ! Dimensionless transitional parameter: H̃_tr = Htr / Hrms.
        Htr_tilde = (Htr_dim / Hrms)

    ! Solve for H1_Hrms and H2_Hrms simultaneously using the Newton-Raphson matrix solver.
        IF (.NOT. newtonRaphsonSystemSolver(Htr_tilde, H1_Hrms, H2_Hrms, EPSILON, 500)) THEN
            WRITE(*,*) "ERROR: Solver failed to converge for Htr_tilde =", Htr_tilde, &
                         ". Cannot generate full report."
            report_str = "ERROR: Solver failed to converge. Cannot generate full report." // NEW_LINE('A')
            CALL writeReportToFile(report_str)
            RETURN
        END IF

    ! Store H1/Hrms and H2/Hrms directly at the beginning of the vector.
        interp_Hrms(1) = H1_Hrms
        interp_Hrms(2) = H2_Hrms

    ! Calculate H1/N values using the `calculate_H1N` function.
        DO j = 1, SIZE(N_values_for_H1N)
            interp_Hrms(j + 2) = calculate_H1N(REAL(N_values_for_H1N(j), KIND=dp), H1_Hrms, H2_Hrms, k1, k2, Htr_tilde)
        END DO

    ! Convert dimensionless values to dimensional wave heights (in meters)
        DO i_temp = 1, SIZE(interp_Hrms)
            dimensional(i_temp) = interp_Hrms(i_temp) * Hrms
        END DO

    ! Calculate diagnostic ratios.
        ratio_110_13 = (interp_Hrms(4) / interp_Hrms(3))
        ratio_150_13 = (interp_Hrms(5) / interp_Hrms(3))
        ratio_1100_13 = (interp_Hrms(6) / interp_Hrms(3))
        ratio_1250_13 = (interp_Hrms(7) / interp_Hrms(3))
        ratio_11000_13 = (interp_Hrms(8) / interp_Hrms(3))

    ! Build the report string with formatted output using internal write
        WRITE(temp_line, '(A)') "======================"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   INPUT PARAMETERS"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "======================"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        ! Adjusted to include colon and match spacing
        WRITE(temp_line, '(A17, F7.4)') "Hm0 (m)         :", Hm0
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A17, F8.4)') "d (m)           :", d
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        ! Adjusted to include colon and match spacing
        WRITE(temp_line, '(A17, F8.4, A17, F6.4, A1)') "Beach slope (m) :", slopeM, "   (tan(alpha) = ", tanAlpha, ")"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        report_str = TRIM(report_str) // NEW_LINE('A')! Blank line

    ! Section: CALCULATED PARAMETERS
        WRITE(temp_line, '(A)') "==========================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   CALCULATED PARAMETERS"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "==========================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        ! Adjusted to include colon and match spacing
        WRITE(temp_line, '(A34, F7.4)') "Free surface variance m0         :", m0
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A34, F7.4)') "Mean square wave height Hrms (m) :", Hrms
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A34, F7.4)') "Dimensionless H~_tr (Htr/Hrms)   :", Htr_tilde
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A34, F7.4)') "Transitional wave height Htr (m) :", Htr_dim
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        report_str = TRIM(report_str) // NEW_LINE('A')! Blank line

    ! Section: DIMENSIONLESS WAVE HEIGHTS (H/Hrms)
        WRITE(temp_line, '(A)') "========================================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   DIMENSIONLESS WAVE HEIGHTS (H/Hrms)"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "========================================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/Hrms       : ", interp_Hrms(1)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H2/Hrms       : ", interp_Hrms(2)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/3 / Hrms   : ", interp_Hrms(3)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/10 / Hrms  : ", interp_Hrms(4)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/50 / Hrms  : ", interp_Hrms(5)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/100 / Hrms : ", interp_Hrms(6)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/250 / Hrms : ", interp_Hrms(7)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/1000 /Hrms : ", interp_Hrms(8)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        report_str = TRIM(report_str) // NEW_LINE('A')! Blank line

    ! Section: DIMENSIONAL WAVE HEIGHTS (m)
        WRITE(temp_line, '(A)') "=================================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   DIMENSIONAL WAVE HEIGHTS (m)"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "=================================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1 (m)        : ", dimensional(1)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H2 (m)        : ", dimensional(2)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/3 (m)      : ", dimensional(3)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/10 (m)     : ", dimensional(4)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/50 (m)     : ", dimensional(5)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/100 (m)    : ", dimensional(6)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/250 (m)    : ", dimensional(7)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A15, F7.4)') "H1/1000 (m)   : ", dimensional(8)
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        report_str = TRIM(report_str) // NEW_LINE('A')! Blank line

    ! Section: DIAGNOSTIC RATIOS
        WRITE(temp_line, '(A)') "======================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "   DIAGNOSTIC RATIOS"
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')
        WRITE(temp_line, '(A)') "======================="
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        ! Adjusted to include colon and match spacing
        WRITE(temp_line, '(A18, F7.4)') "(H1/10)/(H1/3)   :", ratio_110_13
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A18, F7.4)') "(H1/50)/(H1/3)   :", ratio_150_13
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A18, F7.4)') "(H1/100)/(H1/3)  :", ratio_1100_13
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A18, F7.4)') "(H1/250)/(H1/3)  :", ratio_1250_13
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A18, F7.4)') "(H1/1000)/(H1/3) :", ratio_11000_13
        report_str = TRIM(report_str) // TRIM(temp_line) // NEW_LINE('A')

        WRITE(temp_line, '(A)') "End of Report"
        report_str = TRIM(report_str) // TRIM(temp_line) ! Removed NEW_LINE('A')

        CALL writeReportToFile(TRIM(report_str))
        WRITE(*,'(A)') TRIM(report_str)
    END SUBROUTINE buildReportAndDisplay

!---------------------------------------------------------------------
! Subroutine: writeReportToFile
! Purpose:
!   Writes the generated report string to a file named "report.txt".
!---------------------------------------------------------------------
    SUBROUTINE writeReportToFile(report)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: report
        INTEGER :: iostat_val

        OPEN(UNIT=REPORT_UNIT, FILE="report.txt", STATUS="REPLACE", ACTION="WRITE", IOSTAT=iostat_val)
        IF (iostat_val /= 0) THEN
            WRITE(*,*) "Error: Could not open report.txt for writing."
            RETURN
        END IF
        WRITE(REPORT_UNIT, '(A)') report
        CLOSE(REPORT_UNIT)
    END SUBROUTINE writeReportToFile

END PROGRAM ShallowWaterWaves
