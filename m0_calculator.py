# -*- coding: utf-8 -*-
"""
===============================================================================
 COMPREHENSIVE WAVE PARAMETER CALCULATOR
===============================================================================

 SCIENTIFIC AND PHYSICAL MODEL DESCRIPTION:
 ------------------------------------------
 This calculator provides a robust estimation of the free-surface variance (m0)
 of ocean waves in the nearshore zone. The free-surface variance, or zeroth
 spectral moment, is a fundamental measure of the total energy in a sea state.
 Its accurate prediction in shallow water is critical for coastal engineering,
 sediment transport, and morphological modeling.

 The model implements a hybrid methodology, transitioning between two distinct
 physical regimes based on the local wave environment, which is primarily
 dictated by the ratio of water depth to wavelength (d/L).

 1. LINEAR, ENERGY-CONSERVING REGIME (Deeper Water):
    Outside the surf zone, where wave breaking is not the dominant process,
    the model assumes energy conservation. Wave transformation is governed by
    linear wave theory. The change in wave height is calculated through the
    principle of shoaling, derived from the conservation of energy flux
    (d/dx(E*Cg) = 0). In this regime, the sea surface is assumed to be a
    Gaussian process, and consequently, the wave heights are described by a
    Rayleigh distribution. This leads to the well-known relationships:
    - Hm0 = 4 * sqrt(m0)
    - Hrms = sqrt(8 * m0)

 2. NON-LINEAR, DISSIPATIVE REGIME (Surf Zone):
    As waves propagate into shallow water (typically d < L/20), they steepen
    and break, leading to significant energy dissipation. The assumptions of
    linear theory and simple energy conservation are no longer valid. In this
    regime, the calculator employs a more sophisticated two-step approach
    grounded in established surf zone physics:

    a) Hrms Calculation via Energy Flux Balance with Dissipation:
       The local root-mean-square wave height (Hrms) is calculated using an
       energy dissipation model based on the work of Thornton & Guza (1983).
       This model extends the energy flux balance equation to include a term
       for energy dissipation due to breaking: d/dx(E*Cg) = -<epsilon_b>.
       The dissipation rate, <epsilon_b>, is modeled by analogy to a periodic
       bore, where dissipation is proportional to H^3. By integrating this
       balance from an offshore boundary, the model predicts the decay of
       Hrms across the surf zone.

    b) m0 Calculation via Advanced Statistical Model:
       Within the surf zone, depth-induced breaking clips the upper tail of
       the wave height distribution, causing it to deviate from the Rayleigh
       distribution. To account for this, the model uses the empirical
       relationship proposed by Battjes & Groenendijk (2000). This formula
       was derived from extensive laboratory data and provides a more accurate
       link between Hrms, m0, and local depth (d) in a breaking-dominated
       environment. The calculator solves this empirical formula to find the
       local m0 from the Hrms value calculated in the previous step.

 This hybrid approach ensures that the most appropriate physical models are
 applied, providing a more accurate and physically realistic estimation of
 wave parameters across the entire nearshore profile.

 KEY LITERATURE:
 ---------------------
 - Battjes, J. A., & Groenendijk, H. W. (2000). Wave height distributions on
   shallow foreshores. Coastal Engineering, 40(3), 161-182.
   (Provides the empirical Hrms-m0 relationship for the surf zone).

 - Thornton, E. B., & Guza, R. T. (1983). Transformation of wave height
   distribution. Journal of Geophysical Research, 88(C10), 5925-5938.
   (Provides the foundational model for wave energy dissipation used to
   calculate Hrms in the surf zone).

 - Thornton, E. B., & Guza, R. T. (1982). Energy saturation and phase speeds
   measured on a natural beach. Journal of Geophysical Research, 87(C12),
   9499-9508.
   (Provides the physical basis for the concept of energy saturation, where
   Hrms becomes proportional to local depth in the inner surf zone).
"""

import math
import sys

# NEW: A helper class to direct print statements to multiple outputs.
class Logger(object):
    """
    A file-like object that writes to both a terminal and a log file.
    """
    def __init__(self, terminal, logfile):
        self.terminal = terminal
        self.logfile = logfile

    def write(self, message):
        self.terminal.write(message)
        self.logfile.write(message)

    def flush(self):
        # This flush method is needed for compatibility with Python 3.
        self.terminal.flush()
        self.logfile.flush()

def get_validated_input(prompt, validation_func=None):
    """
    A helper function to get and validate numeric input from the user.
    It shows the prompt in the console and writes the input value to a file.

    Args:
        prompt (str): The message to display to the user.
        validation_func (function, optional): A function that takes the user's
                                              input as an argument and returns
                                              True if it's valid, False otherwise.
                                              Defaults to None.

    Returns:
        float: The validated numeric value entered by the user.
    """
    while True:
        try:
            # Write the prompt directly to the console and flush to ensure it's visible.
            # We use __stdout__ to bypass our Logger and ensure this only goes to the console.
            sys.__stdout__.write(prompt)
            sys.__stdout__.flush()

            # Read the line directly from standard input.
            value_str = sys.stdin.readline()
            # Convert the input string to a floating-point number.
            value = float(value_str)

            # The repetitive print statement that echoed user input has been removed.
            
            # Check if the input is valid using the provided validation function.
            if validation_func is None or validation_func(value):
                # We also write the value to the log file here so it's recorded.
                sys.stdout.logfile.write(f"{prompt.strip()} {value}\n")
                return value
            else:
                # If validation fails, notify the user on the console.
                sys.stderr.write("Input is not valid (e.g., must be > 0). Please try again.\n")

        except (ValueError, TypeError):
            # Handle cases where the input is not a valid number.
            sys.stderr.write("Invalid input. Please enter a numeric value.\n")
        except Exception as e:
            # Handle any other unexpected errors during input.
            sys.stderr.write(f"An unexpected error occurred: {e}\n")

def shoaling_coefficient(k, depth, tolerance=1e-9):
    """
    Calculates the shoaling coefficient (Ks) based on linear wave theory.

    Physical Justification:
    The shoaling coefficient quantifies the change in wave height due to
    variations in water depth, assuming no energy dissipation. It is derived
    from the principle of conservation of energy flux (P = E * Cg = constant).
    The derivation leads to: Ks = H / H_deep = sqrt(Cg_deep / Cg_local).
    This function computes Ks using the local group velocity (Cg_local).

    Args:
        k (float): The local wave number (rad/m), defined as 2 * pi / L.
        depth (float): The local water depth (m).
        tolerance (float, optional): A small number to prevent division by zero.

    Returns:
        float: The dimensionless shoaling coefficient (Ks).
    """
    if k <= 0.0 or depth <= 0.0:
        return 1.0

    kh = k * depth

    try:
        # The group velocity Cg is given by n*C, where n = 0.5 * (1 + 2*kh / sinh(2*kh)).
        # Cg_deep is 0.5 * C_deep.
        # Ks^2 = Cg_deep / Cg = (0.5 * C_deep) / (n * C)
        # Using C/C_deep = tanh(kh), this simplifies to:
        # Ks^2 = 1 / (2 * n * tanh(kh))
        # Ks^2 = 1 / ( (1 + 2*kh / sinh(2*kh)) * tanh(kh) )
        denominator = math.tanh(kh) * (1.0 + 2.0 * kh / math.sinh(2.0 * kh))

        if denominator < tolerance:
            return 1.0

        return 1.0 / math.sqrt(denominator)

    except OverflowError:
        # For very large kh (deep water), sinh(2*kh) overflows. In this physical
        # limit, Cg approaches Cg_deep, so the ratio approaches 1, and Ks is 1.
        return 1.0

def solve_for_L_newton(L0, d, tolerance=1e-7, max_iterations=100):
    """
    Solves the linear wave theory dispersion relation for the local wavelength (L)
    using the Newton-Raphson iterative method.

    Physical Justification:
    The dispersion relation, omega^2 = g * k * tanh(k*h), is the fundamental
    equation of linear wave theory. It links a wave's temporal property (period T,
    via omega = 2*pi/T) to its spatial property (wavelength L, via k = 2*pi/L)
    at a given water depth (d). This equation is transcendental and must be solved
    numerically. This function solves the form: L = L0 * tanh(2*pi*d/L).

    Args:
        L0 (float): The deep-water wavelength (m), given by g*T^2/(2*pi).
        d (float): The local water depth (m).
        tolerance (float, optional): The convergence tolerance for the solution.
        max_iterations (int, optional): The maximum number of iterations.

    Returns:
        float: The calculated local wavelength (L) in meters.
    """
    # Define the function f(L) = 0 that we want to solve: f(L) = L - L0*tanh(2*pi*d/L)
    def f(L):
        if L <= 0: return float('inf')
        return L - L0 * math.tanh(2 * math.pi * d / L)

    # Define the derivative of f(L) for the Newton-Raphson method.
    def f_prime(L):
        if L <= 0: return float('inf')
        try:
            cosh_val = math.cosh(2 * math.pi * d / L)
            # f'(L) = 1 + (2*pi*d*L0/L^2) * sech^2(2*pi*d/L)
            return 1 + (2 * math.pi * d * L0 / (L**2)) / (cosh_val**2)
        except OverflowError:
            # For large arguments, cosh overflows, f_prime approaches 1.
            return 1.0

    try:
        # Start with a good initial guess for L. The Eckart (1952) approximation
        # provides a robust starting point that is explicit and reasonably accurate.
        # L_eckart = L0 * sqrt(tanh((2*pi/L0)^2 * d * (g*T^2/(2*pi)))
        # This simplifies to L_eckart = L0 * tanh( (2*pi*d/L0)^0.5 )
        # Note: The original code had a slight variation, this is a common form.
        L_current = L0 * math.tanh((2 * math.pi * d / L0)**0.5)
    except (ValueError, OverflowError):
        # Fallback to deep water wavelength if guess fails.
        L_current = L0

    # Iterate using Newton-Raphson: L_next = L_current - f(L)/f'(L)
    for i in range(max_iterations):
        f_val = f(L_current)
        f_prime_val = f_prime(L_current)

        # Avoid division by zero if derivative is flat.
        if abs(f_prime_val) < 1e-10:
            break

        L_next = L_current - f_val / f_prime_val

        # Check for convergence.
        if abs(L_next - L_current) < tolerance:
            return L_next

        L_current = L_next

    return L_current

def calculate_m0_with_full_TG_formula(Hs0, d, d0_boundary):
    """
    Calculates m0 in the saturated surf zone using the hybrid T&G/B&G method.

    Physical Justification:
    This function models the physics of a breaking-dominated surf zone. It
    combines two key theories:

    1. Hrms Calculation (Thornton & Guza, 1983):
       It first calculates the local root-mean-square wave height (Hrms) by
       solving the energy flux balance equation that includes a term for
       energy dissipation due to breaking. The dissipation is modeled by
       analogy with a periodic bore, resulting in a predictive equation for
       the decay of Hrms as waves move from an offshore boundary (d0_boundary)
       to the local depth (d).

    2. m0 Calculation (Battjes & Groenendijk, 2000):
       In the surf zone, the relationship Hrms = sqrt(8*m0) (valid for a
       Rayleigh distribution) breaks down due to depth-induced breaking.
       This function uses the empirical formula from B&G (2000), which was
       derived from extensive lab data on shallow foreshores:
           Hrms = (2.69 + 3.24 * sqrt(m0)/d) * sqrt(m0)
       This equation is solved as a quadratic for sqrt(m0) to find the local
       free-surface variance (m0) from the Hrms calculated in step 1.

    Args:
        Hs0 (float): The significant wave height at the shallow water boundary (d0).
        d (float): The local water depth where m0 is being calculated.
        d0_boundary (float): The water depth at the offshore boundary of the surf zone.

    Returns:
        tuple[float, float] or tuple[None, None]: A tuple containing (m0, Hrms),
                                                  or (None, None) if an error occurs.
    """
    try:
        # --- Step 1: Calculate Hrms using Thornton & Guza (1983) model ---
        # This is an integrated form of d(E*Cg)/dx = -<epsilon_b>.
        # The parameter 'a' is an empirical coefficient related to the breaking
        # process, and 'yd' is the offshore forcing parameter.
        a = (0.51 * d**0.1)**5
        yd_offshore_forcing = Hs0**2 * math.sqrt(d0_boundary)

        if yd_offshore_forcing == 0:
            print("Error: Offshore forcing parameter 'yd' is zero.")
            return None, None

        inner_term_d0 = 1 / d0_boundary**(23/4)
        inner_term_yd = a / yd_offshore_forcing**(5/2)

        full_term = 1 - d**(23/4) * (inner_term_d0 - inner_term_yd)

        if full_term <= 0:
            # This indicates that the model predicts full dissipation before
            # reaching the target depth 'd', which can happen if the beach
            # is very dissipative or the offshore forcing is low.
            print(f"Error: Calculation resulted in a non-positive term ({full_term:.4f}) inside the root.")
            return None, None

        Hrms = a**(1/5) * d**(9/10) * (full_term)**(-1/5)

        # --- Step 2: Calculate m0 from Hrms using Battjes & Groenendijk (2000) ---
        # The B&G formula is: Hrms = (2.69 + 3.24 * sqrt(m0)/d) * sqrt(m0)
        # Let x = sqrt(m0). This gives a quadratic equation for x:
        # (3.24/d) * x^2 + 2.69 * x - Hrms = 0

        # Coefficients for the quadratic formula: ax^2 + bx + c = 0
        quad_a = 3.24 / d
        quad_b = 2.69
        quad_c = -Hrms

        # Calculate the discriminant of the quadratic equation.
        discriminant = quad_b**2 - 4 * quad_a * quad_c

        if discriminant < 0:
            # A negative discriminant means no real solution exists, which is
            # physically unrealistic and points to an issue with the inputs or
            # the applicability of the formula for the given parameters.
            print(f"Error: Negative discriminant ({discriminant:.4f}) in quadratic solver for m0.")
            return None, None

        # Solve for x = sqrt(m0). We take the positive root as sqrt(m0) must be positive.
        sqrt_m0 = (-quad_b + math.sqrt(discriminant)) / (2 * quad_a)

        # Final calculation for m0.
        m0 = sqrt_m0**2

        return m0, Hrms

    except (ValueError, ZeroDivisionError) as e:
        print(f"\nAn error occurred during the m0 calculation: {e}")
        return None, None

def comprehensive_wave_calculator():
    """
    Main function to drive the wave parameter calculations. It gathers user
    inputs, performs intermediate calculations, decides which physical model
    regime to apply, and prints a comprehensive report.
    """
    print("--- Comprehensive Wave Parameter Calculator ---")
    print("Please provide the following initial values.")

    g = 9.81  # Acceleration due to gravity (m/s^2)
    pi = math.pi

    # --- Get User Inputs ---
    Hm0_d_initial = get_validated_input("Enter Wave Height at depth d (Hm0_d): ", lambda x: x > 0)
    Tp = get_validated_input("Enter Peak Wave Period (Tp): ", lambda x: x > 0)
    d = get_validated_input("Enter Local Water Depth (d): ", lambda x: x > 0)
    m = get_validated_input("Enter Beach Slope (1:m): ", lambda x: x > 0)

    # --- Initial Calculations based on Linear Wave Theory ---
    # Deep water wavelength (L0) depends only on wave period.
    L0 = (g * Tp**2) / (2 * pi)
    # Solve the dispersion relation to find the local wavelength at depth d.
    L_local = solve_for_L_newton(L0, d)
    # Calculate the local wave number.
    k_local = 2 * pi / L_local if L_local > 0 else float('inf')

    # --- Wave Breaking Check (Miche Criterion) at Local Depth (d) ---
    # Miche (1944) proposed a theoretical limit for wave steepness (H/L) before
    # breaking. This is a geometric constraint.
    # [H/L]_max = 0.142 * tanh(k*d)
    Hb_Miche_d = 0.142 * L_local * math.tanh(k_local * d)
    is_capped_d = False
    Hm0_d = Hm0_d_initial
    if Hm0_d > Hb_Miche_d:
        print(f"\n(!) NOTE: Initial wave height {Hm0_d:.4f}m exceeds breaking limit {Hb_Miche_d:.4f}m at local depth d.")
        Hm0_d = Hb_Miche_d
        print(f"    Adjusted Hm0_d to {Hm0_d:.4f}m for subsequent calculations.")
        is_capped_d = True

    # --- Back-calculate deep water parameters from local conditions ---
    # Use the shoaling coefficient to find the equivalent deep-water wave height.
    Ks_d = shoaling_coefficient(k_local, d)
    Hm0_deep = Hm0_d / Ks_d

    # --- Conditional Logic for Model Selection ---
    m0 = None
    Hrms = None
    m0_method_str = ""
    use_tg_formula = False

    L_at_d0, d0, k_at_d0, Ks_at_d0, Hm0_d0 = 0, 0, 0, 0, 0
    is_capped_d0 = False

    # Determine the physical regime. The transition to a "shallow water"
    # dissipative regime is often defined by d/L < 1/20. We use this to
    # define the boundary depth 'd0' where the surf zone model should be applied.
    # First, find the wavelength at d/L = 1/20 by solving L = L0*tanh(2*pi*(L/20)/L).
    d0_check_L = L0 * math.tanh(pi / 10)
    d0_check = d0_check_L / 20 if d0_check_L > 0 else 0

    # Check if the local depth 'd' is inside this shallow water boundary.
    if d <= d0_check:
        use_tg_formula = True
        print(f"\nCondition met: d ({d:.4f}m) <= d0_check ({d0_check:.4f}m).")
        print(f"Physical Regime: Shallow Water / Surf Zone.")
        print(f"Applying Model: Thornton & Guza (1983) / Battjes & Groenendijk (2000) hybrid.")

        # These are the wave parameters at the boundary of the surf zone.
        d0 = d0_check
        L_at_d0 = d0_check_L
        k_at_d0 = 2 * pi / L_at_d0 if L_at_d0 > 0 else float('inf')
        Ks_at_d0 = shoaling_coefficient(k_at_d0, d0)
        # Shoal the deep water wave height to the boundary depth d0.
        Hm0_d0 = Hm0_deep * Ks_at_d0

        # Perform another wave breaking check at this new shallow depth (d0).
        Hb_Miche_d0 = 0.142 * L_at_d0 * math.tanh(k_at_d0 * d0)
        if Hm0_d0 > Hb_Miche_d0:
            print(f"\n(!) NOTE: Calculated wave height at d0 ({Hm0_d0:.4f}m) exceeds breaking limit {Hb_Miche_d0:.4f}m.")
            Hm0_d0 = Hb_Miche_d0
            print(f"    Adjusted Hm0_d0 to {Hb_Miche_d0:.4f}m.")
            is_capped_d0 = True

        # Calculate m0 and Hrms using the hybrid T&G/B&G formula for the surf zone.
        m0, Hrms = calculate_m0_with_full_TG_formula(Hs0=Hm0_d0, d=d, d0_boundary=d0)

        if m0 is None:
            m0_method_str = "(T&G/B&G method failed, no result)"
        else:
            m0_method_str = "(T&G/B&G method)"
    else:
        # If outside the surf zone, use the simpler Rayleigh method.
        print(f"\nCondition for surf zone model not met (d ({d:.4f}m) > d0_check ({d0_check:.4f}m)).")
        print(f"Physical Regime: Intermediate/Deep Water.")
        print(f"Applying Model: Linear Theory with Rayleigh Statistics.")

        # For a Rayleigh distribution, m0 is related to Hm0 by Hm0 = 4*sqrt(m0).
        m0 = (Hm0_d / 4.0)**2
        # Hrms is then calculated from m0, assuming Hrms = sqrt(8 * m0).
        Hrms = math.sqrt(8 * m0)
        m0_method_str = "(Rayleigh method: (Hm0_d/4)^2)"

    if m0 is None:
        print("\nCould not complete calculations due to an error in the m0 computation.")
        return

    # The ratio Hrms/d is a key parameter in surf zone dynamics. In the inner
    # surf zone, it approaches a constant value (e.g., ~0.42 in Thornton & Guza, 1982).
    Hrms_over_d = Hrms / d if d > 0 else 0

    # --- Display Unified Report ---
    print("\n" + "="*65)
    print("--- Comprehensive Wave Parameter Report ---")
    print("="*65)

    print("\n[Input Parameters]")
    print(f"{'Initial Wave Height at d (Hm0_d)':<45}: {Hm0_d_initial:.4f} m" + (f" (capped to {Hm0_d:.4f} m)" if is_capped_d else ""))
    print(f"{'Peak Wave Period (Tp)':<45}: {Tp:.4f} s")
    print(f"{'Local Water Depth (d)':<45}: {d:.4f} m")
    print(f"{'Beach Slope (m)':<45}: {m:.4f}")

    print("\n[Derived Deep Water Parameters (d -> inf)]")
    print(f"{'Deep Water Wave Height (Hm0_deep)':<45}: {Hm0_deep:.4f} m")
    print(f"{'Deep Water Wavelength (L0)':<45}: {L0:.4f} m")

    if use_tg_formula:
        print("\n[Derived Parameters at Surf Zone Boundary (d0)]")
        print(f"{'Wave Height (Hm0_d0)':<45}: {Hm0_d0:.4f} m" + (" (breaking)" if is_capped_d0 else ""))
        print(f"{'Shallow Water Depth (d0)':<45}: {d0:.4f} m")
        print(f"{'Wavelength at d0 (L_d0)':<45}: {L_at_d0:.4f} m")
        print(f"{'Shoaling Coefficient (Ks_d0)':<45}: {Ks_at_d0:.4f}")

    print("\n[Derived Parameters at Local Depth (d)]")
    print(f"{'Local Water Depth (d)':<45}: {d:.4f} m")
    if is_capped_d:
        print(f"{'Breaking Wave Height Limit (Hb_Miche)':<45}: {Hb_Miche_d:.4f} m")
    else:
        print(f"{'Local Wave Height (Hm0_d)':<45}: {Hm0_d:.4f} m")
    if is_capped_d:
        print(f"{'Wave Height to Depth ratio (Hb_Miche/d)':<45}: {Hb_Miche_d/d:.4f}")
    else:
        print(f"{'Wave Height to Depth ratio (Hm0_d/d)':<45}: {Hm0_d/d:.4f}")
    print(f"{'Local Wavelength (L_d)':<45}: {L_local:.4f} m")
    print(f"{'Shoaling Coefficient (Ks_d)':<45}: {Ks_d:.4f}")
    if use_tg_formula:
        offshore_forcing_yd = Hm0_d0**2 * math.sqrt(d0) if d0 > 0 else 0
        print(f"{'Offshore Forcing Parameter (yd)':<45}: {offshore_forcing_yd:.4f} m^2.5")
        print(f"{'T&G Empirical Coefficient (a_eq^5)':<45}: {(0.51 * d**0.1):.4f}") # This is related to B in T&G83
    print(f"{'Free-surface variance (m0)':<45}: {m0:.4f} m^2 {m0_method_str}")
    print(f"{'Root mean square wave height (Hrms)':<45}: {Hrms:.4f} m")
    print(f"{'Hrms / d ratio (Saturation Index)':<45}: {Hrms_over_d:.4f}")
    print("\n" + "="*65)

if __name__ == '__main__':
    """
    This block executes when the script is run directly.
    It now redirects all 'print' statements to both the console and 'output.txt'.
    """
    # Keep track of the original standard output (the console).
    original_stdout = sys.stdout
    
    try:
        # Open the output file in write mode.
        with open('output.txt', 'w') as f:
            # Create an instance of our Logger class, passing it the console and the file.
            sys.stdout = Logger(original_stdout, f)

            # Run the main calculator function. All print statements inside it
            # will now go to both the console and the file.
            comprehensive_wave_calculator()

            # Restore the original standard output.
            sys.stdout = original_stdout

        print("\nCalculation complete. Output has been displayed and saved to 'output.txt'.")
    except KeyboardInterrupt:
        sys.stderr.write("\n\nProgram interrupted by user. Exiting.\n")
    except Exception as e:
        # Restore stdout before printing the error to avoid potential issues.
        sys.stdout = original_stdout
        sys.stderr.write(f"\nA critical error occurred: {e}\n")
