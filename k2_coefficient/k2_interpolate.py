import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.stats import linregress # Import for linear regression
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages # Import PdfPages

def f(target_sat: float, target_slope: float, data_df: pd.DataFrame, plot_results: bool = True) -> float:
    """
    Interpolates the k2 value for a given saturation (sat) and slope based on a two-step process.

    This interpolation method is inspired by approaches used in wave height statistics,
    particularly in the context of the Composite Weibull distribution, where parameters
    like k2 are related to physical wave characteristics.

    Reference:
    * **Groenendijk, H. W., & Van Gent, M. R. A. (1998).** *Shallow foreshore wave height statistics;
      A predictive model for the probability of exceedance of wave heights*. Technical Report H3351,
      WL | delft hydraulics, The Netherlands. <http://dx.doi.org/10.13140/RG.2.2.14180.68486>

    The process proceeds in two steps:
    1. For each unique slope in the provided dataset, it calculates a linear regression line
       (k2 vs sat) using all available data points for that slope. It then uses this
       regression line to predict k2 for the `target_sat` value. This creates a set of
       (slope, k2_predicted_at_target_sat) pairs. This step allows for extrapolation
       of 'k2' with respect to 'sat' for each individual slope.
    2. It then interpolates k2 for the `target_slope` value using PCHIP (Piecewise Cubic Hermite
       Interpolating Polynomial) interpolation on the (slope, k2_predicted_at_target_sat)
       pairs generated in the first step. PCHIP is chosen for its shape-preserving properties,
       ensuring that the interpolated curve does not overshoot or undershoot the data points.
       However, PCHIP does not extrapolate, so if `target_slope` is outside the range of
       available slopes in the dataset, a ValueError will be raised.

    Args:
        target_sat (float): The saturation value (e.g., Psi = sqrt(m0)/d) for which to interpolate k2.
                            This is the x-axis variable in the first interpolation step.
        target_slope (float): The foreshore slope value (e.g., 1:X) for which to interpolate k2.
                              This is the x-axis variable in the second interpolation step.
        data_df (pd.DataFrame): A DataFrame containing 'sat', 'slope', and 'k2' columns.
                                This DataFrame should contain the raw data points for interpolation.
                                It is crucial that this data accurately represents the relationships
                                between these parameters.
        plot_results (bool): If True, generates plots of the interpolation steps. These plots
                             are saved to a PDF if the main script is configured to do so.

    Returns:
        float: The interpolated k2 value corresponding to the `target_sat` and `target_slope`.

    Raises:
        ValueError:
            - If there are not enough unique points (less than 2) for linear regression for a given slope.
            - If there are less than 2 unique 'slope' values after the first interpolation step,
              which prevents PCHIP interpolation.
            - If the `target_slope` is outside the range of available slopes for PCHIP interpolation
              (as PCHIP does not extrapolate).
            - If no valid interpolated k2 values could be generated across slopes for the given target_sat.
        KeyError: If the input DataFrame does not contain the required columns ('sat', 'slope', 'k2').
    """
    # --- Input Validation ---
    # Ensure the input DataFrame contains all necessary columns.
    if not all(col in data_df.columns for col in ['sat', 'slope', 'k2']):
        raise KeyError("Input DataFrame must contain 'sat', 'slope', and 'k2' columns.")

    # Get unique slope values from the dataset and sort them.
    unique_slopes = sorted(data_df['slope'].unique())
    # List to store intermediate k2 values predicted for the target_sat at each unique slope.
    interpolated_k2_for_slopes = []

    # --- Consistent X-axis Scaling for Saturation Plot ---
    # Calculate the overall minimum and maximum 'sat' values from the entire dataset.
    overall_min_sat = data_df['sat'].min()
    overall_max_sat = data_df['sat'].max()
    # Determine the x-axis limits for the saturation plot.
    # A small buffer is added to ensure the target_sat is visible even if it's an extrapolation.
    x_axis_min_limit = min(overall_min_sat, target_sat) - 0.01
    x_axis_max_limit = max(overall_max_sat, target_sat) + 0.01

    # --- Color Mapping for Plots ---
    # Retrieve the default color cycle from Matplotlib's style parameters.
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    # Create a dictionary to map each unique slope value to a consistent color from the cycle.
    # This ensures that data points and their corresponding regression lines for the same slope
    # have matching colors.
    color_map = {slope: colors[i % len(colors)] for i, slope in enumerate(unique_slopes)}

    # Initialize figure object to None. It will be created only if plot_results is True.
    fig = None
    if plot_results:
        # --- Figure Setup for A3 Landscape PDF ---
        # Define A3 landscape dimensions in inches (1 inch = 25.4 mm).
        # A3 paper size is 420mm x 297mm. For landscape, width is 420mm, height is 297mm.
        # Conversion: 420/25.4 = 16.53 inches, 297/25.4 = 11.69 inches.
        fig = plt.figure(figsize=(16.53, 11.69)) # Set figure size for A3 landscape.

        # --- Subplot 1: k2 vs Saturation (Linear Regression) ---
        plt.subplot(1, 2, 1) # Create the first subplot (1 row, 2 columns, first plot).
        plt.title('Step 1: k2 vs Saturation (Linear Regression)', fontsize=16) # Increased font size
        plt.xlabel('Saturation (sat)', fontsize=14) # Increased font size
        plt.ylabel('k2', fontsize=14) # Increased font size
        plt.grid(True) # Add a grid for better readability.
        # Apply the consistent x-axis limits for saturation across all plots.
        plt.xlim(x_axis_min_limit, x_axis_max_limit)
        plt.tick_params(axis='both', which='major', labelsize=12) # Increased tick label font size


    # --- Step 1: Linear Regression for Each Unique Slope ---
    # Iterate through each unique slope value to perform linear regression.
    for slope_val in unique_slopes:
        # Get the consistent color for the current slope value.
        current_color = color_map[slope_val]

        # Filter the DataFrame to get data points only for the current slope.
        # Sort by 'sat' for proper plotting and regression.
        slope_data = data_df[data_df['slope'] == slope_val].sort_values(by='sat')
        sats_for_slope = slope_data['sat'].values
        k2s_for_slope = slope_data['k2'].values

        # Check if there are enough data points for linear regression (minimum 2 points).
        if len(sats_for_slope) < 2:
            print(f"Warning: Not enough unique 'sat' points ({len(sats_for_slope)}) for slope {slope_val} to perform linear regression. Skipping this slope.")
            continue # Skip to the next slope if insufficient data.

        try:
            # Perform linear regression: k2 as a function of sat for the current slope.
            # linregress returns slope, intercept, r-value, p-value, and standard error.
            slope_reg, intercept_reg, r_value, p_value, std_err = linregress(sats_for_slope, k2s_for_slope)

            # Predict the k2 value at the `target_sat` using the derived regression line.
            k2_at_target_sat = slope_reg * target_sat + intercept_reg
            # Store the slope and the predicted k2 value for the next interpolation step.
            interpolated_k2_for_slopes.append((slope_val, k2_at_target_sat))

            if plot_results:
                # Plot the original data points for the current slope.
                plt.plot(sats_for_slope, k2s_for_slope, 'o', color=current_color, markersize=8, label=f'Slope {slope_val} Data') # Increased markersize

                # Plot the linear regression line.
                # The line range is set to the consistent x-axis limits to ensure it spans the entire plot.
                sat_reg_range = np.linspace(x_axis_min_limit, x_axis_max_limit, 100)
                plt.plot(sat_reg_range, slope_reg * sat_reg_range + intercept_reg, '-', color=current_color, linewidth=2.5, label=f'Slope {slope_val} Regression')

                # Plot a red 'X' marker at the predicted k2 for the target_sat on this regression line.
                # The label is now an empty string to remove it from the legend.
                plt.plot(target_sat, k2_at_target_sat, 'X', markersize=12, color='red', mew=2, label='') # Removed label
                # Add the k2 value as text annotation near the red cross, slightly above it.
                plt.text(target_sat, k2_at_target_sat + 0.1, f'{k2_at_target_sat:.2f}', color='red', ha='center', va='bottom', fontsize=10) # Increased font size

        except ValueError as e:
            # Catch potential errors during linear regression (e.g., if data is all identical).
            print(f"Error during linear regression for slope {slope_val} at target_sat {target_sat}: {e}")
            continue # Continue to the next slope.

    if plot_results:
        plt.legend(fontsize=12) # Increased legend font size


    # --- Validation for Second Interpolation Step ---
    # Ensure that intermediate interpolated values were generated.
    if not interpolated_k2_for_slopes:
        raise ValueError("No valid interpolated k2 values could be generated across slopes for the given target_sat. Check data range and density.")

    # Prepare data for the second interpolation step (PCHIP).
    # Extract slopes and their corresponding predicted k2 values.
    slopes_for_final_interp = np.array([item[0] for item in interpolated_k2_for_slopes])
    k2_values_for_final_interp = np.array([item[1] for item in interpolated_k2_for_slopes])

    # Check if there are enough points for PCHIP interpolation (minimum 2 unique points).
    if len(slopes_for_final_interp) < 2:
        raise ValueError(f"Not enough unique slope points ({len(slopes_for_final_interp)}) after first interpolation step to perform PCHIP interpolation for target_slope.")

    # --- Step 2: PCHIP Interpolation for Target Slope ---
    # PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) is used.
    # IMPORTANT: PCHIP does NOT support extrapolation. The target_slope must be within
    # the range of the `slopes_for_final_interp` values.
    try:
        # Create the PCHIP interpolator object.
        pchip_interp_slope = PchipInterpolator(slopes_for_final_interp, k2_values_for_final_interp)
        # Perform the final interpolation to get k2 for the `target_slope`.
        final_k2 = pchip_interp_slope(target_slope)

        if plot_results:
            # --- Subplot 2: k2 vs Slope (PCHIP Interpolation) ---
            plt.subplot(1, 2, 2) # Create the second subplot.
            plt.title('Step 2: k2 vs Slope (PCHIP Interpolation)', fontsize=16) # Increased font size
            plt.xlabel('Slope', fontsize=14) # Increased font size
            plt.ylabel('k2 (predicted at target sat)', fontsize=14) # Increased font size
            plt.grid(True) # Add a grid.
            plt.tick_params(axis='both', which='major', labelsize=12) # Increased tick label font size
            # Plot the intermediate points (slope, k2_predicted_at_target_sat).
            plt.plot(slopes_for_final_interp, k2_values_for_final_interp, 'o', color='blue', markersize=8, label='Intermediate Points') # Increased markersize
            
            # Generate a finer range of slope values for plotting the smooth PCHIP curve.
            slope_interp_range = np.linspace(slopes_for_final_interp.min(), slopes_for_final_interp.max(), 100)
            # Plot the PCHIP interpolation curve.
            plt.plot(slope_interp_range, pchip_interp_slope(slope_interp_range), '-', color='green', linewidth=2.5, label='PCHIP Interpolation')
            # Plot a red 'X' marker at the final interpolated point.
            plt.plot(target_slope, final_k2, 'X', markersize=12, color='red', mew=2, label='Final Interpolated Point') # Increased markersize
            # Add the final k2 value as text annotation near the red cross, slightly above it.
            plt.text(target_slope, final_k2 + 0.1, f'{final_k2:.2f}', color='red', ha='center', va='bottom', fontsize=10) # Increased font size
            plt.legend(fontsize=12) # Increased legend font size

            # --- Update Super Title with Final k2 Value ---
            # Set the main title for the entire figure, including the calculated final k2.
            # This is done here because final_k2 is known at this point.
            fig.suptitle(f'Interpolation Steps for Target Sat={target_sat:.4f}, Target Slope={target_slope:.2f}, Final k2={final_k2:.4f}', fontsize=18) # Increased font size
            
            plt.tight_layout() # Adjust subplot parameters for a tight layout.
            # Return the final k2 value and the figure object. The figure object will be saved to PDF.
            return final_k2, fig
        else:
            # If plotting is not requested, just return the final k2 value.
            return final_k2

    except ValueError as e:
        # Catch errors if target_slope is out of range for PCHIP or other PCHIP issues.
        raise ValueError(f"Error during final PCHIP interpolation for target_slope {target_slope}. "
                         f"Ensure target_slope is within the range [{slopes_for_final_interp.min():.2f}, {slopes_for_final_interp.max():.2f}]: {e}")

    # This line should ideally not be reached if all paths return or raise an exception.
    return final_k2

# --- Example Usage and Data Setup ---

# 1. Prepare your data from the image into a pandas DataFrame.
#    This data represents empirical relationships between saturation (sat),
#    foreshore slope, and the k2 parameter of a Composite Weibull distribution.
#    The data is manually transcribed from Figure 5 of the Groenendijk & Van Gent (1998) report.
data = {
    'sat': [
        0.076, 0.095, 0.117, 0.119, 0.121,
        0.130, 0.134, 0.140, 0.146, 0.153,
        0.160, 0.167, 0.180,
        0.045, 0.050, 0.067, 0.072, 0.093,
        0.114, 0.114, 0.123, 0.138, 0.145,
        0.150,
        0.038, 0.039, 0.042, 0.065, 0.072,
        0.085, 0.085, 0.110, 0.116, 0.114,
        0.133, 0.144, 0.145, 0.158,
        0.030, 0.030, 0.041, 0.020, 0.062,
        0.070, 0.092, 0.106, 0.113, 0.130,
        0.152
    ],
    'slope': [
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
        100, 100, 100, 100,
        250, 250, 250, 250, 250, 250, 250, 250, 250, 250,
        250
    ],
    'k2': [
        4.733, 2.498, 3.979, 3.194, 3.684,
        2.970, 3.188, 3.286, 2.987, 3.178,
        2.380, 2.976, 3.190,
        2.473, 2.544, 2.389, 2.583, 3.380,
        4.030, 3.174, 3.561, 3.284, 3.574,
        3.229,
        1.988, 3.189, 6.977, 3.883, 4.148,
        4.411, 3.481, 3.861, 4.225, 3.467,
        4.363, 4.211, 3.305, 3.470,
        7.471, 5.319, 2.084, 3.937, 3.937,
        5.250, 3.038, 3.251, 4.128, 4.193,
        3.742
    ]
}
df = pd.DataFrame(data)

# --- Generate Multiple Cases for PDF Output ---

print("--- Generating Multiple Interpolation Examples (Linear Regression + PCHIP with Plots) ---")

# Define 3 different saturation values to be used for interpolation.
sat_values = [0.06, 0.12, 0.16]

# Define 4 specific slope values for interpolation, as requested.
slope_values = [10, 50, 100, 300]

# Create a list to store all combinations of saturation and slope values.
# This will result in 3 saturation values * 4 slope values = 12 total cases.
examples = []
for sat in sat_values:
    for slope in slope_values:
        examples.append((sat, slope))

# Print the generated examples for user reference.
print(f"Generated {len(examples)} examples:")
for i, (sat, slope) in enumerate(examples):
    print(f"  Example {i+1}: target_sat={sat}, target_slope={slope}")

# --- PDF Generation Loop ---
# Create a PDF file to save all charts. The file will be named 'k2_interpolate.pdf'.
# The 'with' statement ensures the PDF file is properly closed after all operations.
with PdfPages('k2_interpolate.pdf') as pdf:
    # Iterate through each generated example.
    for i, (target_sat, target_slope) in enumerate(examples):
        print(f"\n--- Running Example {i+1} ---")
        try:
            # Call the 'f' function with plot_results=True.
            # The function will return the calculated k2 value and the Matplotlib figure object.
            result_k2, fig = f(target_sat, target_slope, df, plot_results=True)
            print(f"Result for sat={target_sat} and slope={target_slope}: k2 = {result_k2:.4f}")
            
            # Save the current figure to the PDF.
            # 'bbox_inches='tight'' ensures that the bounding box of the figure is tight
            # around the plotted content, helping the chart fit the page width when opened.
            pdf.savefig(fig, bbox_inches='tight')
            # Close the figure to free up memory, which is important when generating many plots.
            plt.close(fig)

        except ValueError as e:
            # Handle specific ValueError exceptions that might occur during interpolation.
            print(f"Example {i+1}: Error calculating for sat={target_sat} and slope={target_slope}: {e}")
        except KeyError as e:
            # Handle KeyError if DataFrame columns are missing.
            print(f"Example {i+1}: Error: {e}")

# Inform the user that the PDF generation is complete.
print("\nAll charts have been saved to 'k2_interpolate.pdf'.")
