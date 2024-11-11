import numpy as np
import argparse
from scipy.optimize import curve_fit
from given_functions import gaussian_initial_estimates
from datetime import datetime
from fitting_functions import combined_model_gaussian, combined_model_double_gaussian, resolution_function, linear_model
from plotting_functions import plot_spectrum, plot_fit_result, plot_calibration_curve, plot_detector_efficiencies, plot_resolution
from given_data import BGO_spectra, source_data, NaI_spectra, CdTe_spectra, detector_geometry
from data_loading_functions import load_background, load_and_extract_roi

def parse_arguments():
    """ parse the arguments """

    parser = argparse.ArgumentParser()
    parser.add_argument('detector', type=str, help='Input name of the detector (bgo, nai, cdte)')
    return parser.parse_args()

# Fit a peak within a region of interest (ROI)
def fit_peak(roi_channels, roi_net_count_rates, model):
    """
    Fits a Gaussian or double Gaussian model to a spectral peak in the ROI.
    Returns: Fitted parameters, the fit function, and formatted parameter text.
    """

    background_guess = (0, 0, min(roi_net_count_rates))
    peak_guess = gaussian_initial_estimates(roi_channels, roi_net_count_rates)
    
    if model == 'gaussian':
        initial_guess = background_guess + peak_guess
        popt, covariance_matrix = curve_fit(combined_model_gaussian, roi_channels, roi_net_count_rates, p0=initial_guess, maxfev=5000)
        fit_func = combined_model_gaussian
        param_names = ["a", "b", "c", "mu", "sigma", "amplitude"]
    elif model == 'double_gaussian':
        initial_guess = background_guess + 2 * peak_guess
        popt, covariance_matrix = curve_fit(combined_model_double_gaussian, roi_channels, roi_net_count_rates, p0=initial_guess, maxfev=5000)
        fit_func = combined_model_double_gaussian
        param_names = ["a", "b", "c", "mu1", "sigma1", "amp1", "mu2", "sigma2", "amp2"]
    
    perr = np.sqrt(np.diag(covariance_matrix))
    
    # Prepare text summary of parameters with uncertainties
    param_text = f"Fit parameters for ({model.capitalize()}):\n"
    for name, param, error in zip(param_names, popt, perr):
        param_text += f"{name}: {param:.2f} ± {error:.2f}\n"

    return popt, fit_func, param_text


# Calibrate detector based on known isotope data
def calibrate_detector(isotope_data):
    """
    Performs linear calibration based on spectral centroids and known energies.
    Returns: Parameters of linear fit and their uncertainties.
    """
    
    all_centroids = np.array([centroid for entry in isotope_data for centroid in entry["centroids"]])
    all_energies = np.array([energy for entry in isotope_data for energy in entry["energies"]])
    
    # Using curve_fit for linear calibration
    popt, pcov = curve_fit(linear_model, all_centroids, all_energies)
    slope, intercept = popt
    slope_err, intercept_err = np.sqrt(np.diag(pcov))
    
    print("\nFitted calibration parameters:")
    print(f"Slope: {slope:.2f} ± {slope_err:.2f}")
    print(f"Intercept: {intercept:.2f} ± {intercept_err:.2f}")
    
    # Optional: Check for non-linearity using a quadratic fit
    quadratic_fit = np.polyfit(all_centroids, all_energies, 2)
    quadratic_coefficient = quadratic_fit[0]
    
    if abs(quadratic_coefficient) < 0.001:  # Adjustable tolerance 
        print("Quadratic coefficient is close to 0; a linear fit is likely sufficient.\n")
    else:
        print("Quadratic coefficient is significant; consider using a quadratic calibration.\n")
    
    return slope, intercept

# Calculate absolute and intrinsic detector efficiencies
def calculate_detector_efficiency(spectra, isotope_data, source_data, CdTe, detector_geometry):
    """
    Calculates the absolute and intrinsic efficiencies of the detector for each peak.
    Returns: Efficiencies, Intrinsic Efficiencies and Energies used for calculation.
    """
    
    # Detector configuration
    detector_radius = detector_geometry["detector_radius"]
    source_distance = detector_geometry["source_distance"]

    # Adjust area based on detector type
    if CdTe:
        detector_area = 0.25  #CdTe Area
    else:
        detector_area = np.pi * (detector_radius ** 2)   #BGO/NaI Area
    
    geometry_factor = detector_area / (4 * np.pi * (source_distance ** 2))

    efficiencies, intrinsic_efficiencies, energies_list = [], [], []
    
    print("Isotopes and Energies:\n")

    for spectrum, roi_counts in zip(spectra, [entry["amplitudes"] for entry in isotope_data]):
        _, isotope, _, _, energies = spectrum

        # Match source data for decay and emission details
        source_info = next((source for source in source_data if source['isotope'] == isotope), None)
        if not source_info:
            continue

        # Calculate current source activity with decay correction
        calibration_date = source_info['calibration_date']
        initial_activity = source_info['initial_activity'] * 37000  # Convert to Bq 
        half_life = source_info['half_life']
        
        elapsed_time_days = (datetime.now() - calibration_date).days
        decay_constant = np.log(2) / half_life
        current_activity = initial_activity * np.exp(-decay_constant * elapsed_time_days)

        # Calculate efficiencies for each peak energy
        for energy, roi_count in zip(energies, roi_counts):
            emission_fraction = source_info['emission_fractions'].get(energy, 1.0)
            absolute_efficiency = roi_count / (current_activity * emission_fraction)
            intrinsic_efficiency = absolute_efficiency / geometry_factor

            efficiencies.append(absolute_efficiency)
            intrinsic_efficiencies.append(intrinsic_efficiency)
            energies_list.append(energy)
            print(f"Isotope: {isotope}, Energy: {energy} keV, Absolute Efficiency: {absolute_efficiency:.4e}, Intrinsic Efficiency: {intrinsic_efficiency:.4e}") #printing for understanding if needed
    
    print(f"\nGeometry Factor (G): {geometry_factor:.4e}\n")
    return efficiencies, intrinsic_efficiencies, energies_list

# Calculate the resolution of the detector
def calculate_resolution(isotope_data, slope, intercept):
    """
    Computes detector resolution as a percentage of the Full Width at Half Maximum (FWHM).
    Returns: Resolutions, Energies and the fitted parameters with their uncertainties.
    """
    # Skip Barium data points
    all_centroids = [centroid for entry in isotope_data if entry["isotope"] != "133-Ba" for centroid in entry["centroids"]]
    all_sigmas = [sigma for entry in isotope_data if entry["isotope"] != "133-Ba" for sigma in entry["sigmas"]]

    # Calculate FWHM and resolutions
    fwhm = 2.355 * np.array(all_sigmas)
    centroid_energies = slope * np.array(all_centroids) + intercept
    resolutions = (fwhm / centroid_energies) * 100
    resolutions_squared = resolutions ** 2
    
    # Fit the resolution data to the resolution function model
    res_popt, res_pcov = curve_fit(resolution_function, centroid_energies, resolutions_squared, p0=[1, 1, 1], bounds=(0, np.inf), maxfev=5000)
    
    # Extract parameters and uncertainties
    a, b, c = res_popt
    a_err, b_err, c_err = np.sqrt(np.diag(res_pcov))

    print("\nFitted resolution function parameters:")
    print(f"a: {a:.2f} ± {a_err:.2f}")
    print(f"b: {b:.2f} ± {b_err:.2f}")
    print(f"c: {c:.2f} ± {c_err:.2f}")

    return resolutions, centroid_energies, res_popt

# Process a single spectrum, including ROI extraction and background subtraction, plotting and fitting
def process_spectrum(spectrum, count_rates_bg, isotope_data, CdTe):
    """
    Processes spectrum data, fits ROIs, plots raw and fitted spectra 
    and appends results to isotope_data.
    """
    
    file_path, isotope, rois, measurement_time, energies = spectrum
    if isotope == "Background":
        return

    channels, count_rates, count_rates_bgsubtracted, roi_data = load_and_extract_roi(file_path, rois, measurement_time, count_rates_bg, CdTe)
    plot_spectrum(channels, count_rates, count_rates_bgsubtracted, roi_data, isotope, measurement_time)
    
    total_roi_count_rates, centroids, sigmas, amplitudes = [], [], [], []

    for model, roi_channels, roi_net_count_rates in roi_data:
        popt, fit_func, perr = fit_peak(roi_channels, roi_net_count_rates, model)
        plot_fit_result(channels, count_rates_bgsubtracted, roi_channels, roi_net_count_rates, popt, perr, fit_func, isotope, measurement_time, model)
        
        if model == 'gaussian':
            centroids.append(popt[3])
            sigmas.append(popt[4])
            amplitudes.append(popt[5])
        elif model == 'double_gaussian':
            centroids.extend([popt[3], popt[6]])
            sigmas.extend([popt[4], popt[7]])
            amplitudes.extend([popt[5], popt[8]])

    isotope_data.append({
        "isotope": isotope,
        "total_roi_count_rates": total_roi_count_rates,
        "centroids": centroids,
        "sigmas": sigmas,
        "energies": energies,
        "amplitudes": amplitudes
    })

# Main function for entire analysis process
def main(detector):
    """
    Main execution function for spectrum processing and analysis.
    """
    if detector=="bgo":
        spectra = BGO_spectra #input detector name
        CdTe = False
    elif detector=="nai":
        spectra = NaI_spectra #input detector name
        CdTe = False
    elif detector=="cdte":
        spectra = CdTe_spectra #input detector name
        CdTe = True
    else:
        raise ValueError("Please insert correct detector name")


    count_rates_bg = load_background(spectra, CdTe)
    isotope_data = []

    #iterates process function through all spectra
    for spectrum in spectra:
        process_spectrum(spectrum, count_rates_bg, isotope_data, CdTe)
    

    slope, intercept = calibrate_detector(isotope_data)
    plot_calibration_curve(isotope_data, slope, intercept)
    
    efficiencies, intrinsic_efficiencies, energies_list = calculate_detector_efficiency(spectra, isotope_data, source_data, CdTe, detector_geometry)
    plot_detector_efficiencies(efficiencies, intrinsic_efficiencies, energies_list)

    resolutions, centroid_energies, res_popt = calculate_resolution(isotope_data, slope, intercept)
    plot_resolution(centroid_energies, resolutions, res_popt)

if __name__ == '__main__':
    args = parse_arguments()
    main(args.detector)