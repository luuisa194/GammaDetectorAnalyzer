import numpy as np
import argparse
from scipy.optimize import curve_fit
from given_functions import gaussian_initial_estimates
from datetime import datetime
from given_data import BGO_off_axis_spectra, source_data, NaI_off_axis_spectra, CdTe_off_axis_spectra
from fitting_functions import combined_model_gaussian, linear_model
from plotting_functions import plot_detector_efficiencies_off_axis, plot_resolution_off_axis, plot_countrate_off_axis
from data_loading_functions import load_background, load_and_extract_roi

def parse_arguments():
    """ parse the arguments """

    parser = argparse.ArgumentParser()
    parser.add_argument('detector', type=str, help='Input name of the detector (bgo, nai, cdte)')
    return parser.parse_args()



# Data processing functions
def fit_peak(roi_channels, roi_net_count_rates, model):
    """
    Fits a Gaussian or double Gaussian model to a spectral peak in the ROI.
    Returns: Fitted parameters, the fit function, and formatted parameter text.
    """

    background_guess = (0, 0, min(roi_net_count_rates))
    peak_guess = gaussian_initial_estimates(roi_channels, roi_net_count_rates)
    initial_guess = background_guess + peak_guess
    
    # Fit Gaussian model to the ROI data
    popt, covariance_matrix = curve_fit(combined_model_gaussian, roi_channels, roi_net_count_rates, p0=initial_guess, maxfev=5000)
    fit_func = combined_model_gaussian
    perr = np.sqrt(np.diag(covariance_matrix))
    return popt, fit_func, perr

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

def calculate_detector_efficiency(spectra, isotope_data, Am_data):
    """
    Calculates the absolute efficiency of the detector for each angle.
    Returns: Absolute efficiencies calculated for each energy.
    """
    all_amplitudes = [entry["amplitudes"] for entry in isotope_data]  
    efficiencies = []

    for spectrum, roi_counts in zip(spectra, all_amplitudes):
        _, angle, _, _, energies = spectrum

        # Calculate current activity accounting for decay since calibration date
        calibration_date = Am_data['calibration_date']
        initial_activity = Am_data['initial_activity'] * 37000  # Convert to Bq
        half_life = Am_data['half_life']
        
        elapsed_time_days = (datetime.now() - calibration_date).days
        decay_constant = np.log(2) / half_life
        current_activity = initial_activity * np.exp(-decay_constant * elapsed_time_days)

        # Calculate efficiency for each peak energy
        for energy, roi_count in zip(energies, roi_counts):
            emission_fraction = Am_data['emission_fractions'].get(energy, 1.0)
            absolute_efficiency = roi_count / (current_activity * emission_fraction)
            efficiencies.append(absolute_efficiency)
            print(f"Angle: {angle} degrees, Absolut Efficiency: {absolute_efficiency:.4e}")
        
    return efficiencies

def calculate_resolution(isotope_data, slope, intercept):
    """
    Calculates the resolution of the detector based on FWHM.
    Returns: Ccalculated resolutions.
    """
    all_centroids = [entry["centroids"] for entry in isotope_data]
    all_sigmas = [entry["sigmas"] for entry in isotope_data]

    fwhm = 2.355 * np.array(all_sigmas)
    centroid_energies = slope * np.array(all_centroids) + intercept
    resolutions = (fwhm / centroid_energies) * 100
    
    return resolutions

def process_spectrum(spectrum, count_rates_bg, isotope_data,CdTe):
    """
    Processes a single spectrum, extracting and fitting ROIs and appends results to isotope_data.
    """
    file_path, degrees, rois, measurement_time, energies = spectrum
    if degrees == "Background":
        return

    # Load and extract ROI data with background subtraction
    _, _, _, roi_data = load_and_extract_roi(file_path, rois, measurement_time, count_rates_bg,CdTe)
    
    total_roi_count_rates, centroids, sigmas, amplitudes = [], [], [], []

    for model, roi_channels, roi_net_count_rates in roi_data:
        popt, _, _ = fit_peak(roi_channels, roi_net_count_rates, model)
        
        # Sum ROI count rates for efficiency calculations
        total_roi_count_rate = np.sum(roi_net_count_rates)
        total_roi_count_rates.append(total_roi_count_rate)

        # Extract centroid (mean) and sigma
        centroids.append(popt[3])
        sigmas.append(popt[4])
        amplitudes.append(popt[5])

    # Append processed data to isotope_data
    isotope_data.append({
        "degrees": degrees,
        "total_roi_count_rates": total_roi_count_rates,
        "centroids": centroids,
        "sigmas": sigmas,
        "energies": energies,
        "amplitudes": amplitudes
    })

def main(detector):
    """
    Main function for the analysis of off-axis spectra.
    """
    if detector=="bgo":
        spectra = BGO_off_axis_spectra #input detector name
        CdTe = False
    elif detector=="nai":
        spectra = NaI_off_axis_spectra #input detector name
        CdTe = False
    elif detector=="cdte":
        spectra = CdTe_off_axis_spectra #input detector name
        CdTe = True
    else:
        raise ValueError("Please insert correct detector name")

    # Retrieve 241-Am source data
    Am_data = next(item for item in source_data if item["isotope"] == "241-Am")

    count_rates_bg = load_background(spectra, CdTe)
    isotope_data = []

    # Process each spectrum in the dataset
    for spectrum in spectra:
        process_spectrum(spectrum, count_rates_bg, isotope_data,CdTe)

    # Perform calibration and calculate efficiencies
    slope, intercept = calibrate_detector(isotope_data)
    efficiencies = calculate_detector_efficiency(spectra, isotope_data, Am_data)
    plot_detector_efficiencies_off_axis(efficiencies, isotope_data)

    # Calculate and plot resolutions and count rates
    resolutions = calculate_resolution(isotope_data, slope, intercept)
    plot_resolution_off_axis(isotope_data, resolutions)
    plot_countrate_off_axis(isotope_data)

if __name__ == '__main__':
    args = parse_arguments()
    main(args.detector)