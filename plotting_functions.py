import matplotlib.pyplot as plt
import numpy as np
from fitting_functions import resolution_function

def plot_spectrum(channels, count_rates, count_rates_bgsubtracted, roi_data, isotope, measurement_time):
    """
    Plots raw and background-subtracted spectra with highlighted regions of interest (ROIs).
    """

    plt.figure(figsize=(10, 6))
    plt.plot(channels, count_rates, label="Raw Spectrum", color="gray")
    plt.plot(channels, count_rates_bgsubtracted, label="Background-Subtracted Spectrum", color="blue")

    # Highlight each ROI
    for model, roi_channels, roi_net_count_rates in roi_data:
        plt.plot(roi_channels, roi_net_count_rates, label=f"ROI ({model})", linewidth=2, color="orange")

    plt.xlabel("Channel")
    plt.ylabel("Count Rate (counts/sec)")
    plt.title(f"Spectrum Data for {isotope} ({measurement_time} s)")
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.show()

def plot_fit_result(channels, counts, roi_channels, roi_net_count_rates, popt, param_text, fit_func, isotope, measurement_time, model):
    """
    Plots the fitted model over the ROI with raw spectrum data for comparison.
    """

    plt.figure(figsize=(10, 6))
    plt.scatter(channels, counts, label="Raw Spectrum (without background)", marker='+', color="gray")
    plt.scatter(roi_channels, roi_net_count_rates, marker='+', color="green", label="Region of Interest")
    
    # Generate fit data
    fit_x = np.linspace(min(roi_channels), max(roi_channels), 500)
    fit_y = fit_func(fit_x, *popt)
    plt.plot(fit_x, fit_y, color="red", label="Fitted Model")
    
    # Display fit parameters within the plot area
    plt.text(0.95, 0.6, param_text, transform=plt.gca().transAxes, fontsize=6,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="lightyellow"))

    plt.xlabel("Channel")
    plt.ylabel("Count Rate (counts/sec)")
    plt.legend()
    plt.title(f"{isotope} {model.capitalize()} Spectrum Fit ({measurement_time} s)")
    plt.grid(True)
    plt.show()

def plot_calibration_curve(isotope_data, slope, intercept):
    """
    Plots calibration curve using centroids and known energies with a fitted linear model.
    """

    all_centroids = [centroid for entry in isotope_data for centroid in entry["centroids"]]
    all_energies = [energy for entry in isotope_data for energy in entry["energies"]]

    fit_x = np.linspace(min(all_centroids), max(all_centroids), 100)
    fit_y = slope * fit_x + intercept
    plt.figure(figsize=(8, 6))
    plt.scatter(all_centroids, all_energies, color="blue", label="Calculated Calibrations")
    plt.plot(fit_x, fit_y, color="pink", linestyle="--", label="Fitted Model")
    plt.xlabel("Channel (centroid)")
    plt.ylabel("Energy (keV)")
    plt.title("Calibration Curve")
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_detector_efficiencies(efficiencies, intrinsic_efficiencies, energies):
    """
    Plots detector absolute and intrinsic efficiency against energy.
    """

    # Absolute efficiency plot
    plt.figure()
    plt.scatter(energies, efficiencies, label='Absolute Efficiency')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Absolute Efficiency')
    plt.title('Detector Absolute Efficiency vs. Energy')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Intrinsic efficiency plot
    plt.figure()
    plt.scatter(energies, intrinsic_efficiencies, label='Intrinsic Efficiency')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Intrinsic Efficiency')
    plt.title('Detector Intrinsic Efficiency vs. Energy')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Fit a quadratic model on a log-log scale
    log_energies = np.log(energies)
    log_int_efficiencies = np.log(intrinsic_efficiencies)
    quadratic_fit, cov_matrix = np.polyfit(log_energies, log_int_efficiencies, 2, cov=True)

    # Extract fit parameters and their uncertainties
    a, b, c = quadratic_fit
    a_err, b_err, c_err = np.sqrt(np.diag(cov_matrix))

    print("Fitted Intrinsic Efficiency parameters (log-log scale):")
    print(f"a: {a:.2f} ± {a_err:.2f}")
    print(f"b: {b:.2f} ± {b_err:.2f}")
    print(f"c: {c:.2f} ± {c_err:.2f}")

    # Generate fit curve
    fit_x = np.linspace(min(log_energies), max(log_energies), 500)
    fit_y = np.polyval(quadratic_fit, fit_x)

    # Plot the log-log data and the fitted model
    plt.figure(figsize=(9, 5))
    plt.scatter(energies, intrinsic_efficiencies, label="Calculated Efficiencies")
    plt.plot(np.exp(fit_x), np.exp(fit_y), color="pink", linestyle="--", label="Fitted model")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Intrinsic Efficiency")
    plt.title("Intrinsic Efficiency and Fitted Model (Log-Log Scale)")
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_resolution(centroid_energies, resolutions, res_popt):
    """
    Plots calculated resolutions with a fitted resolution model.
    """

    plt.figure(figsize=(8, 6))
    plt.scatter(centroid_energies, resolutions, label="Calculated Resolutions", color="blue")
    
    # Fit resolution curve
    fit_energies = np.linspace(min(centroid_energies), max(centroid_energies), 500)
    fit_resolutions = np.sqrt(resolution_function(fit_energies, *res_popt))
    plt.plot(fit_energies, fit_resolutions, label="Fitted Resolution Model", color="pink", linestyle="--")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Resolution (%)")
    plt.title("Detector Energy Resolution")
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_detector_efficiencies_off_axis(efficiencies, isotope_data):
    """
    Plots absolute efficiency of the detector as a function of angle.
    """

    all_degrees = [entry["degrees"] for entry in isotope_data]
    
    plt.figure()
    plt.plot(all_degrees, efficiencies, label='Absolut Efficiency', marker="o")
    plt.xlabel('Angle [degrees]')
    plt.ylabel('Absolut Efficiency')
    plt.title('Detector Absolut Efficiency at Different Angles')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_resolution_off_axis(isotope_data, resolutions):
    """
    Plots energy resolution of the detector as a function of angle.
    """

    all_degrees = [entry["degrees"] for entry in isotope_data]

    plt.figure(figsize=(8, 6))
    plt.plot(all_degrees, resolutions, label="Calculated Resolutions", marker="o")
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Resolution (%)")
    plt.title("Detector Energy Resolution at Different Angles")
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_countrate_off_axis(isotope_data):
    """
    Plots total count rate as a function of angle.
    """

    all_degrees = [entry["degrees"] for entry in isotope_data]
    total_roi_count_rates = [entry["total_roi_count_rates"] for entry in isotope_data]

    plt.figure(figsize=(8, 6))
    plt.plot(all_degrees, total_roi_count_rates, label="Count Rate", marker="o")
    plt.xlabel("Angle [degrees]")
    plt.ylabel("Count Rate (count/sec)")
    plt.title("Total Count Rate of the Peak at Different Angles")
    plt.legend()
    plt.grid(True)
    plt.show()
