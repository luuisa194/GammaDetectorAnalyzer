Gamma Spectroscopy Analysis 

Overview:
This repository includes tools and data for analyzing gamma-ray spectra using different detectors. The data consists of spectra measurements for each detector with different isotopes, as well as spectra for Americium at various angles. The main scripts to run are `GammaDetAnalyzer.py` and `OffAxisAnalyzer.py`.

Depending on the analysis to perform, the scripts can be executed from a terminal or command prompt. First, navigate to the directory containing all required files, including the scripts and data folders. When running the scripts, specify the detector type (options: “bgo”, “nai”, or “cdte”) as an argument. For example, use the command `python GammaDetAnalyzer.py bgo` to analyze the BGO detector data.

Purpose and Usage of Scripts:
- `GammaDetAnalyzer.py`: Analyzes spectra for different isotopes across the whole detector. This script generates output plots for each raw spectrum with highlighted ROIs, fits peaks, and creates plots for calibration, efficiency (absolute and intrinsic), and resolution. 
- `OffAxisAnalyzer.py`: Specifically analyzes the off-axis response of Americium for each detector. This script produces plots of absolute efficiency, resolution, and count rate as functions of the angle.

Both scripts use data for calibration, efficiency, and resolution, iterating through spectra to perform peak fitting for each one. Results are saved in dictionaries that support subsequent calibration, efficiency, and resolution calculations.

Supporting Modules:
The analysis relies on several supporting modules that provide essential functionality:
- `plotting_functions`: Generates plots for raw spectra, fitted spectra with parameter uncertainties, and efficiency and resolution.
- `loading_functions`: Loads data from various file types, extracts relevant information, stores background data separately, and performs background subtraction. This module also handles the extraction and storage of Regions of Interest (ROIs).
- `given_functions`: Includes specific functions provided in the instructions.
- `given_data`: Contains key data for each detector, such as file paths, isotopes, fitting models, ROIs, measurement times, and energy values. This module also provides off-axis spectra for each detector with angles, fitting models, and calibration details, allowing for adjustments if different files or geometries are used.


Output Details:
GammaDetAnalyzer.py
- Outputs detailed plots for each spectrum, highlighting ROIs and showing fitted peaks with uncertainties.
- Produces a calibration plot, absolute and intrinsic efficiency plots over energy, efficiency comparison plots in log-log scale, and a resolution plot. Calibration, energies, and resolution also include fit parameters.
  
OffAxisAnalyzer.py
- Provides plots showing absolute efficiency, resolution, and count rate as functions of angle.
  
The console output for both scripts includes fitted calibration parameters, absolute and intrinsic efficiency values, and notes on any fitting issues encountered. 

Errors:

-fitting of peaks: adjust the ROIs , adjust (maxfev), or remove ROI if it cant be fit

-fitting of energies/resolution: may happen because of faulty data/peaks try increasing (maxfev)


Notes:
The `cdte` detector data shows limitations in fitting, particularly for Cobalt, due to poor data quality. If cobalt is included in the analysis, the program will plot the raw cobalt data but may not proceed with further outputs. The Resolution cannot be fit properly for this detector, as there are not enough data points.
