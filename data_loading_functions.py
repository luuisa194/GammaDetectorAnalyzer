from given_functions import filter_in_interval

def load_spe_file(file_path, measurement_time):
    """
    Loads data from an SPE file, extracting count rates per channel.
    Returns: Raw channels and count rates.
    """
    counts = []
    with open(file_path, 'r') as file:
        reading_data = False
        for line in file:
            # Start reading data after encountering "$DATA:"
            if "$DATA:" in line:
                reading_data = True
                continue
            if reading_data:
                # Collect counts while they are numeric, stop on encountering "$" again
                if line.strip().isdigit():
                    counts.append(int(line.strip()))
                elif line.startswith("$"):
                    break
    # Convert counts to count rates by dividing by the measurement time
    count_rates = [count / measurement_time for count in counts]
    channels = list(range(len(count_rates))) 
    return channels, count_rates

def load_mca_file(file_path, measurement_time):
    """
    Loads data from an MCA file, extracting count rates per channel.
    Returns: Raw channels and count rates.
    """
    counts = []
    with open(file_path, 'r') as file:
        reading_data = False
        for line in file:
            # Start reading data after encountering "<<DATA>>"
            if "<<DATA>>" in line:
                reading_data = True
                continue
            if reading_data:
                # Collect counts until non-numeric line or line starting with "$" is encountered
                if line.strip().isdigit():
                    counts.append(int(line.strip()))
                elif line.startswith("<<"):
                    break
    count_rates = [count / measurement_time for count in counts]
    channels = list(range(len(count_rates))) 
    return channels, count_rates

def load_background(spectra, CdTe):
    """
    Loads background data from a spectrum designated as background.
    Returns: Background count rates per channel.
    """
    for spectrum in spectra:
        file_path, isotope, _, measurement_time, _ = spectrum
        if isotope == "Background":
            # Choose appropriate file loader based on detector type
            if CdTe:
                _, bg_count_rates = load_mca_file(file_path, measurement_time)
            else:
                _, bg_count_rates = load_spe_file(file_path, measurement_time)
            return bg_count_rates
    return None

def load_and_extract_roi(file_path, rois, measurement_time, count_rates_bg, CdTe):
    """
    Loads spectrum data and extracts regions of interest (ROIs) with background subtraction.
    Returns: Channels, count rates, background-subtracted count rates, and ROI information.
    """
    # Load channels and count rates using the appropriate file loader
    if CdTe:
        channels, count_rates = load_mca_file(file_path, measurement_time)
    else:
        channels, count_rates = load_spe_file(file_path, measurement_time)
    
    # Perform background subtraction, ensuring no negative rates
    count_rates_bgsubtracted = [max(a - b, 0) for a, b in zip(count_rates, count_rates_bg)]
    
    roi_data = []
    for model, (roi_start, roi_end) in rois:
        # Extract channels and count rates within each ROI
        roi_channels, roi_net_count_rates = filter_in_interval(channels, count_rates_bgsubtracted, roi_start, roi_end)
        roi_data.append((model, roi_channels, roi_net_count_rates))
    
    return channels, count_rates, count_rates_bgsubtracted, roi_data
