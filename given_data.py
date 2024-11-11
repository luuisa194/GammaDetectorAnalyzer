from datetime import datetime

#Spectrum data for different detectors and sources
#Each entry in the list contains: 
#(file_path, isotope, [(model_type, ROI_range)], measurement_time, energies)

BGO_spectra = [
    (("BGO/Americium_0degrees_200s_BGO_uncal.Spe"), "241-Am", [("gaussian", (5, 80))], 200, [59.5]),
    (("BGO/Barium_0degrees_300s_BGO_uncal.Spe"), "133-Ba", [("gaussian", (25, 60)), ("double_gaussian", (120, 220))], 300, [53.2, 356.0, 383.8]),
    (("BGO/Caesium_0degrees_200s_BGO_uncal.Spe"), "137-Cs", [("gaussian", (280, 400))], 200, [661.7]),
    (("BGO/Cobalt_0degrees_600s_BGO_uncal.Spe"), "60-Co", [("double_gaussian", (560, 770))], 600, [1173.2, 1332.5]),
    (("BGO/Background_0degrees_600s_BGO_uncal.Spe"), "Background", [("-", (0,0))], 600, [])
]

NaI_spectra = [
    (("NaI/Americium_0degrees_200s_NaI_uncal.Spe"), "241-Am", [("gaussian", (20, 55))], 200, [59.5]),
    (("NaI/Barium_0degrees_300s_NaI.Spe"), "133-Ba", [("gaussian", (25, 60)), ("double_gaussian", (100, 220))], 300, [53.2, 356.0, 383.8]),
    (("NaI/Caesium_0degrees_200s_NaI_uncal.Spe"), "137-Cs", [("gaussian", (240, 390))], 200, [661.7]),
    (("NaI/Cobalt_0degrees_600s_NaI_uncal.Spe"), "60-Co", [("double_gaussian", (450, 600))], 600, [1173.2, 1332.5]),
    (("NaI/Background_0degrees_600s_NaI.Spe"), "Background", [("-", (0,0))], 600, [])
]

CdTe_spectra = [
    (("CdTe/Americium_0degrees_200s_CdTe_uncal.mca"), "241-Am", [("gaussian", (300,380))], 200, [59.5]),
    (("CdTe/Barium_0degrees_300s_CdTe_uncal.mca"), "133-Ba", [("gaussian", (50, 200)), ("double_gaussian", (550, 740))], 300, [53.2, 356.0, 383.8]),
    (("CdTe/Caesium_0degrees_200s_CdTe_uncal.mca"), "137-Cs", [("gaussian", (550, 790))], 200, [661.7]),
    #(("CdTe/Cobalt_0degrees_600s_CdTe_uncal.mca"), "60-Co", [("double_gaussian", (450, 600))], 600, [1173.2, 1332.5]),
    #cobalt not readable
    (("CdTe/Background_0degrees_600s_CdTe.mca"), "Background", [("-", (0,0))], 600, [])
]

#Spectrum data for different angles at detectors
#Each entry in the list contains: 
#(file_path, angle, [(model_type, ROI_range)], measurement_time, energies)

BGO_off_axis_spectra =[
    (("BGO_off_axis/Americium_0degrees_200s_BGO_uncal.Spe"), 0, [("gaussian", (0, 100))], 200, [59.5]),
    (("BGO_off_axis/Americium_15degrees_200s_BGO_uncal.Spe"), 15, [("gaussian", (0, 100))], 200, [59.5]),
    (("BGO_off_axis/Americium_30degrees_200s_BGO_uncal.Spe"), 30, [("gaussian", (0, 100))], 200, [59.5]),
    (("BGO_off_axis/Americium_45degrees_200s_BGO_uncal.Spe"), 45, [("gaussian", (0, 100))], 200, [59.5]),
    (("BGO_off_axis/Americium_60degrees_200s_BGO_uncal.Spe"), 60, [("gaussian", (0, 100))], 200, [59.5]),
    (("BGO_off_axis/Americium_75degrees_200s_BGO_uncal.Spe"), 75, [("gaussian", (0, 100))], 200, [59.5]),
    (("BGO_off_axis/Americium_90degrees_200s_BGO_uncal.Spe"), 90, [("gaussian", (0, 100))], 200, [59.5]),
    (("BGO_off_axis/Background_0degrees_600s_BGO_uncal.Spe"), "Background", [("-", (0,0))], 600, [])
]

NaI_off_axis_spectra = [
    (("NaI_off_axis/Americium_0degrees_200s_NaI.Spe"), 0, [("gaussian", (20, 55))], 200, [59.5]),
    (("NaI_off_axis/Americium_15degrees_200s_NaI.Spe"), 15, [("gaussian", (20, 55))], 200, [59.5]),
    (("NaI_off_axis/Americium_30degrees_200s_NaI.Spe"), 30, [("gaussian", (20, 55))], 200, [59.5]),
    (("NaI_off_axis/Americium_45degrees_200s_NaI.Spe"), 45, [("gaussian", (20, 55))], 200, [59.5]),
    (("NaI_off_axis/Americium_60degrees_200s_NaI.Spe"), 60, [("gaussian", (20, 55))], 200, [59.5]),
    (("NaI_off_axis/Americium_75degrees_200s_NaI.Spe"), 75, [("gaussian", (20, 55))], 200, [59.5]),
    (("NaI_off_axis/Americium_90degrees_200s_NaI.Spe"), 90, [("gaussian", (20, 55))], 200, [59.5]),
    (("NaI_off_axis/Background_0degrees_600s_NaI.Spe"), "Background", [("-", (0,0))], 600, [])
]

CdTe_off_axis_spectra = [
    (("CdTe_off_axis/Americium_0degrees_200s_CdTe_uncal.mca"), 0, [("gaussian", (1150,1250))], 200, [59.5]),
    (("CdTe_off_axis/Americium_15degrees_200s_CdTe_uncal.mca"), 15, [("gaussian", (1150,1250))], 200, [59.5]),
    (("CdTe_off_axis/Americium_30degrees_200s_CdTe_uncal.mca"), 30, [("gaussian", (1150,1250))], 200, [59.5]),
    (("CdTe_off_axis/Americium_45degrees_200s_CdTe_uncal.mca"), 45, [("gaussian", (1150,1250))], 200, [59.5]),
    (("CdTe_off_axis/Americium_60degrees_200s_CdTe_uncal.mca"), 60, [("gaussian", (1150,1250))], 200, [59.5]),
    (("CdTe_off_axis/Americium_75degrees_200s_CdTe_uncal.mca"), 75, [("gaussian", (1150,1250))], 200, [59.5]),
    (("CdTe_off_axis/Americium_90degrees_200s_CdTe_uncal.mca"), 90, [("gaussian", (1150,1250))], 200, [59.5]),
    (("CdTe/Background_0degrees_600s_CdTe.mca"), "Background", [("-", (0,0))], 600, [])
]

#Source data for different isotopes, including calibration and decay information
source_data = [
    {
        "isotope": "241-Am",
        "calibration_date": datetime(1979,12,1),  
        "initial_activity": 5.0,  
        "half_life": 432 * 365.25,  
        "emission_fractions": {59.5: 0.36} 
    },
    {
        "isotope": "137-Cs",
        "calibration_date": datetime(1979,12,1),
        "initial_activity": 3.7,  
        "half_life": 30.09 * 365.25,  
        "emission_fractions": {661.6: 0.85} 
    },
    {
        "isotope": "133-Ba",
        "calibration_date": datetime(1979,12,1),
        "initial_activity": 4.8, 
        "half_life": 10.537 * 365.25,  
        "emission_fractions": {53.2:0.02, 356.0:0.62, 383.8:0.09} 
    },
    {
        "isotope": "60-Co",
        "calibration_date": datetime(1979,12,1),
        "initial_activity": 1.9, 
        "half_life": 5.2714 * 365.25,  
        "emission_fractions": {1173.2:1, 1332.5:1} 
    },
]

#detector radius and distance to source
detector_geometry={"detector_radius" : 2.54,  "source_distance" : 15}