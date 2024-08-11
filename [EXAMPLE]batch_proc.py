# Import module needed
import os, time

from lfpspy import Spectral, WindowReject, Polarization, Result
from lfpspy.utils import *

st = time.time()

# Define file (per station) or folder (all station), window length, and other parameters that will be used
dir = "C:/Users/PC-MSI/Downloads/five_raw_data_example"
window, overlap = 60, 0.5 # window length in seconds; how much we want our data overlap per window
n_input_vh, n_input_psd, max_iter = 1.8, 2, 50 # n parameter to decide how much we want to reject; depending on variance of the data
lb_vh, hb_vh, lb_psd, hb_psd = 1, 6, 1, 4 #in Hz; lower bound (each) and upper bound (each) for limiting the zone of interest in vhsr and psd

# Create data frame to save the calculation later
df = create_df_spectral()

# Listing for each file in the directory given
for fname in os.listdir(dir):
    if fname.endswith('.mseed'):
        # Instantiating object
        fname = os.path.join(dir, fname)
        spc, wr, res = Spectral(fname, win_length=window), WindowReject(fname, window), Result(fname, window)

        # Performing spectral ratio and power spectral density processing
        vh_res, vh_avg, vh_std, vh_int_bfr, psd_res, psd_avg, psd_std, psd_iz_bfr = spc.calculate(
            lb_vh, hb_vh, lb_psd, hb_psd, width=0.1, lc=0.5, hc=10, fs=100, order=5, d=overlap, overlap=True)

        # Applying window rejection to spectral ratio and power spectral density result; Achieving better result
        # Performing Frequency-domain Window-rejection algorithm (Cox et.al., 2020)
        print('---------------------------{}---------------------------'.format(fname.split("\\")[-1]))
        vh_reject, vh_med_curve, vh_std_curve, vh_med_f0_bfr, vh_std_f0_bfr, vh_med_f0_aft, vh_std_f0_aft, vh_n, vh_acc = wr.calculate(
            vh_res, n=n_input_vh, max_iter=max_iter)
        vh_int_aft = spc.integral_area_vhsr(vh_med_curve, lb_vh, hb_vh)
        psd_reject, psd_med_curve, psd_std_curve, psd_med_f0_bfr, psd_std_f0_bfr, psd_med_f0_aft, psd_std_f0_aft, psd_n, psd_acc = wr.calculate(
            psd_res, n=n_input_psd, max_iter=max_iter)
        psd_iz_aft = spc.integral_area_psd(psd_med_curve, lb_psd, hb_psd)
        print('------------------------------------------------------------------------------')

        # Plotting and saving the result for each iteration
        res.plot_spectral(vh_res, vh_reject, vh_avg, vh_med_curve, vh_std, vh_std_curve,
                          vh_med_f0_bfr, vh_std_f0_bfr, vh_med_f0_aft, vh_std_f0_aft, vh_int_bfr,
                          vh_int_aft, vh_n, psd_res, psd_reject, psd_avg, psd_med_curve,
                          psd_std, psd_std_curve, psd_med_f0_bfr, psd_std_f0_bfr, psd_med_f0_aft,
                          psd_std_f0_aft, psd_iz_bfr, psd_iz_aft, psd_n, low_bound_vh=lb_vh, high_bound_vh=hb_vh,
                          low_bound_psd=lb_psd, high_bound_psd=hb_psd, xmin=0.5, xmax=10, ymin=0, ymax_vh=10,
                          ymax_psd=1, batch_proc=True) # Change batch_proc to True for batch processing
        list_spectral(df, spc.Z, spc.frq, spc.frq.tolist().index(lb_vh), spc.frq.tolist().index(hb_vh),
                      spc.frq.tolist().index(lb_psd), spc.frq.tolist().index(hb_psd), vh_int_aft,
                      vh_med_curve, vh_acc, psd_iz_aft, psd_med_curve, psd_acc)
        
# Exporting data frame containing calculation result from all iteration to .xlsx file
to_excel_spectral(df, "Spectral Analysis Calculation-raw python")

########################################   POLARIZATION ANALYSIS   ########################################

# Create data frame to save the calculation later
df = create_df_polar()

# Listing for each file in the directory given
for fname in os.listdir(dir):
    if fname.endswith('.mseed'):
        # Instantiating object
        fname = os.path.join(dir, fname)
        plr, res = Polarization(fname, win_length=window), Result(fname, window)

        # Performing polarization analysis on raw passive seismic data
        azimuth, incidence, rectilinearity, planarity, max_eig, xtime = plr.calculate(lc=0.5, hc=10, fs=100, order=5, d=overlap, overlap=True)

        # Plotting and saving the result of each iteration
        res.plot_polar(azimuth, incidence, rectilinearity, planarity, max_eig, xtime, batch_proc=True) # Change batch_proc to True for batch processing
        list_polar(df, plr.Z, azimuth, incidence, rectilinearity, planarity, max_eig)

# Exporting data frame containing calculation result from all iteration to .xlsx file
to_excel_polar(df, "batch_test/Polarization Analysis Calculation-raw python")

elapsed_time = time.time() - st
print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))