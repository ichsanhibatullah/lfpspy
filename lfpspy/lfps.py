"""Class definition for LFPS, a LFPS processing module."""

import obspy
import numpy as np

from scipy.signal import butter, filtfilt, detrend
from scipy.signal.windows import tukey
import scipy.fftpack as fftpack
from obspy.signal import konnoohmachismoothing as sko

__all__ = ["LFPS"]


class LFPS():
    """Class for creating and manipulating 3-component sensor objects.

    Attributes
    ----------
    ns : TimeSeries
        North-south component, time domain.
    soon : xx Object
        xx soon, xx domain.

    """

    def __init__(self, traces, win_length):
        #region full traces
        self.traces = obspy.read(traces)
        #region per trace
        self.Z = self.traces[0]
        self.N = self.traces[1]
        self.E = self.traces[2]
        #region params
        self.win_length = win_length
        self.dt = self.Z.stats.delta
        self.npts = int(win_length*self.Z.stats.sampling_rate)
        self.nfrqs = int(self.npts/2)
        self.frq = np.abs(np.fft.fftfreq(self.npts, self.dt))[:self.nfrqs]
        #end of region
    
    #region pre-processing
    def norm_amp(self):
        temp_lis=[]
        data=[]
        temp_lis.append(self.Z.data), temp_lis.append(self.N.data), temp_lis.append(self.E.data)
        for x in range(len(temp_lis)):
            for y in range(len(temp_lis[x])):
                data.append(temp_lis[x][y])
        xmin = min(data)
        xmax = max(data)
        Z = 2*((self.Z-xmin)/(xmax-xmin))-1
        N = 2*((self.N-xmin)/(xmax-xmin))-1
        E = 2*((self.E-xmin)/(xmax-xmin))-1
        return Z, N, E
    
    def norm_one_amp(self):
        return 50*((self.Z.data-min(self.Z.data))/(max(self.Z.data)-min(self.Z.data)))-25
    
    @staticmethod
    def bp_filter(normalized_amp, lc=0.5, hc=10, fs=100, order=5):
        fnyq = 0.5*fs
        
        # Create the filter coefficients
        b, a = butter(order,[lc/fnyq, hc/fnyq], btype='bandpass')
        # Apply the filter to your signal (replace 'signal' with your actual data)
        filter = filtfilt(b, a, normalized_amp, padlen=3*(max(len(b), len(a))-1))
        return filter
    
    def split(self, filtered_signal, d=0.5, overlap=True):
        if overlap == False:
            return np.split(filtered_signal, self.Z.stats.npts/(self.npts))
        else:
            w = []
            for i in range(0, len(filtered_signal), round(self.npts*(1-d))):
                w.append(filtered_signal[i:i+self.npts])
            while len(w[-1]) != self.npts:
                del w[-1]
            return w

    @staticmethod
    def cosine_taper(de_trended, width=0.1): # Geopsy default of 0.05 is equal to 0.1
        return de_trended * tukey(len(de_trended[0]), alpha=width) 
        # width `0` is equal to a rectangular and `1` a Hann window. 

    def preprocessing(self, width=0.1, lc=0.5, hc=10, fs=100, order=5, d=0.5, overlap=True):
        normZ, normN, normE = self.norm_amp()
        filterZ, filterN, filterE = self.bp_filter(normZ, lc, hc, fs, order), self.bp_filter(normN, lc, hc, fs, order), self.bp_filter(normE, lc, hc, fs, order)
        windowZ, windowN, windowE = self.split(filterZ, d), self.split(filterN, d), self.split(filterE, d)
        de_trendZ, de_trendN, de_trendE = detrend(windowZ), detrend(windowN), detrend(windowE)
        taperZ, taperN, taperE = self.cosine_taper(de_trendZ, width), self.cosine_taper(de_trendN, width), self.cosine_taper(de_trendE, width)

        one_normZ = self.norm_one_amp()
        one_filterZ = self.bp_filter(one_normZ, lc, hc, fs, order)
        one_windowZ = self.split(one_filterZ, d)
        one_de_trendZ = detrend(one_windowZ)
        one_taperZ = self.cosine_taper(one_de_trendZ, width)
        return taperZ, taperN, taperE, one_taperZ
    #end of region: pre-processing

    #region processing: fft and psd
    def make_fft(self, preprocessed):
        fft, sliced_fft = np.empty(shape=[len(preprocessed), self.npts]), np.empty(shape=[len(preprocessed), self.nfrqs])

        for i in range(len(preprocessed)):
            fft[i] = 2/self.npts*np.abs(fftpack.fft(preprocessed[i]))
            sliced_fft[i] = fft[i][:self.nfrqs]
        return sliced_fft
    
    def make_psd(self, preprocessed):
        psd, sliced_psd = np.empty(shape=[len(preprocessed), self.npts]), np.empty(shape=[len(preprocessed), self.nfrqs])
        
        for iter in range(len(preprocessed)):
            psd[iter] = (np.abs(fftpack.fft(preprocessed[iter])))**2/(self.npts*1/self.dt)
            sliced_psd[iter] = psd[iter][:self.nfrqs]
        return sliced_psd
    
    @staticmethod
    def signal_stack(sliced_list):
        copy_list = sliced_list.copy()
        base, w_start = 0, 1
        
        for w_start in range(len(copy_list)):
            copy_list[base] += copy_list[w_start]
            copy_list[base] / 2
            w_start += 1
        return copy_list[base]
    
    @staticmethod
    def signal_stat(sliced_list, mode="median"):  
        if mode == "median":
            median_list = np.median(sliced_list, axis=0)
            return median_list
        elif mode == "mean":
            mean_list = np.mean(sliced_list, axis=0)
            return mean_list
        else:
            msg = "Method isn't implemented yet."
            raise Exception(msg)
    
    def fft_psd(self, preprocessed, preprocessed2=None, make_psd=False):
        if make_psd == False:
            fft = self.make_fft(preprocessed)
            avg = self.signal_stat(fft, mode="median")
            return fft, avg
        else: # then make psd
            # fft
            fft = self.make_fft(preprocessed)
            avg_fft = self.signal_stat(fft, mode="median")

            # psd
            std_psd = []

            psd = self.make_psd(preprocessed2)
            smooth_psd = sko.konno_ohmachi_smoothing(psd, self.frq, normalize=True)
            for i in range(len(smooth_psd[0])):
                std_c = np.std(smooth_psd[:,i])
                std_psd.append(std_c)
            avg_psd = self.signal_stat(psd, mode="median")
            return fft, smooth_psd, avg_fft, avg_psd, std_psd
    #end of region: fft and psd

    #region processing: spectral ratio
    @staticmethod
    def combine_hor(n, e, method="geom-mean"):
        if method == "geom-mean":
            hor_gm = np.sqrt(n * e)
            return hor_gm
        elif method == "sqd-avg":
            hor_sa = np.sqrt((n * n + e * e)/2)
            return hor_sa
        else:
            msg = "Method isn't implemented yet."
            raise Exception(msg)
        
    def spectral_ratio(self, procZ, procN, procE, combine_method="geom-mean", make_vhsr=True, smooth_ko=True):
        horizontal = self.combine_hor(procN, procE, combine_method)
        vertical = procZ.copy()
        smooth_ver = sko.konno_ohmachi_smoothing(vertical, self.frq, normalize=True)
        smooth_hor = sko.konno_ohmachi_smoothing(horizontal, self.frq, normalize=True)
        std = []
        
        if make_vhsr == True:
            if smooth_ko == True:
                smooth_vhsr = smooth_ver / smooth_hor
                for i in range(len(smooth_vhsr[0])):
                    std_c = np.std(smooth_vhsr[:,i])
                    std.append(std_c)
                return smooth_vhsr, std
            else:
                vhsr = vertical / horizontal
                for i in range(len(vhsr[0])):
                    std_c = np.std(vhsr[:,i])
                    std.append(std_c)
                return vhsr, std
        else: # then make hvsr
            if smooth_ko == True:
                smooth_hvsr = smooth_hor / smooth_ver
                for i in range(len(smooth_hvsr[0])):
                    std_c = np.std(smooth_hvsr[:,i])
                    std.append(std_c)
                return smooth_hvsr, std
            else:
                hvsr = horizontal / vertical
                for i in range(len(hvsr[0])):
                    std_c = np.std(hvsr[:,i])
                    std.append(std_c)            
                return hvsr, std
    
    def make_zero(self, spectrum, low_cut, high_cut):
        copy_list = spectrum.copy()
        
        frq_low = self.frq.tolist().index(low_cut)
        frq_high = self.frq.tolist().index(high_cut)
        
        for i in range(0, int(frq_low)):
            copy_list[:,i] = 0
        for j in range(int(frq_high)+1, len(copy_list[0])):
            copy_list[:,j] = 0
        return copy_list
    #end of region: spectral ratio

    #region integral area
    @staticmethod
    def f(x):
        return x/x
    
    def frq_min_f(self, min_f=1):
        return self.frq.tolist().index(min_f)
    
    def frq_max_f_vhsr(self, max_f=6):
        return self.frq.tolist().index(max_f)
    
    def frq_max_f_psd(self, max_f=4):
        return self.frq.tolist().index(max_f)

    def custom_f(self, spectrum, min_f=1, max_f=6):
        empty_list = []

        for i in range(self.frq_min_f(min_f), self.frq_max_f_vhsr(max_f)):
            if spectrum[i] < 1.0:
                empty_list.append(1.0)
            else:
                empty_list.append(spectrum[i])
        return empty_list
    
    def integral_area_vhsr(self, spectrum, min_f=1, max_f=6):
        integral = np.trapz(self.custom_f(spectrum),self.frq[self.frq_min_f(min_f):self.frq_max_f_vhsr(max_f)]) - np.trapz(self.f(self.frq[self.frq_min_f(min_f):self.frq_max_f_vhsr(max_f)]), self.frq[self.frq_min_f(min_f):self.frq_max_f_vhsr(max_f)])
        return integral
    
    def integral_area_psd(self, spectrum, min_f=1, max_f=4):
        integral = np.trapz(spectrum[self.frq_min_f(min_f):self.frq_max_f_psd(max_f)], self.frq[self.frq_min_f(min_f):self.frq_max_f_psd(max_f)])
        return integral
    #end of region: integral area

    #region whole process
    def calculate(self, width=0.1, lc=0.5, hc=10, fs=100, order=5, d=0.5, overlap=True):
        # Preprocess data
        prepZ, prepN, prepE, one_norm_prepZ = self.preprocessing(width, lc, hc, fs, order, d, overlap)

        # Performing FFT on each component and PSD only for Z component (PSD with smoothing and std available, FFT later in spectral ratio)
        fftZ, smooth_psd, fftZ_avg, avg_psd, std_psd = self.fft_psd(prepZ, one_norm_prepZ, make_psd=True)
        fftN, fftN_avg = self.fft_psd(prepN)
        fftE, fftE_avg = self.fft_psd(prepE)

        # Performing VHSR with smoothing konno ohmachi
        vhsr_smooth, std_vhsr = self.spectral_ratio(fftZ, fftN, fftE)
        vhsr_smooth_zero = self.make_zero(vhsr_smooth, lc, hc)
        avg_vhsr = self.signal_stat(vhsr_smooth_zero, mode="median")

        # Integral area calculation of desired range; before rejection
        vhsr_integral_bfr = self.integral_area_vhsr(avg_vhsr)
        psd_iz_bfr = self.integral_area_psd(avg_psd)
        return vhsr_smooth_zero, avg_vhsr, std_vhsr, vhsr_integral_bfr, smooth_psd, avg_psd, std_psd, psd_iz_bfr #, vhsr_bfr_mean_f0, vhsr_bfr_std_f0, psd_bfr_mean_f0, psd_bfr_std_f0
