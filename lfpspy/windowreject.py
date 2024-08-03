"""Class definition for WindowReject, a statistical window rejection algorithm module proposed by Cox et. al. (2020)."""

import numpy as np

from .spectral import Spectral

__all__ = ["WindowReject"]


class WindowReject(Spectral):
    """Class for applying window rejection on frequency domain to a passive seismic data.

    Attributes
    ----------
    traces : .seed/.mseed/.miniseed
        3 component sensor, time domain.
    win_length : integer
        Length of each rectangular window, second.

    """

    def __init__(self, traces, win_length):
        super().__init__(traces, win_length)

    #region little function
    @staticmethod
    def skip_zero(list):
        # Skip zero values
        return list[~np.all(list == 0, axis=1)]

    @staticmethod
    def mean_factory(values, cond):
        # Factory: Calculate mean values
        return np.mean(values, where=cond)

    @staticmethod
    def std_factory(values, cond):
        # Factory: Calculate standard deviation values
        return np.std(values, ddof=1, where=cond)

    @staticmethod
    def nth_std_factory(n, mean, std):
        # Factory: Calculate standard deviation values from the mean of each peak window/spectra
        return (mean + n*std)

    def peak_search_frq(self, val):
        # Searching peak frequency for each window/spectra
        peak=[]
        for i in range(len(val)):
            peak_w = self.frq[val[i].tolist().index(max(val[i]))]
            peak.append(peak_w)
        return peak

    def mean_f0_frq(self, val, cond):
        # Calculate mean values
        peak = self.peak_search_frq(val)
        return self.mean_factory(peak, cond)

    def std_f0_frq(self, val, cond):
        # Calculate standard deviation values
        peak = self.peak_search_frq(val)
        return self.std_factory(peak, cond)

    def mc_peak_frq(self, all_window):
        # Searching peak frequency of mean curve
        cop = []
        skip = self.skip_zero(all_window)

        if len(skip) == 0:
            raise Exception('All window is rejected. Please check your data again, or raise the value of n.')
        else:
            for i in range(len(skip[0])):
                mean = np.mean(skip[:,i])
                cop.append(mean)
        
        max_mean = max(cop)
        max_mean_frq = self.frq[cop.index(max_mean)]
        return max_mean_frq
    #end of region: little function

    #region whole process
    def calculate(self, spectrum, n=2, max_iter=50):
        # Frequency-domain window-rejection algorithm
        start_window = spectrum.copy()
        sum_start = len(start_window)

        mean_f0_bfr = self.mean_f0_frq(start_window, np.ones(len(start_window), dtype=bool))
        std_f0_bfr = self.std_f0_frq(start_window, np.ones(len(start_window), dtype=bool))
        
        for c_iter in range(1, max_iter+1):        
            iter_window = self.skip_zero(start_window)
            valid_indices = np.ones(len(iter_window), dtype=bool)
            
            # defining beginning params
            mean_f0_b = self.mean_f0_frq(iter_window, valid_indices)
            std_f0_b = self.std_f0_frq(iter_window, valid_indices)
            mc_peak_frq_b = self.mc_peak_frq(iter_window)
            d_b = abs(mean_f0_b - mc_peak_frq_b)

            lower_bound = self.nth_std_factory(-n, mean_f0_b, std_f0_b)
            upper_bound = self.nth_std_factory(+n, mean_f0_b, std_f0_b)
            
            
            peak = self.peak_search_frq(iter_window)
            for i in range(len(valid_indices)):
                if not peak[i] > lower_bound or not peak[i] < upper_bound:
                    valid_indices[i] = False
                if valid_indices[i] == False:
                    iter_window[i] = np.zeros(len(iter_window[i]))
            

            mean_f0_aft = self.mean_f0_frq(iter_window, valid_indices)
            std_f0_aft = self.std_f0_frq(iter_window, valid_indices)
            mc_peak_frq_aft = self.mc_peak_frq(iter_window)
            d_aft = abs(mean_f0_aft - mc_peak_frq_aft)
            
            d_diff = abs(d_aft - d_b)/d_b
            s_diff = abs(std_f0_aft - std_f0_b)

            # del start_window
            start_window = iter_window
            # print("Iteration-{}: Initial Windows = {}, Remaining Windows = {}, Windows Rejected = {}".format(c_iter, len(start_window), np.count_nonzero(valid_indices),
            #                                                                                                 len(start_window)-np.count_nonzero(valid_indices)))
            
            final_window = self.skip_zero(iter_window)
            final_mean = []
            final_std = []
            
            for i in range(len(final_window[0])):
                mean_c = np.mean(final_window[:,i])
                std_c = np.std(final_window[:,i])
                final_mean.append(mean_c)
                final_std.append(std_c)
            
            if d_b == 0 or std_f0_b == 0 or std_f0_aft == 0:
                print("Performed {} iterations, returning b/c 0 values".format(c_iter))
                print("Initial Windows = {}, Remaining Windows = {}, Windows Rejected = {}".format(sum_start, np.count_nonzero(valid_indices), sum_start-np.count_nonzero(valid_indices)))
                return final_window, np.array(final_mean), final_std, mean_f0_bfr, std_f0_bfr, mean_f0_aft, std_f0_aft, n, np.count_nonzero(valid_indices)
            
            if (d_diff < 0.01) and (s_diff < 0.01):
                print("Performed {} iterations, returning b/c rejection converged".format(c_iter))
                print("Initial Windows = {}, Remaining Windows = {}, Windows Rejected = {}".format(sum_start, np.count_nonzero(valid_indices), sum_start-np.count_nonzero(valid_indices)))                
                return final_window, np.array(final_mean), final_std, mean_f0_bfr, std_f0_bfr, mean_f0_aft, std_f0_aft, n, np.count_nonzero(valid_indices)

            if c_iter == max_iter:
                print("Ouch! didn't pass the stopping criteria 'til the end of iteration. Check your data or adjust the parameters!")
                print("Initial Windows = {}, Remaining Windows = {}, Windows Rejected = {}".format(sum_start, np.count_nonzero(valid_indices), sum_start-np.count_nonzero(valid_indices)))                
                return final_window, np.array(final_mean), final_std, mean_f0_bfr, std_f0_bfr, mean_f0_aft, std_f0_aft, n, np.count_nonzero(valid_indices)
    #end of region: whole process