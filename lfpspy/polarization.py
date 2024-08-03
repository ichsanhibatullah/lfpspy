"""Class definition for Polarization, a polarization analysis processing module."""

import numpy as np
import math

from .spectral import Spectral

__all__ = ["Polarization"]


class Polarization(Spectral):
    """Class for performing polarization analysis on 3-component passive seismic data.

    Attributes
    ----------
    traces : .seed/.mseed/.miniseed
        3 component sensor, time domain.
    win_length : integer
        Length of each rectangular window, second.

    """
    
    #region additional pre-processing
    @staticmethod
    def merge(traceZ, traceN, traceE):
        zero_list = np.zeros((3, len(traceZ)), dtype=np.float64)
        # East
        zero_list[0] = traceE
        # North
        zero_list[1] = traceN
        # Z
        zero_list[2] = traceZ
        return zero_list
    #end of region: additional pre-processing

    #region polarization analysis
    def pol_ans(self, traceZ, traceN, traceE):
        azi_list, inc_list, rect_list, plan_list, max_eig_list = [], [], [], [], []
        
        for i in range(len(traceZ)):
            merged = self.merge(traceZ[i], traceN[i], traceE[i])
            cov_mat = np.cov(merged)
            eig_vec, eig_val, v = np.linalg.svd(cov_mat)
            
            # Rectilinearity (Montalbetti & Kanasewich,1970)
            rect = 1.0 - np.sqrt(eig_val[1] / eig_val[0])
            
            # Planarity (Jurkevics,1988)
            plan = 1.0 - (2.0 * eig_val[2] / (eig_val[1] + eig_val[0]))

            # Azimuth
            azimuth = math.degrees(math.atan2(eig_vec[0][0], eig_vec[1][0]))

            # Incidence
            eve = np.sqrt(eig_vec[0][0] ** 2 + eig_vec[1][0] ** 2)
            incidence = math.degrees(math.atan2(eve, eig_vec[2][0]))

            # Max Eigenvalue
            max_eigen_val = max(eig_val)
            
            if azimuth < 0.0:
                azimuth = 360.0 + azimuth
            if incidence < 0.0:
                incidence += 180.0
            if incidence > 90.0:
                incidence = 180.0 - incidence
                if azimuth > 180.0:
                    azimuth -= 180.0
                else:
                    azimuth += 180.0
            incidence = 90 - abs(incidence)

            # Appending value of each window
            azi_list.append(azimuth), inc_list.append(incidence), rect_list.append(rect), plan_list.append(plan), max_eig_list.append(max_eigen_val)
        return azi_list, inc_list, rect_list, plan_list, max_eig_list
    #end of region: polarization analysis

    #region whole process
    def calculate(self, lc=1, hc=3.7, fs=100, order=5, d=0.5, overlap=True):
        # Pre-Processing
        Z_norm, N_norm, E_norm = self.norm_amp()
        Z_bp, N_bp, E_bp = self.bp_filter(Z_norm, lc, hc, fs, order), self.bp_filter(N_norm, lc, hc, fs, order), self.bp_filter(E_norm, lc, hc, fs, order)

        # Processing - Polarization Analysis
        split_Z, split_N, split_E = self.split(Z_bp, d, overlap), self.split(N_bp, d, overlap), self.split(E_bp, d, overlap)
        az, inc, rec, plan, max_eig = self.pol_ans(split_Z, split_N, split_E)

        # Making x for plotting
        time = np.linspace(0, int(len(self.Z.data)/((1/self.dt)*self.win_length)), len(split_Z)) #in minute

        return az, inc, rec, plan, max_eig, time
    #end of region: whole process