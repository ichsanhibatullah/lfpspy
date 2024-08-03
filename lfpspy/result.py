"""Class definition for Result, an output generator module."""

import numpy as np
import matplotlib.pyplot as plt

from .spectral import Spectral

__all__ = ["Result"]


class Result(Spectral):
    """Class for plotting and saving the result."""
    
    def __init__(self, traces, win_length):
        super().__init__(traces, win_length)

    #region plotting
    def plot_spectral(self, vh_res, vh_reject, vh_avg, vh_med_curve, vh_std, vh_std_curve,
                      vh_med_f0_bfr, vh_std_f0_bfr, vh_med_f0_aft, vh_std_f0_aft, vh_int_bfr,
                      vh_int_aft, vh_n, psd_res, psd_reject, psd_avg, psd_med_curve,
                      psd_std, psd_std_curve, psd_med_f0_bfr, psd_std_f0_bfr, psd_med_f0_aft,
                      psd_std_f0_aft, psd_iz_bfr, psd_iz_aft, psd_n, low_bound=1,
                      high_bound_vh=6, high_bound_psd=4, xmin=0.5, xmax=10, ymin=0, ymax_vh=10,
                      ymax_psd=1, batch_proc:bool=False):
        # Limiting the Peak, IF the peak is out of the interest range
        frq_lb, frq_hb_vh, frq_hb_psd = self.frq.tolist().index(low_bound), self.frq.tolist().index(high_bound_vh), self.frq.tolist().index(high_bound_psd) #hz      

        # Plotting
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))

        #region plot1
        # Plot Before Rejection
        # Plot Lower Boundary of Integral VHSR 
        ax1.hlines([1], xmin=xmin, xmax=xmax, color="r")

        # Plot Retained
        for i in range(len(vh_res)):
            ax1.plot(self.frq, vh_res[i], color="darkgray", alpha=0.75, lw=0.4) #, label="Retained Windows"

        # Plot f_0 Retained
        for i in range(len(vh_res)):
            ax1.plot(self.frq[vh_res[i].tolist().index(max(vh_res[i]))], max(vh_res[i]), linestyle="", zorder=2, marker='o',
                    markersize=3, markerfacecolor='w', markeredgewidth=0.5, markeredgecolor='k') #, label='$f_{0,i}$'

        # Plot Mean f_0, Mean f_0 ± STD
        ax1.plot([vh_med_f0_bfr]*2, [ymin, ymax_vh], 'k-.')
        ax1.fill([vh_med_f0_bfr-(2*vh_std_f0_bfr)]*2 + [vh_med_f0_bfr+(2*vh_std_f0_bfr)]*2, [ymin, ymax_vh, ymax_vh, ymin], color = "#ff8080") #, label="$μ_{f0}$ ± 2 STD"

        # Plot Mean Curve, Mean Curve ± STD
        ax1.plot(self.frq, vh_avg, 'k') #, label="Mean Curve"
        ax1.plot(self.frq, vh_avg-vh_std, 'k--', alpha=0.75) #, label="Mean Curve ± 1 STD"
        ax1.plot(self.frq, vh_avg+vh_std, 'k--', alpha=0.75)

        # Plot Mean Curve f_0
        ax1.plot(self.frq[vh_avg.tolist().index(max(vh_avg))], max(vh_avg), linestyle="", marker='D', markersize=6,
                markerfacecolor='lime', markeredgewidth=1, markeredgecolor='k',
                label="Peak$_{}$ at ({:.2f} $Hz$, {:.2f})".format("{mc}", self.frq[vh_avg.tolist().index(max(vh_avg))], max(vh_avg))) 

        if vh_avg.tolist().index(max(vh_avg))<frq_lb or vh_avg.tolist().index(max(vh_avg))>frq_hb_vh:
            ax1.plot(self.frq[vh_avg.tolist().index(max(vh_avg[frq_lb:frq_hb_vh]))], max(vh_avg[frq_lb:frq_hb_vh]), linestyle="", marker='D', markersize=6,
                    markerfacecolor='gold', markeredgewidth=1, markeredgecolor='k',
                    label="Peak$_{}$ at ({:.2f} $Hz$, {:.2f})".format("{mc2}", self.frq[vh_avg.tolist().index(max(vh_avg[frq_lb:frq_hb_vh]))], max(vh_avg[frq_lb:frq_hb_vh])))

        # Plot Integral VHSR
        ax1.fill_between(self.frq, 1, vh_avg, color="b", label='Integral-VHSR = {:.5f}'.format(vh_int_bfr), where=(vh_avg>=1)&(self.frq>=1)&(self.frq<=6), alpha=0.75)

        # Plot Legend for this axis
        ax1.legend(loc="upper right")
        #end region plot1

        #region plot2
        # Plot After Rejection
        # Plot Lower Boundary of Integral VHSR 
        ax2.hlines([1], xmin=xmin, xmax=xmax, color="r")

        # Plot Rejected
        for i in range(len(vh_res)):
            if i == 0:
                fl2, = ax2.plot(self.frq, vh_res[i], "c", alpha=0.75, lw=0.3) #color="silver" , label="Rejected Windows"
            else:
                ax2.plot(self.frq, vh_res[i], "c", alpha=0.75, lw=0.3)

        # Plot Retained
        for i in range(len(vh_reject)):
            if i == 0:
                fl1, = ax2.plot(self.frq, vh_reject[i], color="darkgray", alpha=0.75, lw=0.4) #, label="Retained Windows"
            else:
                ax2.plot(self.frq, vh_reject[i], color="darkgray", alpha=0.75, lw=0.4)

        # Plot f_0 Retained
        for i in range(len(vh_reject)):
            if i == 0:
                fl5, = ax2.plot(self.frq[vh_reject[i].tolist().index(max(vh_reject[i]))], max(vh_reject[i]), linestyle="", zorder=2, marker='o',
                        markersize=3, markerfacecolor='w', markeredgewidth=0.5, markeredgecolor='k') #, label='$f_{0,i}$'
            else:
                ax2.plot(self.frq[vh_reject[i].tolist().index(max(vh_reject[i]))], max(vh_reject[i]), linestyle="", zorder=2, marker='o',
                        markersize=3, markerfacecolor='w', markeredgewidth=0.5, markeredgecolor='k')

        # Plot Mean f_0, Mean f_0 ± STD
        ax2.plot([vh_med_f0_aft]*2, [ymin, ymax_vh], 'k-.')
        fl7, = ax2.fill([vh_med_f0_aft-(2*vh_std_f0_aft)]*2 + [vh_med_f0_aft+(2*vh_std_f0_aft)]*2, [ymin, ymax_vh, ymax_vh, ymin], color = "#ff8080")

        # Plot Mean Curve, Mean Curve ± STD
        fl3, = ax2.plot(self.frq, vh_med_curve, 'k') #, label="Mean Curve"
        fl4, = ax2.plot(self.frq, vh_med_curve-vh_std_curve, 'k--', alpha=0.75) #, label="Mean Curve ± 1 STD"
        ax2.plot(self.frq, vh_med_curve+vh_std_curve, 'k--', alpha=0.75)

        # Plot Mean Curve f_0
        fl6, = ax2.plot(self.frq[vh_med_curve.tolist().index(max(vh_med_curve))], max(vh_med_curve), linestyle="", marker='D', markersize=6,
                        markerfacecolor='lime', markeredgewidth=1, markeredgecolor='k',
                        label="Peak$_{}$ at ({:.2f} $Hz$, {:.2f})".format("{mc}", self.frq[vh_med_curve.tolist().index(max(vh_med_curve))], max(vh_med_curve)))

        if vh_med_curve.tolist().index(max(vh_med_curve))<frq_lb or vh_med_curve.tolist().index(max(vh_med_curve))>frq_hb_vh:
            ax2.plot(self.frq[vh_med_curve.tolist().index(max(vh_med_curve[frq_lb:frq_hb_vh]))], max(vh_med_curve[frq_lb:frq_hb_vh]), linestyle="", marker='D', markersize=6,
                    markerfacecolor='gold', markeredgewidth=1, markeredgecolor='k',
                    label="Peak$_{}$ at ({:.2f} $Hz$, {:.2f})".format("{mc2}", self.frq[vh_med_curve.tolist().index(max(vh_med_curve[frq_lb:frq_hb_vh]))], max(vh_med_curve[frq_lb:frq_hb_vh])))

        # Plot Integral VHSR
        ax2.fill_between(self.frq, 1, vh_med_curve, color="b", label='Integral-VHSR = {:.5f}'.format(vh_int_aft), where=(vh_med_curve>=1)&(self.frq>=1)&(self.frq<=6), alpha=0.75)

        # Plot Legend for this axis
        ax2.legend(loc="upper right")
        #end region plot2

        #region plot3
        # Plot Before Rejection
        # Plot Retained
        for i in range(len(psd_res)):
            ax3.plot(self.frq, psd_res[i], color="darkgray", alpha=0.75, lw=0.4) #, label="Retained Windows"

        # Plot f_0 Retained
        for i in range(len(psd_res)):
            ax3.plot(self.frq[psd_res[i].tolist().index(max(psd_res[i]))], max(psd_res[i]), linestyle="", zorder=2, marker='o',
                    markersize=3, markerfacecolor='w', markeredgewidth=0.5, markeredgecolor='k') #, label='$f_{0,i}$'

        # Plot Mean f_0, Mean f_0 ± STD
        ax3.plot([psd_med_f0_bfr]*2, [ymin, ymax_psd], 'k-.')
        ax3.fill([psd_med_f0_bfr-(2*psd_std_f0_bfr)]*2 + [psd_med_f0_bfr+(2*psd_std_f0_bfr)]*2, [ymin, ymax_psd, ymax_psd, ymin], color = "#ff8080") #, label="$μ_{f0}$ ± 2 STD"

        # Plot Mean Curve, Mean Curve ± STD
        ax3.plot(self.frq, psd_avg, 'k') #, label="Mean Curve"
        ax3.plot(self.frq, psd_avg-psd_std, 'k--', alpha=0.75) #, label="Mean Curve ± 1 STD"
        ax3.plot(self.frq, psd_avg+psd_std, 'k--', alpha=0.75)

        # Plot Mean Curve f_0
        ax3.plot(self.frq[psd_avg.tolist().index(max(psd_avg))], max(psd_avg), linestyle="", marker='D', markersize=6,
                markerfacecolor='lime', markeredgewidth=1, markeredgecolor='k',
                label="Peak$_{}$ at ({:.2f} $Hz$, {:.2f} $Count^2/Hz$)".format("{mc}", self.frq[psd_avg.tolist().index(max(psd_avg))], max(psd_avg)))

        if psd_avg.tolist().index(max(psd_avg))<frq_lb or psd_avg.tolist().index(max(psd_avg))>frq_hb_psd:
            ax3.plot(self.frq[psd_avg.tolist().index(max(psd_avg[frq_lb:frq_hb_psd]))], max(psd_avg[frq_lb:frq_hb_psd]), linestyle="", marker='D', markersize=6,
                    markerfacecolor='gold', markeredgewidth=1, markeredgecolor='k',
                    label="Peak$_{}$ at ({:.2f} $Hz$, {:.2f} $Count^2/Hz$)".format("{mc2}", self.frq[psd_avg.tolist().index(max(psd_avg[frq_lb:frq_hb_psd]))], max(psd_avg[frq_lb:frq_hb_psd])))

        # Plot Integral PSD-Z
        ax3.fill_between(self.frq, 0, psd_avg, color="b", label='PSD-IZ = {:.5f}'.format(psd_iz_bfr), where=(psd_avg>=0)&(self.frq>=1)&(self.frq<=4), alpha=0.75)

        # Plot Legend for this axis
        ax3.legend(loc="upper right")
        #end region plot3

        #region plot4
        # Plot After Rejection
        # Plot Rejected
        for i in range(len(psd_res)):
            ax4.plot(self.frq, psd_res[i], "c", alpha=0.75, lw=0.3) #color="silver" , label="Rejected Windows"
            
        # Plot Retained
        for i in range(len(psd_reject)):
            ax4.plot(self.frq, psd_reject[i], color="darkgray", alpha=0.75, lw=0.4) #, label="Retained Windows"

        # Plot f_0 Retained
        for i in range(len(psd_reject)):
            ax4.plot(self.frq[psd_reject[i].tolist().index(max(psd_reject[i]))], max(psd_reject[i]), linestyle="", zorder=2, marker='o',
                    markersize=3, markerfacecolor='w', markeredgewidth=0.5, markeredgecolor='k') #, label='$f_{0,i}$'

        # Plot Mean f_0, Mean f_0 ± STD
        ax4.plot([psd_med_f0_aft]*2, [ymin, ymax_psd], 'k-.')
        fl8, = ax4.fill([psd_med_f0_aft-(2*psd_std_f0_aft)]*2 + [psd_med_f0_aft+(2*psd_std_f0_aft)]*2, [ymin, ymax_psd, ymax_psd, ymin], color = "#ff8080")

        # Plot Mean Curve, Mean Curve ± STD
        ax4.plot(self.frq, psd_med_curve, 'k') #, label="Mean Curve"
        ax4.plot(self.frq, psd_med_curve-psd_std_curve, 'k--', alpha=0.75) #, label="Mean Curve ± 1 STD"
        ax4.plot(self.frq, psd_med_curve+psd_std_curve, 'k--', alpha=0.75)

        # Plot Mean Curve f_0
        ax4.plot(self.frq[psd_med_curve.tolist().index(max(psd_med_curve))], max(psd_med_curve), linestyle="", marker='D', markersize=6,
                markerfacecolor='lime', markeredgewidth=1, markeredgecolor='k',
                label="Peak$_{}$ at ({:.2f} $Hz$, {:.2f} $Count^2/Hz$)".format("{mc}", self.frq[psd_med_curve.tolist().index(max(psd_med_curve))], max(psd_med_curve)))

        if psd_med_curve.tolist().index(max(psd_med_curve))<frq_lb or psd_med_curve.tolist().index(max(psd_med_curve))>frq_hb_psd:
            ax4.plot(self.frq[psd_med_curve.tolist().index(max(psd_med_curve[frq_lb:frq_hb_psd]))], max(psd_med_curve[frq_lb:frq_hb_psd]), linestyle="", marker='D', markersize=6,
                    markerfacecolor='gold', markeredgewidth=1, markeredgecolor='k',
                    label="Peak$_{}$ at ({:.2f} $Hz$, {:.2f} $Count^2/Hz$)".format("{mc2}", self.frq[psd_med_curve.tolist().index(max(psd_med_curve[frq_lb:frq_hb_psd]))], max(psd_med_curve[frq_lb:frq_hb_psd])))

        # Plot Integral PSD-Z
        # Plot Integral PSD-Z
        ax4.fill_between(self.frq, 0, psd_med_curve, color="b", label='PSD-IZ = {:.5f}'.format(psd_iz_aft), where=(psd_med_curve>=0)&(self.frq>=1)&(self.frq<=4), alpha=0.75)

        # Plot Legend for this axis
        ax4.legend(loc="upper right")
        #end region plot4

        #region limiting plot
        ax1.set_xlim(xmin, xmax), ax2.set_xlim(xmin, xmax), ax3.set_xlim(xmin, xmax), ax4.set_xlim(xmin, xmax)
        ax1.set_ylim(ymin, ymax_vh), ax2.set_ylim(ymin, ymax_vh), ax3.set_ylim(ymin, ymax_psd), ax4.set_ylim(ymin, ymax_psd)
        #end region limiting plot

        #region label axis
        ax1.set_xlabel("Frequency ($Hz$)"), ax2.set_xlabel("Frequency ($Hz$)"), ax3.set_xlabel("Frequency ($Hz$)"), ax4.set_xlabel("Frequency ($Hz$)")
        ax1.set_ylabel("VHSR Amplitude"), ax2.set_ylabel("VHSR Amplitude"),
        ax3.set_ylabel("Power Spectral Density ($Count^2/Hz$)"), ax4.set_ylabel("Power Spectral Density ($Count^2/Hz$)")
        #end region label axis

        #region title axis
        ax1.set_title('VHSR - Before Rejection'), ax2.set_title('VHSR - After Rejection'), ax3.set_title('PSD-Z - Before Rejection'), ax4.set_title('PSD-Z - After Rejection')
        #end region title axis

        fig.suptitle("Spectral Analysis: Station {}".format(self.Z.stats.station), fontsize = 14)
        fig.legend([fl1, fl2, fl3, fl4, fl5, fl6, fl7, fl8], ["Retained Windows", "Rejected Windows", "Median Curve", "Median Curve ± 1 STD",
                                                        "Peak Each Window", "Peak Average Curve", "$μ$ Peak ± {} STD (VHSR)".format(vh_n), "$μ$ Peak ± {} STD (PSD-Z)".format(psd_n)],
                loc='upper center', bbox_to_anchor=(0.5, 0.01), shadow=True, ncol=4)
        #fig.gca().add_artist(mainlegend)
        fig.tight_layout()
        
        plt.savefig("result/Station {}-Spectral Analysis.png".format(self.Z.stats.station), bbox_inches='tight')

        if batch_proc == False:
            plt.show()
        else:
            plt.close()
    
    def plot_polar(self, azimuth, incidence, rectilinearity, planarity, max_eig, time, batch_proc:bool=False):        
        # Plotting
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

        ax1.plot(time, incidence, ".k")
        ax1.axhline(np.sum(incidence)/len(incidence), color="r")        
        ax1.set_xlabel("Time (min)")
        ax1.set_ylabel("Dip (°)")
        ax1.set_ylim((0, 90))

        ax2.plot(time, azimuth, ".b")
        ax2.set_xlabel("Time (min)")
        ax2.set_ylabel("Azimuth (°)")
        ax2.set_ylim((0, 360))

        ax3.plot(time, rectilinearity, "-r")
        ax3.axhline(np.sum(rectilinearity)/len(rectilinearity), color="k")
        ax3.set_xlabel("Time (min)")
        ax3.set_ylabel("Rectilinearity")
        ax3.set_ylim((0, 1))

        ax4.plot(time, max_eig, "-g")
        ax4.set_xlabel("Time (min)")
        ax4.set_ylabel("Largest Eigenvalue")
        ax4.set_ylim(0)

        fig.suptitle("Polarization Analysis: Station {}".format(self.Z.stats.station), fontsize = 14)
        fig.tight_layout()

        plt.savefig("result/Station {}-Polarization Analysis.png".format(self.Z.stats.station))
        
        if batch_proc == False:
            plt.show()
        else:
            plt.close()   
    #end of region: plotting