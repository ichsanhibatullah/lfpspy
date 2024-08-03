import pandas as pd
import numpy as np

#region saving calculation
def create_df_spectral():
    df = pd.DataFrame(columns=['Station','V/H Integral', 'V/H Max Ratio', 'Frequency Max Amplitude VHSR (Hz)',
                                'PSD-IZ', 'PSD-Z Max Amplitude (Count^2/Hz)', 'Frequency Max Amplitude PSD-Z (Hz)',
                                'VHSR window accepted', 'PSD-Z window accepted'])
    return df


def create_df_polar():
    df = pd.DataFrame(columns=['Station', 'Azimuth: Stability', 'Incidence Angle (°)',
                                'Rectilinearity', 'Planarity', 'Max Eigenvalue'])
    return df

def list_spectral(df, comp, frq, frq_bb, frq_ba_vh, frq_ba_psd, vh_int_aft,
                    vh_med_curve, vh_acc, psd_iz_aft, psd_med_curve, psd_acc):
    df.loc[int(comp.stats.station),
                            'Station'] = '{}'.format(comp.stats.station)
    df.loc[int(comp.stats.station),
                            'V/H Integral'] = vh_int_aft
    df.loc[int(comp.stats.station),
                            'V/H Max Ratio'] = max(vh_med_curve[frq_bb:frq_ba_vh])
    df.loc[int(comp.stats.station),
                            'Frequency Max Amplitude VHSR (Hz)'] = frq[vh_med_curve.tolist().index(max(vh_med_curve[frq_bb:frq_ba_vh]))]        
    df.loc[int(comp.stats.station),
                            'PSD-IZ'] = psd_iz_aft
    df.loc[int(comp.stats.station),
                            'PSD-Z Max Amplitude (Count^2/Hz)'] = max(psd_med_curve[frq_bb:frq_ba_psd])
    df.loc[int(comp.stats.station),
                            'Frequency Max Amplitude PSD-Z (Hz)'] = frq[psd_med_curve.tolist().index(max(psd_med_curve[frq_bb:frq_ba_psd]))]
    df.loc[int(comp.stats.station),
                            'VHSR window accepted'] = vh_acc
    df.loc[int(comp.stats.station),
                            'PSD-Z window accepted'] = psd_acc

def list_polar(df, comp, azimuth, incidence, rectilinearity, planarity, max_eig):
    df.loc[int(comp.stats.station),
                            'Station'] = '{}'.format(comp.stats.station)
    df.loc[int(comp.stats.station),
                            'Azimuth: Stability'] = np.var(azimuth)
    df.loc[int(comp.stats.station),
                            'Incidence Angle (°)'] = np.mean(incidence)
    df.loc[int(comp.stats.station),
                            'Rectilinearity'] = np.mean(rectilinearity)
    df.loc[int(comp.stats.station),
                            'Planarity'] = np.mean(planarity)
    df.loc[int(comp.stats.station),
                            'Max Eigenvalue'] = np.mean(max_eig)        

def to_excel_spectral(df, filename:str="Spectral Analysis Calculation-raw python"):
    df.reset_index(drop=True, inplace=True)
    df.to_excel("result/" + filename + ".xlsx", index=False)

# def to_csv_spectral(df, filename:str="Spectral Analysis Calculation-raw python"):
#     df.reset_index(drop=True, inplace=True)
#     df.to_csv("result/" + filename + ".csv", index=False)

def to_excel_polar(df, filename="Polarization Analysis Calculation-raw python"):
    df.reset_index(drop=True, inplace=True)
    df.to_excel("result/" + filename + ".xlsx", index=False)

# def to_csv_polar(df, filename="Polarization Analysis Calculation-raw python"):
#     df.reset_index(drop=True, inplace=True)
#     df.to_excel("result/" + filename + ".csv", index=False, engine='python')             
#end of region: saving calculation