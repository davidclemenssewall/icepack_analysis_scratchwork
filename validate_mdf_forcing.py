"""
Tools for comparing Icepack runs with each other and observations

"""

from pathlib import Path
import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime

from icepacktools import load_icepack_hist, plot_hist_var, plot_forc_var
from icepacktools import plot_ice_var, plot_handler

# Load history output
ip_dirs_path = os.path.join(Path.home(), "icepack-dirs")
run_dict = {"mos_v2_kazr":None,
            "mos_v2_init":None,
            }
trcr_dict = {19: 'apnd',
             20: 'hpnd',
             21: 'ipnd'}
trcrn_dict = {19: 'apndn',
              20: 'hpndn',
              21: 'ipndn'}

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value, 
                                       trcr_dict=trcr_dict, 
                                       trcrn_dict=trcrn_dict)

# Create ice thickness change
for value in hist_dict.values():
    value['vice_chg'] = value['vice'] - value['vice'].isel(time=0)
    value['pndaspect'] = value['hpnd'] / value['apnd']
    value['hin'] = value['vicen'] / value['aicen']

# Load forcing
forcing_path = os.path.join(ip_dirs_path, "input", "Icepack_data",
                            "forcing", "MDF")
atm_filename = 'MOSAiC_atm_drift1_stakes_snow_syi_MDF_20191015_20200731.nc'
ds_atm = xr.open_dataset(os.path.join(forcing_path, atm_filename))
ocn_filename = 'MOSAiC_ocn_drift1_MDF_20191006_20200731.nc'
ds_ocn = xr.open_dataset(os.path.join(forcing_path, ocn_filename))
ds_forc = xr.merge([ds_atm, ds_ocn], compat='override')
# Create temp above frz
ds_forc['sst_above_frz'] = ds_forc['tos'] - ds_forc['tosf']

# mapping from Icepack names to forcing names
forc_var_map = {'flw': ['rld'],
                'fsw': ['rsd'],
                'Tair': ['tas'],
                'Qa': ['hus'],
                'flwout': ['rlu'],
                #'fsens': ['wthv', 'hfss_ec', 'hfss_bulk'],
                #'flat': ['wqv', 'hfls_bulk'],
                'sst': ['tos'],
                'sss': ['so'],
                'Tf': ['tosf'],
                'sst_above_frz': ['sst_above_frz'],
                }

# Load ice properties
data_path = os.path.join(Path.home(), "data", "mass_balance_data")
data_filename = "ablationStakes_hotwireThicknessGauges_MOSAiC.csv"
df_stak = pd.read_csv(os.path.join(data_path, data_filename), parse_dates=[3,4,7])

# We want to create a dataframe where the row index is site and date and the
# column index is stat_name, and then mean and standard error of mean.
# rename columns to make our life easier
df_renamed = df_stak.rename(columns={'Site name': 'site', 
                                     'Measurement date': 'time',
                                     'Ice thickness (calculated) (cm)': 'hi_hotwire',
                                     'Drilled ice thickness (cm)': 'hi_drill',
                                     'Snow depth (calculated) (cm)': 'hs_gauge',
                                     'Pond depth (cm)': 'hpnd',
                                     'Pond flag': 'apnd'})
# cm to m
for var_name in ['hi_hotwire', 'hi_drill', 'hs_gauge', 'hpnd']:
    df_renamed[var_name] = df_renamed[var_name]/100
stat_names = ['hi_hotwire','hi_drill','hs_gauge','hpnd','apnd']
df_ice = df_renamed.groupby(['site', 'time'])[stat_names].agg(
    ['mean', 'sem'])

ice_var_map = {'vice': ['hi_hotwire', 'hi_drill'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }

# Validate FYI
run_plot_dict = {"mos_v2_kazr": [1],
                 "mos_v2_init": [1]}
var_names = ['vice', 'vsno', 'apnd', 'hpnd', 'sst', 'sst_above_frz']
site_names = ['Ridge Ranch/dart_stakes_clu_6',
              'Stakes 1/dart_stakes_clu_4',
              'Runaway Stakes/dart_stakes_clu_7',
              'Drone Bones/dart_stakes_clu_11',
              'Reunion Stakes/dart_stakes_clu_12',
              ]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice, forc_var_map=forc_var_map,
                 ds_forc=ds_forc)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2019-11-16'),
                  datetime.datetime.fromisoformat('2020-07-31')])
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-04-16T18:00:00'),
#                  datetime.datetime.fromisoformat('2020-04-17')])

axs[0].set_ylabel('Ice thickness (m)')
axs[1].set_ylabel('Snow thickness (m)')
axs[2].set_ylabel('Pond fraction')
axs[3].set_ylabel('Pond depth (m)')
axs[4].set_ylabel('Mixed layer temp. (C)')
axs[5].set_ylabel('Freezing point diff. (C)')
plt.show()

# Validate SYI (note, not exact match due to snow)
run_plot_dict = {"mos_v2_kazr": [2, 3],
                 "mos_v2_init": [2, 3]}
var_names = ['vice', 'vsno', 'apnd', 'hpnd', 'sst', 'sst_above_frz']
site_names = ['Bow Stakes/dart_stakes_clu_1',
              'Stakes 3/dart_stakes_clu_3',
              'MET Stakes/dart_stakes_clu_5',
              'Beanpole Stakes/dart_stakes_clu_13'
              ]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice, forc_var_map=forc_var_map,
                 ds_forc=ds_forc)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2019-11-16'),
                  datetime.datetime.fromisoformat('2020-07-31')])
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-04-16T18:00:00'),
#                  datetime.datetime.fromisoformat('2020-04-17')])

axs[0].set_ylabel('Ice thickness (m)')
axs[1].set_ylabel('Snow thickness (m)')
axs[2].set_ylabel('Pond fraction')
axs[3].set_ylabel('Pond depth (m)')
axs[4].set_ylabel('Mixed layer temp. (C)')
axs[5].set_ylabel('Freezing point diff. (C)')
plt.show()
