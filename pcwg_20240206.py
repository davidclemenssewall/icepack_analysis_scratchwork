"""
PCWG 2024 meeting

"""

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime

from icepacktools import load_icepack_hist, plot_hist_var, plot_forc_var
from icepacktools import plot_ice_var, plot_handler

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mos_raph_fyi_frshbud": None,
            "mos_raph_fyi_pnd_allmods": None,
            }

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value)
    hist_dict[key]['hin'] = hist_dict[key]['vicen']/hist_dict[key]['aicen']

# Create ice thickness change
for value in hist_dict.values():
    value['vice_chg'] = value['vice'] - value['vice'].isel(time=0)
    value['pndaspect'] = value['hpnd'] / value['apnd']

# Load forcing
forcing_path = os.path.join(ip_dirs_path, "input", "Icepack_data",
                            "forcing", "MOSAiC")
atm_filename = "MOSAiC_MODF_20191011-20201001_v0.2.nc"
ds_atm = xr.open_dataset(os.path.join(forcing_path, atm_filename))
# Surface fluxes have opposite sign in Icepack output
for var_name in ['rlu', 'wthv', 'hfss_ec', 'hfss_bulk', 'wqv', 'hfls_bulk']:
    ds_atm[var_name] = -1*ds_atm[var_name]
ds_atm['hus'] = ds_atm['hus']/1000 # fix specific humidity error
ocn_filename = 'MOSAiC_ocn_MDF_20191007-20200920.nc'
ds_ocn = xr.open_dataset(os.path.join(forcing_path, ocn_filename))
ds_forc = xr.merge([ds_atm, ds_ocn])
# Create temp above frz
ds_forc['sst_above_frz'] = ds_forc['tos'] - ds_forc['tosf']

# mapping from Icepack names to forcing names
forc_var_map = {'flw': ['rld'],
                'fsw': ['rsd'],
                'Tair': ['tas'],
                'Qa': ['hus'],
                'flwout': ['rlu'],
                'fsens': ['wthv', 'hfss_ec', 'hfss_bulk'],
                'flat': ['wqv', 'hfls_bulk'],
                'sst': ['tos'],
                'sss': ['so'],
                'Tf': ['tosf'],
                'sst_above_frz': ['sst_above_frz'],
                }

# Load ice properties
data_path = "/home/dcsewall/data/mass_balance_data"
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

# Plot Icepack default
run_plot_dict = {"mos_raph_fyi_frshbud": [2]}
var_names = ['vice', 'apnd', 'hpnd', 'ipnd']
site_names = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
              'Drone Bones/dart_stakes_clu_11',
              'Reunion Stakes/dart_stakes_clu_12']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice, forc_var_map=forc_var_map,
                 ds_forc=ds_forc, mean_only=True)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-01'),
                  datetime.datetime.fromisoformat('2020-07-28')])
axs[0].set_ylabel('Ice thickness (m)')
axs[1].set_ylabel('Pond fraction')
axs[2].set_ylabel('Pond depth (m)')
axs[3].set_ylabel('Pond lid thicknes (m)')
plt.show()

# Compare with Sea level linear ponds
run_plot_dict = {"mos_raph_fyi_frshbud": [2],
                 "mos_raph_fyi_pnd_allmods": [2]}
var_names = ['vice', 'apnd', 'hpnd', 'ipnd']
site_names = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
              'Drone Bones/dart_stakes_clu_11',
              'Reunion Stakes/dart_stakes_clu_12']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice, forc_var_map=forc_var_map,
                 ds_forc=ds_forc, mean_only=True)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-01'),
                  datetime.datetime.fromisoformat('2020-07-28')])
axs[0].set_ylabel('Ice thickness (m)')
axs[1].set_ylabel('Pond fraction')
axs[2].set_ylabel('Pond depth (m)')
axs[3].set_ylabel('Pond lid thicknes (m)')
plt.show()