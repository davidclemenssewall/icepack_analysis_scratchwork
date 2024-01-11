"""
Tools for comparing Icepack runs with each other and observations

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
run_dict = {"mosaic_raphael_syi": 'icepack.h.20191129.nc',
            "mosaic_raphael_fyi": 'icepack.h.20191129.nc',
            "sheba_raphael_fyi": None,
            "sheba_raphael_myi": None,
            "mosaic_raphael_L2": None,
            "mosaic_raphael_L2_pndaspect": None,
            "sheba_raphael_myi_nov28": None,
            "mosaic_raphael_fyi_rfracaice": None,
            "mosaic_raphael_fyi_nopndfbd": None,
            "mosaic_raphael_fyi_nopndfbd_pndaspect4": None,
            "mosaic_raphael_fyi_topo": None,
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

# Explore plotting variables
run_plot_dict = {"mosaic_raphael_syi": [2],
                 "mosaic_raphael_fyi": [2],
                 "sheba_raphael_fyi": [2],
                 "sheba_raphael_myi": [2],
                 "mosaic_raphael_L2": [2],
                 "sheba_raphael_myi_nov28": [2]
                }
var_names = ['aice', 'vice', 'vice_chg', 'vsno']

f = plot_handler(run_plot_dict, var_names, hist_dict)

var_names = ['apnd', 'hpnd', 'ipnd']

f = plot_handler(run_plot_dict, var_names, hist_dict)

# Test fluxes
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['vice', 'vice_chg']
f = plot_handler(run_plot_dict, var_names, hist_dict)
var_names = ['congel', 'frazil', 'snoice', 'meltt', 'meltl', 'meltb']
f = plot_handler(run_plot_dict, var_names, hist_dict, 
                 cumulative=True, resample='D', mult=24)



# Compare with changing pond parameterizations
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "mosaic_raphael_fyi_nopndfbd": [2],
                 "mosaic_raphael_fyi_nopndfbd_pndaspect4": [2],
                 "mosaic_raphael_fyi_rfracaice": [2],
                 "mosaic_raphael_fyi_topo": [2],
                 }
var_names = ['vice', 'apnd', 'hpnd', 'ipnd']
ice_var_map = {'vice': ['hi_hotwire'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }
site_names = ['Reunion Stakes/dart_stakes_clu_12']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice, mean_only=True)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-07-28')])
axs[0].set_ylabel('Ice Thickness (m)')
axs[1].set_ylabel('Pond Area Fraction')
axs[2].set_ylabel('Pond Depth (m)')
plt.show()

# Look at per category variables
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "mosaic_raphael_fyi_topo": [2],
                 }
var_names = ['vicen','aicen','apndn','hpndn','ipndn']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-07-28')])
plt.show()

# Look at per category variables
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "mosaic_raphael_fyi_nopndfbd_pndaspect4": [2],
                 }
var_names = ['vicen','aicen','apndn','hpndn','ipndn']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-07-28')])
plt.show()


# Pond aspect
run_plot_dict = {"mosaic_raphael_L2": [2],
                 "mosaic_raphael_L2_pndaspect": [2],
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'pndaspect']

f, axs = plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_ylim([0, 1.0])
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-08-07')])
plt.show()

# Explore plotting with forcing
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['flwout', 'fsens', 'flat']
f = plot_handler(run_plot_dict, var_names, hist_dict, 
                 forc_var_map=forc_var_map, ds_forc=ds_forc)
# Explore plotting with forcing
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['sst', 'Tf', 'sst_above_frz']
f = plot_handler(run_plot_dict, var_names, hist_dict, 
                 forc_var_map=forc_var_map, ds_forc=ds_forc)

# Explore plotting fyi with ice state
run_plot_dict = {"mosaic_raphael_fyi": [2]}
var_names = ['vice', 'vsno', 'apnd', 'hpnd']
site_names = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
              'Drone Bones/dart_stakes_clu_11',
              'Reunion Stakes/dart_stakes_clu_12']
f = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice)

# Mosaic SYI with ice state
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['vice', 'vsno', 'apnd', 'hpnd']
site_names = ['Bow Stakes/dart_stakes_clu_1',
              'Stakes 3/dart_stakes_clu_3',
              'MET Stakes/dart_stakes_clu_5',
              'Beanpole Stakes/dart_stakes_clu_13']
f = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice)


f, ax = plt.subplots(1,1)
plot_ice_var(df_ice, 'hpnd', 'Ridge Ranch/dart_stakes_clu_6', ax)

for k, value in hist_dict.items():
    print(k)
    print(value.sel(ni=2, time="2020-05-12")['vice_chg'].values[0])

# Explore plotting variables
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "sheba_raphael_myi": [2],
                 }
var_names = ['aice', 'vice', 'apnd', 'hpnd', 'apndn', 'hpndn']

f = plot_handler(run_plot_dict, var_names, hist_dict)
