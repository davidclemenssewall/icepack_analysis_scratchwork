"""
Tools for comparing Icepack runs with each other and observations

"""

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime

from icepacktools import load_icepack_hist, plot_hist_var, plot_forc_var
from icepacktools import plot_ice_var, plot_handler, plot_freshwater_budget

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mos_leg4_lvlpnd": 'icepack.h.20200701.nc',
            "mos_leg4_fixpnd": 'icepack.h.20200701.nc',
            "mos_leg4_seapnd": 'icepack.h.20200701.nc',
            "mos_syi_lvlpnd": 'icepack.h.20191129.nc',
            "mos_syi_fixpnd": 'icepack.h.20191129.nc',
            "mos_syi_seapnd": 'icepack.h.20191129.nc',
            "mos_fyi_lvlpnd": 'icepack.h.20191129.nc',
            "mos_fyi_fixpnd": 'icepack.h.20191129.nc',
            "mos_fyi_seapnd": 'icepack.h.20191129.nc',
            }

trcr_dict = {17: 'alvl',
             18: 'vlvl',}
trcrn_dict = {17: 'alvln',
              18: 'vlvln',}

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value, pnd_budget=True,
                                       volp=True, trcr_dict=trcr_dict,
                                       trcrn_dict=trcrn_dict)

# Create ice thickness change
for value in hist_dict.values():
    value['vice_chg'] = value['vice'] - value['vice'].isel(time=0)
    value['pndaspect'] = value['hpnd'] / value['apnd']
    value['hin'] = value['vicen'] / value['aicen']

# Load forcing
forcing_path = os.path.join(ip_dirs_path, "input", "Icepack_data",
                            "forcing", "MOSAiC")
atm_filename = "MOSAiC_MODF_20191011-20201001_v0.2.nc"
atm_filename = "MOSAiC_kazr_snow_MDF_20191005-20201001.nc"
ds_atm = xr.open_dataset(os.path.join(forcing_path, atm_filename))
ocn_filename = 'MOSAiC_ocn_MDF_20191006-20200919.nc'
ds_ocn = xr.open_dataset(os.path.join(forcing_path, ocn_filename))
ds_forc = xr.merge([ds_atm, ds_ocn])
# Create temp above frz
ds_forc['sst_above_frz'] = ds_forc['tos'] - ds_forc['tosf']
# mapping from Icepack names to forcing names
forc_var_map = {'flw': ['rlds'],
                'fsw': ['rsds'],
                'Tair': ['tas'],
                'Qa': ['hus'],
                'sst': ['tos'],
                'sss': ['so'],
                'Tf': ['tosf'],
                'sst_above_frz': ['sst_above_frz'],
                }

# Load ice properties
# Load stakes
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
# Load transect
data_path = "/home/dcsewall/data/transect_data_melinda"
data_filename = "leg4transect_dcs.csv"
df_tran = pd.read_csv(os.path.join(data_path, data_filename), parse_dates=[0])
df_tran.fillna(value=0, inplace=True)
df_tran['site'] = 'transect'
df_tran = df_tran.rename(columns={'Date time': 'time',
                                  'ice thickness (m)': 'hin_gems',
                                  'snow/SSL thickness (m)': 'hsn_magnaprobe',
                                  'pond fraction': 'apndn',
                                  'pond depth (m)': 'hpndn',
                                  'ice thickness category': 'category'})
df_tran.set_index(['site', 'time', 'category'], inplace=True)
df_agg = df_tran[['hin_gems', 'hsn_magnaprobe', 'apndn', 'hpndn']].mul(
    df_tran['aicen'], axis=0).groupby(level=[0, 1]).agg('sum')
df_agg = df_agg.rename(columns={'hin_gems': 'hi_gems',
                                'hsn_magnaprobe': 'hs_magnaprobe',
                                'apndn': 'apnd',
                                'hpndn': 'hpnd'})
df_sem = df_agg.copy()
df_sem.values[:,:] = np.NaN
#df_agg = pd.concat([df_agg], axis=1, keys=['mean']).swaplevel(0, 1, 1)
#df_sem = pd.concat([df_sem], axis=1, keys=['sem']).swaplevel(0, 1, 1)
df_agg = pd.concat([df_agg, df_sem], axis=1, keys=['mean', 'sem']
                   ).reorder_levels([1, 0], axis=1).sort_index(axis=1)
df_ice = pd.concat([df_ice, df_agg], axis=0)

ice_var_map = {'vice': ['hi_hotwire', 'hi_drill', 'hi_gems'],
               'vsno': ['hs_gauge', 'hs_magnaprobe'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               'hin': ['hin_gems'],
               'hpndn': ['hpndn'],
               'apndn': ['apndn'],
               'aicen': ['aicen']
               }

# Validate Simulation
run_plot_dict = {"mos_leg4_lvlpnd": [3],
                 "mos_leg4_fixpnd": [3],
                 "mos_leg4_seapnd": [3],
                 "mos_syi_lvlpnd": [2],
                 "mos_syi_fixpnd": [2],
                 "mos_syi_seapnd": [2],
                 "mos_fyi_lvlpnd": [1],
                 "mos_fyi_fixpnd": [1],
                 "mos_fyi_seapnd": [1],}
var_names = ['aice', 'vice', 'vsno', 
             'Tair', 'Qa', 'fsw', 'flw', 
             'sst', 'Tf', 'sst_above_frz', 'sss',
             'apnd', 'hpnd', 'ipnd',
             'congel', 'frazil', 'snoice', 'meltt', 'meltb', 'meltl']
ice_sites = ['Ridge Ranch/dart_stakes_clu_6',
             'Drone Bones/dart_stakes_clu_11',
             'Reunion Stakes/dart_stakes_clu_12',
             'Beanpole Stakes/dart_stakes_clu_13',
             'transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, xlim=xlim,
                 ice_sites=ice_sites, mean_only=False)
axs[1].set_ylabel('Ice thickness (m)')
axs[3].set_ylabel('Air temperature (K)')
axs[3].set_ylim([268, 275])
axs[4].set_ylabel('Specific humidity')
axs[4].set_ylim([0.0028, 0.0045])
axs[5].set_ylabel('Shortwave down (W/m2)')
axs[6].set_ylabel('Longwave down (W/m2)')
axs[6].set_ylim([240, 350])
plt.show()

# Examine category variables
run_plot_dict = {"mos_leg4_lvlpnd": [3]}
var_names = ['aicen', 'hin', 'apndn', 'hpndn', 'ipndn']
ice_sites = ['transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, df_icen=df_tran, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
f.suptitle('Level pond')

# Examine category variables
run_plot_dict = {"mos_leg4_fixpnd": [3]}
var_names = ['aicen', 'hin', 'apndn', 'hpndn', 'ipndn']
ice_sites = ['transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, df_icen=df_tran, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
f.suptitle('Fixed aspect')

# Examine category variables
run_plot_dict = {"mos_leg4_seapnd": [3]}
var_names = ['aicen', 'hin', 'apndn', 'hpndn', 'ipndn']
ice_sites = ['transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, df_icen=df_tran, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
f.suptitle('Sealevel')

# Below here is sequence of figures to show in meeting
# Level pond, leg4
run_plot_dict = {"mos_leg4_lvlpnd": [3],
                 }
var_names = ['vice', 
             'apnd', 'hpnd', 'ipnd',
             ]
ice_sites = ['Ridge Ranch/dart_stakes_clu_6',
             'Drone Bones/dart_stakes_clu_11',
             'Reunion Stakes/dart_stakes_clu_12',
             'transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
plt.show()

# Level pond, leg4 and fyi
run_plot_dict = {"mos_leg4_lvlpnd": [3],
                 "mos_fyi_lvlpnd": [1],
                 }
var_names = ['vice', 
             'apnd', 'hpnd',
             ]
ice_sites = ['Ridge Ranch/dart_stakes_clu_6',
             'Drone Bones/dart_stakes_clu_11',
             'Reunion Stakes/dart_stakes_clu_12',
             'transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
plt.show()

# Level pond fixpond, leg4 and fyi
run_plot_dict = {"mos_leg4_lvlpnd": [3],
                 "mos_fyi_lvlpnd": [1],
                 "mos_leg4_fixpnd": [3],
                 "mos_fyi_fixpnd": [1],
                 }
var_names = ['vice', 
             'apnd', 'hpnd',
             ]
ice_sites = ['Ridge Ranch/dart_stakes_clu_6',
             'Drone Bones/dart_stakes_clu_11',
             'Reunion Stakes/dart_stakes_clu_12',
             'transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
plt.show()

# Examine category variables
run_plot_dict = {"mos_leg4_fixpnd": [3]}
var_names = ['aicen', 'hin', 'apndn', 'hpndn']
ice_sites = ['transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, df_icen=df_tran, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
f.suptitle('Fixed aspect')

# Level pond fixpond, leg4 and fyi
run_plot_dict = {"mos_leg4_lvlpnd": [3],
                 "mos_fyi_lvlpnd": [1],
                 "mos_leg4_fixpnd": [3],
                 "mos_fyi_fixpnd": [1],
                 "mos_leg4_seapnd": [3],
                 "mos_fyi_seapnd": [1],
                 }
var_names = ['vice', 
             'apnd', 'hpnd',
             ]
ice_sites = ['Ridge Ranch/dart_stakes_clu_6',
             'Drone Bones/dart_stakes_clu_11',
             'Reunion Stakes/dart_stakes_clu_12',
             'transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
plt.show()

# Examine category variables
run_plot_dict = {"mos_leg4_seapnd": [3]}
var_names = ['aicen', 'hin', 'apndn', 'hpndn']
ice_sites = ['transect']
xlim = [datetime.datetime.fromisoformat('2020-06-15'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, df_icen=df_tran, xlim=xlim,
                 ice_sites=ice_sites, mean_only=True)
f.suptitle('Sealevel')

# Examine freshwater budget
# Plot pond budget
f, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 10))

plot_freshwater_budget(hist_dict['mos_leg4_lvlpnd'], 3, axs[0])
axs[0].legend(loc='lower left')
axs[0].set_ylabel("meltwater (m), lvlpnd")

plot_freshwater_budget(hist_dict['mos_leg4_seapnd'], 3, axs[1])
axs[1].legend(loc='lower left')
axs[1].set_ylabel("meltwater (m), sealvl")

xlim = [datetime.datetime.fromisoformat('2020-06-26'),
                    datetime.datetime.fromisoformat('2020-07-31')]
axs[-1].set_xlim(xlim)
plt.show()