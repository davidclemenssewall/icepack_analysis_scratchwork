# lightning_talk_20250513.py

import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import numpy as np
import datetime

from icepacktools import load_icepack_hist, plot_hist_var, plot_forc_var
from icepacktools import plot_ice_var, plot_handler, plot_freshwater_budget

import matplotlib.units as munits

converter = mdates.ConciseDateConverter()
munits.registry[np.datetime64] = converter
munits.registry[datetime.date] = converter
munits.registry[datetime.datetime] = converter

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mos_leg4_lvlpnd": 'icepack.h.20200701.nc',
            "mos_leg4_seapnd": 'icepack.h.20200701.nc',
            "mos_fyi_lvlpnd": 'icepack.h.20191129.nc',
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
                 "mos_leg4_seapnd": [3],
                 "mos_fyi_lvlpnd": [1],
                 "mos_fyi_seapnd": [1],
                 }
var_names = ['apnd', 'hpnd']
ice_sites = ['transect']
xlim = [datetime.datetime.fromisoformat('2020-06-22'),
        datetime.datetime.fromisoformat('2020-07-31')]
f, axs = plot_handler(run_plot_dict, var_names, hist_dict,
                 forc_var_map=forc_var_map,
                 ice_var_map=ice_var_map,
                 ds_forc=ds_forc, df_ice=df_ice, xlim=xlim,
                 ice_sites=ice_sites, mean_only=False)
plt.show()

# Create figures for presentation
xfont = 14
lfont = 11

f, axs = plt.subplots(2, 1, sharex=True, figsize=(6,6))

# Area fraction
handles = []
handles.append(plot_ice_var(df_ice, 'apnd', 'transect', axs[0],
                            mean_only=True, linestyle='-',
                            color='k'))
handles.append(plot_hist_var(hist_dict["mos_leg4_lvlpnd"],
                             'apnd', 3, axs[0],
                             linestyle='--', color='b'))
handles.append(plot_hist_var(hist_dict["mos_fyi_lvlpnd"],
                             'apnd', 1, axs[0],
                             linestyle='--', color='r'))
handles.append(plot_hist_var(hist_dict["mos_leg4_seapnd"],
                             'apnd', 3, axs[0],
                             linestyle='-', color='b'))
handles.append(plot_hist_var(hist_dict["mos_fyi_seapnd"],
                             'apnd', 1, axs[0],
                             linestyle='-', color='r'))

axs[0].set_ylabel('Pond area fraction', fontsize=xfont)
axs[0].grid(lw=0.5)

# Depth
handles = []
handles.append(plot_ice_var(df_ice, 'hpnd', 'transect', axs[1],
                            mean_only=True, linestyle='-',
                            color='k'))
handles.append(plot_hist_var(hist_dict["mos_leg4_lvlpnd"],
                             'hpnd', 3, axs[1],
                             linestyle='--', color='b'))
handles.append(plot_hist_var(hist_dict["mos_fyi_lvlpnd"],
                             'hpnd', 1, axs[1],
                             linestyle='--', color='r'))
handles.append(plot_hist_var(hist_dict["mos_leg4_seapnd"],
                             'hpnd', 3, axs[1],
                             linestyle='-', color='b'))
handles.append(plot_hist_var(hist_dict["mos_fyi_seapnd"],
                             'hpnd', 1, axs[1],
                             linestyle='-', color='r'))

axs[1].set_ylabel('Pond depth (m)', fontsize=xfont)
axs[1].grid(lw=0.5)
axs[1].legend(['Obs. (Transect)', 'Lvlpnd (Transect)',
            'Lvlpnd (FYI)', 'Sealvl (Transect)',
            'Sealvl (FYI)'], fontsize=lfont, ncol=2,
            loc='lower left')

xlim = [datetime.datetime.fromisoformat('2020-06-22'),
        datetime.datetime.fromisoformat('2020-07-31')]
axs[1].set_xlim(xlim)
axs[1].set_ylim([0, 0.11])
f.savefig('./figures/lightning_sealevel_20240513.png',
          bbox_inches='tight')

f, axs = plt.subplots(2, 1, sharex=True, figsize=(6,6))

# Area fraction
handles = []
handles.append(plot_ice_var(df_ice, 'apnd', 'transect', axs[0],
                            mean_only=True, linestyle='-',
                            color='k'))
handles.append(plot_hist_var(hist_dict["mos_leg4_lvlpnd"],
                             'apnd', 3, axs[0],
                             linestyle='--', color='b'))
handles.append(plot_hist_var(hist_dict["mos_fyi_lvlpnd"],
                             'apnd', 1, axs[0],
                             linestyle='--', color='r'))

axs[0].set_ylabel('Pond area fraction', fontsize=xfont)
axs[0].grid(lw=0.5)

# Depth
handles = []
handles.append(plot_ice_var(df_ice, 'hpnd', 'transect', axs[1],
                            mean_only=True, linestyle='-',
                            color='k'))
handles.append(plot_hist_var(hist_dict["mos_leg4_lvlpnd"],
                             'hpnd', 3, axs[1],
                             linestyle='--', color='b'))
handles.append(plot_hist_var(hist_dict["mos_fyi_lvlpnd"],
                             'hpnd', 1, axs[1],
                             linestyle='--', color='r'))

axs[1].set_ylabel('Pond depth (m)', fontsize=xfont)
axs[1].grid(lw=0.5)
axs[1].legend(['Obs. (Transect)', 'Lvlpnd (Transect)',
            'Lvlpnd (FYI)'], fontsize=lfont, ncol=1,
            loc='lower left')

xlim = [datetime.datetime.fromisoformat('2020-06-22'),
        datetime.datetime.fromisoformat('2020-07-31')]
axs[1].set_xlim(xlim)
axs[1].set_ylim([0, 0.11])

f.savefig('./figures/lightning_level_20240513.png',
          bbox_inches='tight')