# clivar_talk_20240117.py

import xarray as xr
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.dates as mdates
import pandas as pd
import numpy as np
import datetime

import icepacktools as ipt

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"clivar_mosfyi": None,
            "clivar_mossyi": None,
            "clivar_shemyi": None,
            "clivar_chamyi": None,
            "clivar_mosfyi_qdp": None,
            "clivar_mossyi_qdp": None,
            "clivar_shemyi_qdp": None,
            "clivar_chamyi_qdp": None,
            "clivar_mosfyi_ksno": None,
            "clivar_mossyi_ksno": None,
            "clivar_shemyi_ksno": None,
            "clivar_chamyi_ksno": None,
            }

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = ipt.load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value, snhf=True)
    
# Pond aspect
run_plot_dict = {"clivar_mosfyi": [1, 2],
                 "clivar_mossyi": [2],
                 "clivar_shemyi": [2],
                 "clivar_chamyi": [2],
                 "clivar_mosfyi_qdp": [1, 2],
                 "clivar_mossyi_qdp": [2],
                 "clivar_shemyi_qdp": [2],
                 "clivar_chamyi_qdp": [2],
                }
var_names = ['aice', 'vice', 'vsno']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
#                  datetime.datetime.fromisoformat('2020-08-07')])
plt.show()

# Pond aspect
run_plot_dict = {"clivar_mosfyi": [1, 2],
                 "clivar_mossyi": [2],
                 "clivar_shemyi": [2],
                 "clivar_chamyi": [2],
                 "clivar_mosfyi_ksno": [1, 2],
                 "clivar_mossyi_ksno": [2],
                 "clivar_shemyi_ksno": [2],
                 "clivar_chamyi_ksno": [2],
                }
var_names = ['aice', 'vice', 'vsno']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
#                  datetime.datetime.fromisoformat('2020-08-07')])
plt.show()

print("ice thickness mos ow on may 1")
print(hist_dict['clivar_mosfyi'].sel(ni=1, time="2020-05-01").vice.mean().values)
print(hist_dict['clivar_mosfyi_ksno'].sel(ni=1, time="2020-05-01").vice.mean().values)

print("ice thickness cha on may 1")
print(hist_dict['clivar_chamyi'].sel(ni=2, time="2020-05-01").vice.mean().values)
print(hist_dict['clivar_chamyi_ksno'].sel(ni=2, time="2020-05-01").vice.mean().values)

# Color plots
color_plot_dict = {" Open Water": ("clivar_mosfyi", 1),
                   " MOSAiC FYI": ("clivar_mosfyi", 2),
                   " MOSAiC SYI": ("clivar_mossyi", 2),
                   "  SHEBA MYI": ("clivar_shemyi", 2),
                   "CHARLIE MYI": ("clivar_chamyi", 2),
                   }
dt_min = "2019-12-01"
dt_max = "2020-04-30"
resample = 'M'
hi_dict = {}
hi_max = -np.inf
snhf_dict = {}
snhf_min = np.inf
snhf_max = -np.inf

for key, val in color_plot_dict.items():
    df = hist_dict[val[0]].sel(ni=val[1], time=slice(dt_min, dt_max))[
        ['vice', 'snhf']].to_pandas().resample(resample).mean().reset_index()
    td = df.at[1, 'time'] - df.at[0, 'time']
    hi_max = max([hi_max, df.vice.max()])
    hi_dict[key] = df.vice.values
    snhf_min = min([snhf_min, df.snhf.min()])
    snhf_max = max([snhf_max, df.snhf.max()])
    snhf_dict[key] = df.snhf.values

f, axs = plt.subplots(1, 2, sharey=True, figsize=(10, 7), tight_layout=True)

for key, hi in hi_dict.items():
    axs[0].barh([key], width=td, color=cm.cividis(hi/hi_max), left=df.time-td)
for key, snhf in snhf_dict.items():
    axs[1].barh([key], width=td, color=cm.plasma_r((snhf-snhf_min)
                                                   /(snhf_max-snhf_min)), left=df.time-td)

axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
norm = Normalize(vmin=0, vmax=hi_max)
f.colorbar(cm.ScalarMappable(norm=norm, cmap='cividis'), ax=axs[0], 
           orientation='horizontal', label='Ice thickness (m)',
           location='top')
axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
norm = Normalize(vmin=snhf_min, vmax=snhf_max)
f.colorbar(cm.ScalarMappable(norm=norm, cmap='plasma_r'), ax=axs[1], 
           orientation='horizontal', label='Net air-surface heat flux (W/m2)',
           location='top')

# Difference plots
# qdp
diff = 'qdp'
color_plot_dict = {" Open Water": ("clivar_mosfyi", 1),
                   " MOSAiC FYI": ("clivar_mosfyi", 2),
                   " MOSAiC SYI": ("clivar_mossyi", 2),
                   "  SHEBA MYI": ("clivar_shemyi", 2),
                   "CHARLIE MYI": ("clivar_chamyi", 2),
                   }
dt_min = "2019-12-01"
dt_max = "2020-04-30"
resample = 'M'
hi_dict = {}
hi_min = np.inf
hi_max = -np.inf
snhf_dict = {}
snhf_min = np.inf
snhf_max = -np.inf

for key, val in color_plot_dict.items():
    df_base = hist_dict[val[0]].sel(ni=val[1], time=slice(dt_min, dt_max))[
        ['vice', 'snhf']].to_pandas().resample(resample).mean()
    df_mod = hist_dict[val[0]+'_'+diff].sel(ni=val[1], time=slice(dt_min, dt_max))[
        ['vice', 'snhf']].to_pandas().resample(resample).mean()
    df = df_mod - df_base
    df = df.reset_index()
    td = df.at[1, 'time'] - df.at[0, 'time']
    hi_min = min([hi_min, df.vice.min()])
    hi_max = max([hi_max, df.vice.max()])
    hi_dict[key] = df.vice.values
    snhf_min = min([snhf_min, df.snhf.min()])
    snhf_max = max([snhf_max, df.snhf.max()])
    snhf_dict[key] = df.snhf.values

abs_hi_max = max([abs(hi_min), abs(hi_max)])
abs_snhf_max = max([abs(snhf_min), abs(snhf_max)])

f, axs = plt.subplots(1, 2, sharey=True, figsize=(10, 7), tight_layout=True)

for key, hi in hi_dict.items():
    axs[0].barh([key], width=td, color=cm.RdBu((hi+abs_hi_max)/
                                                  (2*abs_hi_max)), left=df.time-td)
for key, snhf in snhf_dict.items():
    axs[1].barh([key], width=td, color=cm.RdBu((snhf+abs_snhf_max)
                                                   /(2*abs_snhf_max)), left=df.time-td)

axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
norm = Normalize(vmin=-abs_hi_max, vmax=abs_hi_max)
f.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'), ax=axs[0], 
           orientation='horizontal', label='Change in ice thickness (m)',
           location='top')
axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
norm = Normalize(vmin=-abs_snhf_max, vmax=abs_snhf_max)
f.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'), ax=axs[1], 
           orientation='horizontal', label='Change in net air-surface heat flux (W/m2)',
           location='top')

# Difference plots
# ksno
diff = 'ksno'
color_plot_dict = {" Open Water": ("clivar_mosfyi", 1),
                   " MOSAiC FYI": ("clivar_mosfyi", 2),
                   " MOSAiC SYI": ("clivar_mossyi", 2),
                   "  SHEBA MYI": ("clivar_shemyi", 2),
                   "CHARLIE MYI": ("clivar_chamyi", 2),
                   }
dt_min = "2019-12-01"
dt_max = "2020-04-30"
resample = 'M'
hi_dict = {}
hi_min = np.inf
hi_max = -np.inf
snhf_dict = {}
snhf_min = np.inf
snhf_max = -np.inf

for key, val in color_plot_dict.items():
    df_base = hist_dict[val[0]].sel(ni=val[1], time=slice(dt_min, dt_max))[
        ['vice', 'snhf']].to_pandas().resample(resample).mean()
    df_mod = hist_dict[val[0]+'_'+diff].sel(ni=val[1], time=slice(dt_min, dt_max))[
        ['vice', 'snhf']].to_pandas().resample(resample).mean()
    df = df_mod - df_base
    df = df.reset_index()
    td = df.at[1, 'time'] - df.at[0, 'time']
    hi_min = min([hi_min, df.vice.min()])
    hi_max = max([hi_max, df.vice.max()])
    hi_dict[key] = df.vice.values
    snhf_min = min([snhf_min, df.snhf.min()])
    snhf_max = max([snhf_max, df.snhf.max()])
    snhf_dict[key] = df.snhf.values

abs_hi_max = max([abs(hi_min), abs(hi_max)])
abs_snhf_max = max([abs(snhf_min), abs(snhf_max)])

f, axs = plt.subplots(1, 2, sharey=True, figsize=(10, 7), tight_layout=True)

for key, hi in hi_dict.items():
    axs[0].barh([key], width=td, color=cm.RdBu((hi+abs_hi_max)/
                                                  (2*abs_hi_max)), left=df.time-td)
for key, snhf in snhf_dict.items():
    axs[1].barh([key], width=td, color=cm.RdBu((snhf+abs_snhf_max)
                                                   /(2*abs_snhf_max)), left=df.time-td)

axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
norm = Normalize(vmin=-abs_hi_max, vmax=abs_hi_max)
f.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'), ax=axs[0], 
           orientation='horizontal', label='Change in ice thickness (m)',
           location='top')
axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
norm = Normalize(vmin=-abs_snhf_max, vmax=abs_snhf_max)
f.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'), ax=axs[1], 
           orientation='horizontal', label='Changie in net air-surface heat flux (W/m2)',
           location='top')

# Difference plots any values
diff = 'qdp'
color_plot_dict = {" Open Water": ("clivar_mosfyi", 1),
                   " MOSAiC FYI": ("clivar_mosfyi", 2),
                   " MOSAiC SYI": ("clivar_mossyi", 2),
                   "  SHEBA MYI": ("clivar_shemyi", 2),
                   "CHARLIE MYI": ("clivar_chamyi", 2),
                   }
dt_min = "2019-12-01"
dt_max = "2020-04-30"
resample = 'M'
var_names = ['vice', 'flwout', 'fsens', 'flat', 'fcondtop']
var_dict = {}
min_dict = {}
max_dict = {}
for var_name in var_names:
    var_dict[var_name] = {}
    min_dict[var_name] = np.inf
    max_dict[var_name] = -np.inf

for key, val in color_plot_dict.items():
    df_base = hist_dict[val[0]].sel(ni=val[1], time=slice(dt_min, dt_max))[
        var_names].to_pandas().resample(resample).mean()
    df_mod = hist_dict[val[0]+'_'+diff].sel(ni=val[1], time=slice(dt_min, dt_max))[
        var_names].to_pandas().resample(resample).mean()
    df = df_mod - df_base
    df = df.reset_index()
    td = df.at[1, 'time'] - df.at[0, 'time']
    for var_name in var_names:
        var_dict[var_name][key] = df[var_name].values
        min_dict[var_name] = min([min_dict[var_name], df[var_name].min()])
        max_dict[var_name] = max([max_dict[var_name], df[var_name].max()])

abs_dict = {}
for var_name in var_names:
    abs_dict[var_name] = max([abs(min_dict[var_name]), abs(max_dict[var_name])])

f, axs = plt.subplots(1, len(var_names), sharey=True, figsize=(10, 7), 
                      tight_layout=True)

for ax, var_name in zip(axs, var_names):
    for key, vals in var_dict[var_name].items():
        ax.barh([key], width=td, color=cm.RdBu((vals+abs_dict[var_name])/
                                                  (2*abs_dict[var_name])), left=df.time-td)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    norm = Normalize(vmin=-abs_dict[var_name], vmax=abs_dict[var_name])
    f.colorbar(cm.ScalarMappable(norm=norm, cmap='RdBu'), ax=ax, 
           orientation='horizontal', label='Change in ' + var_name,
           location='top')
