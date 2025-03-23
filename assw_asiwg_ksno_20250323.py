# figures asiwg talk

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import json

import icepacktools as ipt

# Load model history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mos_raph_fyi_ks030": None,
            "mos_raph_fyi_ks025": None,
            "mos_raph_fyi_ks042": None,
            "mos_raph_syi_ks030": None,
            "mos_raph_syi_ks025": None,
            "mos_raph_syi_ks042": None,
            "mos_buoy_ks030": None,
            "mos_buoy_ks025": None,
            "mos_buoy_ks042": None,
            }

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = ipt.load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value)
    
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
for col in df_ice.columns:
    if col[1] == 'sem':
        df_ice[(col[0], 'semx2')] = df_ice[col]*2

# Load buoy properties
# Load paths
with open(os.path.join('/','home', 'dcsewall', 'code', 'SCMStandz', 'utils', 
                       'paths.json')) as fp:
    paths_dict = json.load(fp)

# Load buoy timeseries
imb_fp = os.path.join(paths_dict['raw_data'], 'imb')
ds_snow = xr.open_dataset(os.path.join(imb_fp, 
                        'MOSAiC_winter_IMB_snow_depths_v2.nc'))
ds_hi = xr.open_dataset(os.path.join(imb_fp, 
                        'MOSAiC_winter_IMB_SIT_v2.nc'))
# Convert to m
ds_snow = ds_snow / 100.0
ds_hi = ds_hi / 100.0
# List of buoys for comparison
subset_buoys = ['SIMBA_2019T58',
                'SIMBA_2019T70',
                'SIMBA_2019T62',
                'SIMBA_2019T72',
                'SIMB3_L2']

ice_var_map = {'vice': ['hi_hotwire', 'hi_drill'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }

# First use standard tools to check output
# Plot FYI stakes
run_plot_dict = {"mos_raph_fyi_ks030": [1],
                 "mos_raph_fyi_ks025": [1],
                 "mos_raph_fyi_ks042": [1],}
var_names = ['vice', 'vsno']
ice_sites = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
            ]
f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, 
                          ice_var_map=ice_var_map, ice_sites=ice_sites, 
                          df_ice=df_ice)

# Plot SYI model and all of the stakes
run_plot_dict = {"mos_raph_syi_ks030": [1,2],
                 "mos_raph_syi_ks025": [1,2],
                 "mos_raph_syi_ks042": [1,2],}
var_names = ['vice', 'vsno']
ice_sites = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
              "Stakes 3/dart_stakes_clu_3",
            'MET Stakes/dart_stakes_clu_5'
            ]
f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, 
                          ice_var_map=ice_var_map, ice_sites=ice_sites, 
                          df_ice=df_ice)

# Plot just stakes
# Plot SYI model and all of the stakes
run_plot_dict = {}
var_names = ['vice', 'vsno']
ice_sites = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
              "Stakes 3/dart_stakes_clu_3",
            'MET Stakes/dart_stakes_clu_5'
            ]
f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, 
                          ice_var_map=ice_var_map, ice_sites=ice_sites, 
                          df_ice=df_ice)

# Plot buoys
run_plot_dict = {"mos_buoy_ks030": [2],
                 "mos_buoy_ks025": [2],
                 "mos_buoy_ks042": [2],
                 }
var_names = ['vice', 'vsno']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, 
                          ice_var_map=ice_var_map,
                          df_ice=df_ice)
for name in subset_buoys:
    ds_hi[name].plot.line(alpha=0.6, ax=axs[0], label=name)
    ds_snow[name].plot.line(alpha=0.6, ax=axs[1], label=name)
axs[0].legend()
axs[1].legend()
axs[0].set_ylabel('vice')
axs[1].set_ylabel('vsno')

f.show()

############################################################################

# Plot ice thickness data only
f, ax = plt.subplots(1, 1, sharex=True, figsize=(6,6))
lfont = 10
axfont = 12

# Plot stakes
var_name = 'vice'
ice_sites = ['Bow Stakes/dart_stakes_clu_1',
       'MET Stakes/dart_stakes_clu_5',
       'Miss Stakes/dart_stakes_clu_8',
       'Ridge Ranch/dart_stakes_clu_6',
       'Runaway Stakes/dart_stakes_clu_7', 'Stakes 1/dart_stakes_clu_4',
       'Stakes 3/dart_stakes_clu_3',
       ]
# Plot ice variables
for name in subset_buoys:
    ds_hi[name].plot.line(linestyle='-.', ax=ax)

for site in ice_sites:
    if var_name in ice_var_map:
        for ice_var_name in ice_var_map[var_name]:
            _ = ipt.plot_ice_var(df_ice, ice_var_name, site, ax, 
                                mean_only=False)

ax.grid()
ax.set_ylabel('Ice thickness (m)', fontsize=axfont)
ax.set_xlim([datetime.datetime.fromisoformat('2019-11-16'),
                  datetime.datetime.fromisoformat('2020-05-15')])
ax.set_xlabel('')
plt.show()

#############################################################################
lfont = 10
axfont = 12
lmodel = ['Icepack, ks=0.3', 'Icepack, ks=0.25', 'Icepack, ks=0.43']

# FYI
f, ax = plt.subplots(1, 1, sharex=True, figsize=(6,6))
lfont = 10
axfont = 12

run_plot_dict = {"mos_raph_fyi_ks030": [1],
                 "mos_raph_fyi_ks025": [1],
                 "mos_raph_fyi_ks042": [1],}
# Plot stakes
var_name = 'vice'
ice_sites = [
       'Ridge Ranch/dart_stakes_clu_6',
       'Runaway Stakes/dart_stakes_clu_7', 'Stakes 1/dart_stakes_clu_4',
       ]
# Plot ice variables
for run_name, nis in run_plot_dict.items():
    # and the desired cell(s) in each run
    for ni in nis:
        _ = ipt.plot_hist_var(hist_dict[run_name], var_name, ni, ax)

for site in ice_sites:
    if var_name in ice_var_map:
        for ice_var_name in ice_var_map[var_name]:
            _ = ipt.plot_ice_var(df_ice, ice_var_name, site, ax, 
                                mean_only=False)

ax.grid()
ax.set_ylabel('First Year Ice thickness (m)', fontsize=axfont)
ax.set_xlim([datetime.datetime.fromisoformat('2019-11-16'),
                  datetime.datetime.fromisoformat('2020-05-15')])
ax.legend(lmodel, fontsize=lfont)
ax.set_xlabel('')
plt.show()