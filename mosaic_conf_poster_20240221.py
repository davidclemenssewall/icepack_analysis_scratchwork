# figures for single column modeling poster

import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import datetime
import json

import icepacktools as ipt

# Load model history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mos_raph_fyi_ks030": None,
            "mos_raph_fyi_ks025": None,
            "mos_raph_fyi_ks042": None,
            "mos_buoy_ks030": None,
            "mos_buoy_ks025": None,
            "mos_buoy_ks042": None,
            "mos_raph_fyi_frshbud": None,
            "mos_raph_fyi_pnd_allmods": None,
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

ice_var_map = {'vice': ['hi_hotwire'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }

# Plot full year for poster
f, axs = plt.subplots(4, 1, sharex=True, figsize=(12.5,12.5))
lfont = 14
axfont = 18
lmodel = ['Icepack, ks=0.3', 'Icepack, ks=0.25', 'Icepack, ks=0.43']

# Plot FYI stakes vice
run_plot_dict = {"mos_raph_fyi_ks030": [1],
                 "mos_raph_fyi_ks025": [1],
                 "mos_raph_fyi_ks042": [1],}
var_name = 'vice'
ice_sites = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
              'Drone Bones/dart_stakes_clu_11',
              'Reunion Stakes/dart_stakes_clu_12',
            ]
for run_name, nis in run_plot_dict.items():
    # and the desired cell(s) in each run
    for ni in nis:
        _ = ipt.plot_hist_var(hist_dict[run_name], var_name, ni, axs[0])
# Plot ice variables
for site in ice_sites:
    if var_name in ice_var_map:
        for ice_var_name in ice_var_map[var_name]:
            _ = ipt.plot_ice_var(df_ice, ice_var_name, site, axs[0], 
                                mean_only=False)
axs[0].legend()
axs[0].grid()
axs[0].legend(lmodel+['Stakes 1', 'Ridge Ranch', 'Runaway Stakes', 'Drone Bones', 'Reunion Stakes'],
                      fontsize=lfont, bbox_to_anchor=(1.0, 1.0), loc='upper left')
axs[0].set_ylabel('Stakes FYI\nIce thickness (m)', fontsize=axfont)
axs[0].yaxis.set_tick_params(labelsize=lfont)

# Stakes snow depth
var_name = 'vsno'
for run_name, nis in run_plot_dict.items():
    # and the desired cell(s) in each run
    for ni in nis:
        _ = ipt.plot_hist_var(hist_dict[run_name], var_name, ni, axs[1])
# Plot ice variables
for site in ice_sites:
    if var_name in ice_var_map:
        for ice_var_name in ice_var_map[var_name]:
            _ = ipt.plot_ice_var(df_ice, ice_var_name, site, axs[1], 
                                mean_only=False)
axs[1].legend()
axs[1].grid()
axs[1].legend(lmodel+['Stakes 1', 'Ridge Ranch', 'Runaway Stakes', 'Drone Bones', 'Reunion Stakes'],
                      fontsize=lfont, bbox_to_anchor=(1.0, 1.0), loc='upper left')
axs[1].set_ylabel('Stakes FYI\nSnow depth (m)', fontsize=axfont)
axs[1].yaxis.set_tick_params(labelsize=lfont)

# Plot Buoys ice thickness
run_plot_dict = {"mos_buoy_ks030": [2],
                 "mos_buoy_ks025": [2],
                 "mos_buoy_ks042": [2],}
var_name = 'vice'
for run_name, nis in run_plot_dict.items():
    # and the desired cell(s) in each run
    for ni in nis:
        _ = ipt.plot_hist_var(hist_dict[run_name], var_name, ni, axs[2])
for name in subset_buoys:
    ds_hi[name].plot.line(linestyle='--', ax=axs[2])
    
axs[2].grid()
axs[2].legend(lmodel+subset_buoys,
              fontsize=lfont, bbox_to_anchor=(1.0, 1.0), loc='upper left')
axs[2].set_ylabel('Buoys ice\nthickness (m)', fontsize=axfont)
axs[2].yaxis.set_tick_params(labelsize=lfont)
axs[2].set_xlabel('')

# Plot Buoys snow depth
var_name = 'vsno'
for run_name, nis in run_plot_dict.items():
    # and the desired cell(s) in each run
    for ni in nis:
        _ = ipt.plot_hist_var(hist_dict[run_name], var_name, ni, axs[3])
for name in subset_buoys:
    ds_snow[name].plot.line(linestyle='--', ax=axs[3])
    
axs[3].grid()
axs[3].legend(lmodel+subset_buoys,
              fontsize=lfont, bbox_to_anchor=(1.0, 1.0), loc='upper left')
axs[3].set_ylabel('Buoys snow\ndepth (m)', fontsize=axfont)
axs[3].yaxis.set_tick_params(labelsize=lfont)
axs[3].set_xlabel('')
axs[3].xaxis.set_tick_params(labelsize=lfont)

plt.show()

f.savefig("./figures/mosaic_conf_poster_full_year.png", bbox_inches='tight')

# Compare with changing pond parameterizations
run_plot_dict = {"mos_raph_fyi_frshbud": [2],
             "mos_raph_fyi_pnd_allmods": [2],
                 }
var_names = ['vice', 'apnd', 'hpnd', 'ipnd']
ice_var_map = {'vice': ['hi_hotwire'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }
ice_sites = [ 'Ridge Ranch/dart_stakes_clu_6',
              'Drone Bones/dart_stakes_clu_11',
              'Reunion Stakes/dart_stakes_clu_12',
            ]
f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=ice_sites, df_ice=df_ice, mean_only=True,
                 figsize=(12.5,12.5), ax_font=axfont, lfont=lfont,
                 xlim=[datetime.datetime.fromisoformat('2020-05-24'),
                  datetime.datetime.fromisoformat('2020-08-01')])

axs[0].set_ylabel('Ice Thickness (m)')
axs[1].set_ylabel('Pond Area Fraction')
axs[2].set_ylabel('Pond Depth (m)')
axs[3].set_ylabel('Pond Lid Thickness (m)')
axs[0].legend(['Icepack default', 'Sea level ponds', 'Ridge Ranch', 
               'Drone Bones', 'Reunion Stakes'],
              fontsize=lfont, bbox_to_anchor=(1.0, 1.0), loc='upper left')
axs[1].legend(['Icepack default', 'Sea level ponds', 'Ridge Ranch', 
               'Drone Bones', 'Reunion Stakes'],
              fontsize=lfont, bbox_to_anchor=(1.0, 1.0), loc='upper left')
axs[2].legend(['Icepack default', 'Sea level ponds', 'Ridge Ranch', 
               'Drone Bones', 'Reunion Stakes'],
              fontsize=lfont, bbox_to_anchor=(1.0, 1.0), loc='upper left')
axs[3].legend(['Icepack default', 'Sea level ponds'],
              fontsize=lfont, bbox_to_anchor=(1.0, 1.0), loc='upper left')
axs[3].xaxis.set_tick_params(rotation=30)
plt.show()

f.savefig("./figures/mosaic_conf_poster_summer.png", bbox_inches='tight')