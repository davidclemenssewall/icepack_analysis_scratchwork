# examine mosaic forcing icepack output

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import json

import icepacktools as ipt

# Load model history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mosaic_forc_clean": None,
            "mosaic_forc_clean_debug": None,
            }
trcr_dict = {19: 'apnd',
             20: 'hpnd',
             21: 'ipnd'}
trcrn_dict = {19: 'apndn',
              20: 'hpndn',
              21: 'ipndn'}

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = ipt.load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value, 
                                       trcr_dict=trcr_dict, 
                                       trcrn_dict=trcrn_dict)
    
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

ice_var_map = {'vice': ['hi_hotwire', 'hi_drill'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }

# First use standard tools to check output
# Plot FYI stakes
run_plot_dict = {"mosaic_forc_clean": [1],
                 "mosaic_forc_clean_debug": [1],
                 }
var_names = ['vice', 'vsno', 'apnd', 'hpnd', 'ipnd',
             'congel', 'frazil', 'snoice', 'meltt', 'meltl', 'meltb']
ice_sites = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
              'Drone Bones/dart_stakes_clu_11',
              'Reunion Stakes/dart_stakes_clu_12',
            ]
f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, 
                          ice_var_map=ice_var_map, ice_sites=ice_sites, 
                          df_ice=df_ice)

# Plot SYI model and stakes
run_plot_dict = {"mosaic_forc_clean": [2],
                 "mosaic_forc_clean_debug": [2],
                 }
var_names = ['vice', 'vsno', 'apnd', 'hpnd', 'ipnd',
             'congel', 'frazil', 'snoice', 'meltt', 'meltl', 'meltb']
ice_sites = [ "Stakes 3/dart_stakes_clu_3",
            'MET Stakes/dart_stakes_clu_5',
            'Beanpole Stakes/dart_stakes_clu_13'
            ]
f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, 
                          ice_var_map=ice_var_map, ice_sites=ice_sites, 
                          df_ice=df_ice)

