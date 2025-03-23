# icepack_init_options_validate.py

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime

import icepacktools as ipt

# Load history output
#ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
ip_dirs_path = "/root/icepack-dirs"
run_dict = {"base": None,
            "test_init_options": None,
            }

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = ipt.load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value)

# The following is useful to check that outputs match
if True:
    # Check whether or not dataarrays are identical
    comp = "test_init_options"
    for ni in [1, 2, 3, 4]:
        print(ni)
        for key, da in hist_dict["base"].sel(ni=ni).data_vars.items():
            if not da.equals(hist_dict[comp].sel(ni=ni)[key]):
                print(key)
                print('max diff: ' + str((hist_dict[comp].sel(ni=ni)[key] - da).max().values))
                print('min diff: ' + str((hist_dict[comp].sel(ni=ni)[key] - da).min().values))
    print("Above are data arrays that do not match.")

# Explore state variables
run_plot_dict = {"base": [1, 2, 3],
                 "test_init_options": [1, 2, 3],
                }
var_names = ['aice', 'vice', 'vsno', 'sst', 'sss']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-01-01'),
#                  datetime.datetime.fromisoformat('2015-01-15')])
plt.show()