import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime

import icepacktools as ipt

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"melinda_mosaic_mtg": "icepack.h.20200617.nc",
            }

trcrn_dict = {19: 'apndn',
              20: 'hpndn',
              21: 'ipndn'}


hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = ipt.load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value,
                                       trcrn_dict=trcrn_dict)
    hist_dict[key]['hin'] = hist_dict[key]['vicen']/hist_dict[key]['aicen']
    hist_dict[key]['hsn'] = hist_dict[key]['vsnon']/hist_dict[key]['aicen']

# Plot
run_plot_dict = {"melinda_mosaic_mtg": [3]
                 }
var_names = ['aicen', 'hin', 'hsn', 'hpndn', 'apndn']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
#                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

run_plot_dict = {"melinda_mosaic_mtg": [3]
                 }
var_names = ['meltt', 'meltb']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict,
                          cumulative=True)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
#                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()


