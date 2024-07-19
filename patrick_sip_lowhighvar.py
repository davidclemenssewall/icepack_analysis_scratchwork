import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime

import icepacktools as ipt

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mos_raph_syi_ks030": None,
            "mos_sip_lowvar": None,
            "mos_sip_highvar": None,
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
run_plot_dict = {"mos_raph_syi_ks030": [2],
            "mos_sip_lowvar": [2],
            "mos_sip_highvar": [2],
                 }
var_names = ['aice', 'vice', 'vsno', 'fsens', 'flat']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
#                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

# Plot
run_plot_dict = {"mos_raph_syi_ks030": [2],
            "mos_sip_lowvar": [2],
            "mos_sip_highvar": [2],
                 }
var_names = ['evap', 'fsnow', 'frain']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, cumulative=True)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
#                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()