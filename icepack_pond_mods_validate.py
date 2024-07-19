# icepack_pond_mods_validate.py

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime

import icepacktools as ipt

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"pndhist": None,
            "pnd_frbd": None,
            "pnd_hyps": None,
            "aicen_issue_pndhist": None,
            "pndhyps_fixed": None,
            "pndhyps_fixed_frbd": None,
            "mos_raph_fyi_frshbud": None,
            "mos_raph_fyi_frshbud_pndhyps_fixed_frbd": None,
            "pndhead": None,
            "pnd_sealevel_head": None,
            "pndmacr": None,
            "pnd_allmods": None,
            "mos_raph_fyi_pnd_allmods": None,
            "pnd_sealevel_log_all": None
            }

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = ipt.load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value)
    hist_dict[key]['hin'] = hist_dict[key]['vicen']/hist_dict[key]['aicen']


# The following is useful to check that outputs match
if False:
    comp = "pndmacr" #"pndhead" #"aicen_issue_pndhist" #"pnd_hyps" #
    # Check whether or not dataarrays are identical
    for key, da in hist_dict["pnd_hyps"].data_vars.items():
        if not da.equals(hist_dict[comp][key]):
            print(key)
            print('max diff: ' + str((hist_dict[comp][key] - da).max().values))
            print('min diff: ' + str((hist_dict[comp][key] - da).min().values))
    print("Above are data arrays that do not match.")

# Explore state variables
run_plot_dict = {#"pndhist": [1],
                 #"aicen_issue_pndhist": [1],
                 "pnd_hyps": [1],
                 "pnd_frbd": [1],
                 "pndhyps_fixed": [1],
                 "pndhyps_fixed_frbd": [1],
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

# Compare pndhead with pndhyps_fixed
run_plot_dict = {"pndhist": [1, 2, 3],
                 "pndhyps_fixed": [1, 2, 3],
                 "pndhead": [1, 2, 3],
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

# Compare sealevel with pndhead
run_plot_dict = {"pndhead": [1, 2, 3],
                 "pnd_sealevel_head": [1, 2, 3]
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

# Compare pndmacr (defaults except pndmacr='head') with pnd_hyps
run_plot_dict = {"pnd_hyps": [1, 2, 3],
                 "pndmacr": [1, 2, 3]
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

# Compare pndmacr (defaults except pndmacr='head') with pnd_hyps
run_plot_dict = {"pnd_hyps": [1, 2, 3],
                 "pnd_allmods": [1, 2, 3]
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

# Compare sealevel_log
run_plot_dict = {"pnd_hyps": [1, 2, 3],
                 "pnd_sealevel_log_all": [1, 2, 3]
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

run_plot_dict = {
                 "pnd_sealevel_log_all": [3]
                }
var_names = ['aicen', 'vicen', 'hin', 'apndn', 'hpndn', 'ipndn']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
#                  datetime.datetime.fromisoformat('2015-08-15')])
#axs[3].set_ylim([0, 1.0])
plt.show()

run_plot_dict = {
                 "pnd_sealevel_log_all": [2]
                }
var_names = ['flpndn', 'expndn', 'frpndn', 'rfpndn']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict, cumulative=True)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()


# Explore state variables
run_plot_dict = {
                 "mos_raph_fyi_frshbud": [2],
                 "mos_raph_fyi_frshbud_pndhyps_fixed_frbd": [2],
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
#axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
#                  datetime.datetime.fromisoformat('2020-08-01')])
plt.show()

# Explore mosaic case
run_plot_dict = {
                 "mos_raph_fyi_frshbud": [2],
                 "mos_raph_fyi_pnd_allmods": [2],
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-08-01')])
plt.show()

# Plots for PCWG talk
# Explore mosaic case
run_plot_dict = {
                 "mos_raph_fyi_frshbud": [2],
                 #"mos_raph_fyi_pnd_allmods": [2],
                }
var_names = ['vice', 'apnd', 'hpnd', 'ipnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-08-01')])
plt.show()

# Explore state variables
run_plot_dict = {"pndhist": [1],
                 "pnd_mods_aicen_issue": [1],
                 #"pnd_hyps": [1],
                }
var_names = ['aicen', 'apndn', 'hpndn', 'flpndn', 'frpndn']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-06-01'),
                  datetime.datetime.fromisoformat('2015-08-15')])
plt.show()

# Time of spike in frpnd for ni=1
spike_time = "2015-07-03T07:00:00.000000000"
time_minus1 = "2015-07-03T06:00:00.000000000"

# Explore state variables
run_plot_dict = {"pndhist": [1],
                 "pnd_mods_aicen_issue": [1],
                 #"pnd_hyps": [1],
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'frpnd', 'flpnd', 'expnd']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-07-02T13:00:00'),
                  datetime.datetime.fromisoformat('2015-07-03T10:00:00')])
plt.show()

# Explore state variables
run_plot_dict = {"pndhist": [1],
                 "pnd_mods_aicen_issue": [1],
                 #"pnd_hyps": [1],
                }
var_names = ['aicen', 'apndn', 'hpndn', 'flpndn', 'frpndn']

f, axs = ipt.plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2015-07-02T13:00:00'),
                  datetime.datetime.fromisoformat('2015-07-03T10:00:00')])
plt.show()
