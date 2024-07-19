# icepack_mosaic_asp_compare.py

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import datetime

import icepacktools as ipt

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mos_raph_fyi_ks030_newforc": None,
            "mos_raph_fyi_ks030_newforc_asp": None,
            }

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = ipt.load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value)
    hist_dict[key]['hin'] = hist_dict[key]['vicen']/hist_dict[key]['aicen']
