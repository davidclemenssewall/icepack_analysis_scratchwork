import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"

# Load met data
forcing_path = os.path.join(icepack_dirs_path, "input", "Icepack_data",
                            "forcing", "MOSAiC")
forc_filename = "MOSAiC_MODF_20191011-20201001_v0.2.nc"
ds_forc = xr.open_dataset(os.path.join(forcing_path, forc_filename))

# Get start and stops of missing data
var = 'hus'

start_stop = []
missing = True
span = np.array([0, 0])

for i in np.arange(ds_forc.time01.size):
    if missing:
        if not np.isnan(ds_forc[var].values[i]):
            span[1] = i
            start_stop.append(span.copy())
            missing = False
    else:
        if np.isnan(ds_forc[var].values[i]):
            span[0] = i
            missing = True

start_stop_arr = np.array(start_stop)
start_stop_dates = ds_forc.time01.values[start_stop_arr]
df_start_stop = pd.DataFrame(start_stop_dates, columns=['start', 'stop'])
df_start_stop['duration'] = df_start_stop.stop - df_start_stop.start
df_start_stop.sort_values('duration', ascending=False, inplace=True)

df_start_stop.to_csv("missing_" + var + ".csv", index=False)

# Plot
f, axs = plt.subplots(5, 1, sharex=True, figsize=(10, 15))

ds_forc['tas'].plot(ax=axs[0])
ds_forc['rld'].plot(ax=axs[1])
ds_forc['uas'].plot(ax=axs[2])
ds_forc['vas'].plot(ax=axs[2])
ds_forc['hus'].plot(ax=axs[3])
ds_forc['rsd'].plot(ax=axs[4])

axs[-1].set_xlim(pd.to_datetime(['2019-11-15', '2019-12-01']))
f.savefig('met_forcing_nov.png', bbox_inches='tight')