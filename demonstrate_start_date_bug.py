import xarray as xr
import os
import matplotlib.pyplot as plt

icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive"
run_name = "netcdf_start_date"

data_path = os.path.join(icepack_dirs_path, "icepack-dirs", "runs", run_name, "history")
filename = os.listdir(data_path)[1]

ds = xr.open_dataset(os.path.join(data_path, filename))

# We need to convert the cftime index to a datetime index for plotting and stuff. Also, note
datetimeindex = ds.indexes['time'].to_datetimeindex()
ds['time'] = datetimeindex

ds['fsw'].plot.line(x='time')
plt.show()