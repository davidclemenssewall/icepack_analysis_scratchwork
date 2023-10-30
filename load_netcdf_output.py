import xarray as xr
import os
import matplotlib.pyplot as plt

icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive"
run_name = "test_netcdf" #"test_start_date"

data_path = os.path.join(icepack_dirs_path, "icepack-dirs", "runs", run_name, "history")
filename = os.listdir(data_path)[0]

ds = xr.open_dataset(os.path.join(data_path, filename))

# We need to convert the cftime index to a datetime index for plotting and stuff. Also, note
# that 2020 was a leap year so the Icepack default no leap year output won't work...
# if we use leap year output this is not needed
datetimeindex = ds.indexes['time'].to_datetimeindex()
ds['time'] = datetimeindex

# Now try plotting a variable
ds['vice'].plot.line(x='time')
plt.show()

ds['fsw'].plot.line(x='time')
plt.show()

