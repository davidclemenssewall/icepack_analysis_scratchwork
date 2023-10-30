import xarray as xr
import os
import matplotlib.pyplot as plt

icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive"

filename = "adv_diff_2.nc"

ds = xr.open_dataset(os.path.join(icepack_dirs_path, filename))

# We need to convert the cftime index to a datetime index for plotting and stuff. Also, note
# that 2020 was a leap year so the Icepack default no leap year output won't work...
# if we use leap year output this is not needed
datetimeindex = ds.indexes['time'].to_datetimeindex()
ds['time'] = datetimeindex

# Which variables to plot
vars_to_plt = [
               'aice', 'vice', 'vsno', 
               #'evap', 'fswabs', 'flwout', 'fsens', 'Tair', 'fcondtop', 
               #'meltt', 'meltb', 'meltl', 'snoice', 'congel', 'frazil',
               #'sst', 'fhocn'
               ]

# Now plot variables
f, axs = plt.subplots(len(vars_to_plt), 1, sharex=True, figsize=(10, 2*len(vars_to_plt)))

for i in range(len(vars_to_plt)):
    ds[vars_to_plt[i]].plot.line(x='time', ax=axs[i])

plt.show()