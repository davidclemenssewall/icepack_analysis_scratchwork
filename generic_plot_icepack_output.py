import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd

#%matplotlib widget

icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive"
run_name = "base_run" #"test_mosaic_forcing" #"test_adv_none" #"test_no_adv" #"test_adv_net" #"test_netcdf" #"test_clos_adv_inline" # #"test_divu_import"#"test_divu_plus" #"test_divu" 

data_path = os.path.join(icepack_dirs_path, "icepack-dirs", "runs", run_name, "history")
filename = os.listdir(data_path)[0]
if run_name == "test_igs":
    filename = os.listdir(data_path)[1] # for test_igs

ds = xr.open_dataset(os.path.join(data_path, filename))

# We need to convert the cftime index to a datetime index for plotting and stuff. Also, note
# that 2020 was a leap year so the Icepack default no leap year output won't work...
# if we use leap year output this is not needed
datetimeindex = ds.indexes['time'].to_datetimeindex()
ds['time'] = datetimeindex

# SHEBA Ice area from Perovich et al 2002
sheba = pd.DataFrame({ "time": pd.to_datetime(
                                ["1998-05-17", "1998-05-20", "1998-06-10",
                                "1998-06-15", "1998-06-18", "1998-06-22",
                                "1998-06-30", "1998-07-08", "1998-07-15",
                                "1998-07-20", "1998-07-25", "1998-08-07",
                                "1998-08-22", "1998-09-11", "1998-10-04"]
                                ),
                       "aice0":[0.03, 0.021, 0.057,
                                0.031, 0.031, 0.028,
                                0.051, 0.038, 0.058,
                                0.053, 0.050, 0.185,
                                0.182, 0.126, 0.005]})
sheba['aice'] = 1 - sheba['aice0']
sheba['time'] += pd.offsets.DateOffset(years=17)

plt_sheba = False

# Which variables to plot
vars_to_plt = [
               'aice', 'vice', 'vsno', 
               'Tair', 'flw', 'fsw', 'Qa',
               'evap', 'fswabs', 'flwout', 'fsens', 'fcondtop', 
               'meltt', 'meltb', 'meltl', 'snoice', 'congel',
               'sst', 'sss', 'frazil', 'fhocn',
               ]

#t_lims = ["2015-04-01", "2015-11-01"]

# Now plot variables
f, axs = plt.subplots(len(vars_to_plt), 1, sharex=True, figsize=(10, 3*len(vars_to_plt)))

for i in range(len(vars_to_plt)):
    if plt_sheba:
        if vars_to_plt[i] in sheba.keys():
            axs[i].plot('time', vars_to_plt[i], 'k+-', data=sheba)
    ds[vars_to_plt[i]].to_pandas().reset_index().plot.line(x='time', ax=axs[i])
    axs[i].set_ylabel(vars_to_plt[i])
    axs[i].grid()



#axs[0].set_xlim(t_lims)
axs[0].set_ylim([0.7, 1.0])

#axs[2].set_ylim([0, 0.005])

plt.show()

# Get a second run to compare
run_name2 = "r_snw_neg"

data_path = os.path.join(icepack_dirs_path, "icepack-dirs", "runs", run_name2, "history")
filename = os.listdir(data_path)[0]
if run_name == "test_igs":
    filename = os.listdir(data_path)[1] # for test_igs

ds2 = xr.open_dataset(os.path.join(data_path, filename))

# We need to convert the cftime index to a datetime index for plotting and stuff. Also, note
# that 2020 was a leap year so the Icepack default no leap year output won't work...
# if we use leap year output this is not needed
datetimeindex = ds2.indexes['time'].to_datetimeindex()
ds2['time'] = datetimeindex

# Which variables to plot
vars_to_plt = [
               'aice', 'vice', 'vsno', 
               #'Tair', 'flw', 'fsw', 'Qa',
               #'evap', 'fswabs', 'flwout', 'fsens', 'fcondtop', 
               #'meltt', 'meltb', 'meltl', 'snoice', 'congel',
               #'sst', 'sss', 'frazil', 'fhocn',
               ]

#t_lims = ["2015-04-01", "2015-11-01"]

ni = 3

# Now plot variables
f, axs = plt.subplots(len(vars_to_plt), 1, sharex=True, figsize=(10, 3*len(vars_to_plt)))

for i in range(len(vars_to_plt)):
    if plt_sheba:
        if vars_to_plt[i] in sheba.keys():
            axs[i].plot('time', vars_to_plt[i], 'k+-', data=sheba)
    ds[vars_to_plt[i]].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs[i])
    ds2[vars_to_plt[i]].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs[i], linestyle='--')
    axs[i].set_ylabel(vars_to_plt[i])
    axs[i].grid()
    axs[i].legend([run_name, run_name2])