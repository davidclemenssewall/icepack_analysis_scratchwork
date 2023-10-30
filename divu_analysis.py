import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd

icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive"
run_name = "test_adv_net" #"test_netcdf" #"test_clos_adv" #"test_divu_import"#"test_divu_plus" #"test_divu" 

data_path = os.path.join(icepack_dirs_path, "icepack-dirs", "runs", run_name, "history")
filename = os.listdir(data_path)[0]

ds = xr.open_dataset(os.path.join(data_path, filename))

# We need to convert the cftime index to a datetime index for plotting and stuff. Also, note
# that 2020 was a leap year so the Icepack default no leap year output won't work...
# if we use leap year output this is not needed
datetimeindex = ds.indexes['time'].to_datetimeindex()
ds['time'] = datetimeindex

# Now try plotting a variable
ds['divu'].plot.line(x='time')
plt.show()

# Convert to a dataframe
df_divu = ds.sel({'ni':1})['divu'].to_pandas()

# Also read in forcing data
data_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs/input/Icepack_data/forcing/SHEBA"
filename = "open_clos_lindsay.dat"

df = pd.read_csv(os.path.join(data_path, filename), names=['days', 'open', 'clos'], sep='\s+')
df.set_index('days', inplace=True)

# Let's take a look at aice over time
df_plot = pd.concat((ds['aice'].to_pandas().reset_index(), df.reset_index()), 
                    axis=1)
df_plot2 = pd.concat((ds['sst'].to_pandas().reset_index(), df.reset_index()), 
                    axis=1)

f, axs = plt.subplots(3, 1, figsize=(20,10), sharex=True)

df_plot[['time', 2]].plot.line(x='time', ax=axs[0], ylabel='aice')
df_plot2[['time', 2]].plot.line(x='time', ax=axs[1], ylabel='vice')
df_plot[['time', 'open', 'clos']].plot.line(x='time', ax=axs[2], ylabel='Deformation (1/s)')

axs[0].set_xlim(["2015-03-16", "2015-03-26"])
if run_name=="test_adv_net":
    #axs[1].set_ylim([2.6, 2.75])
    axs[1].set_ylim([-1.91, -1.9])
elif run_name=="test_netcdf":
    axs[1].set_ylim([2.7, 2.85])

axs[0].set_ylim([0.985, 0.995])
axs[2].set_ylim([-2e-7, 2e-7])

plt.show()

# Let's look at the change in aice vs opening and closing