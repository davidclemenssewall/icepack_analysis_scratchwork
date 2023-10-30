import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd

icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"

# Load Icepack output
run_name = "mosaic_raphael_syi" #"mosaic_atm_ocn" #"gap_fill_atm" #"test_igs"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = os.listdir(hist_path)[0]
hist_filename = 'icepack.h.20191115.nc'
hist_filename = 'icepack.h.20191129.nc'
ds = xr.open_dataset(os.path.join(hist_path, hist_filename))

# Load met data
forcing_path = os.path.join(icepack_dirs_path, "input", "Icepack_data",
                            "forcing", "MOSAiC")
forc_filename = "MOSAiC_MODF_20191011-20201001_v0.2.nc"
ds_forc = xr.open_dataset(os.path.join(forcing_path, forc_filename))
ds_forc['hus'] = ds_forc['hus']/1000

# Load ocean data
ocn_filename = 'MOSAiC_ocn_MDF_20191007-20200920.nc'
ds_ocn = xr.open_dataset(os.path.join(forcing_path, ocn_filename))


# Load stakes data
data_path = "/home/dcsewall/data/mass_balance_data"
data_filename = "ablationStakes_hotwireThicknessGauges_MOSAiC.csv"
df_stak = pd.read_csv(os.path.join(data_path, data_filename), parse_dates=[3,4,7])


# Compare forcing
forc_names = (('flw', 'rld'),
              ('fsw', 'rsd'),
              ('Tair', 'tas'),
              ('Qa', 'hus'),
             )
ni = 2
# Plot comparison
f, axs = plt.subplots(len(forc_names), 1, sharex=True, figsize=(10, 3*len(forc_names)))

for i in range(len(forc_names)):
    ds[forc_names[i][0]].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', 
                                                                        ax=axs[i],
                                                                        label='Icepack')
    ds_forc[forc_names[i][1]].to_pandas().reset_index().plot.line(x='time01', ax=axs[i],
                                                                  label='MET', linestyle=':')
    axs[i].set_ylabel(forc_names[i][0])
    axs[i].grid()
    axs[i].legend(forc_names[i])

axs[0].set_xlim([ds.time.values[0], ds.time.values[-1]])
#axs[0].set_xlim([ds.time.values[0], ds.time.values[23]])

plt.show()
# looks like about a 3 hour offset, need to assess if this is due to time basis errors, 
# interpolation choices, or how Icepack records the time at which things happen

# Compare atmospheric fluxes
flux_names = (#('', 'tsns'), # surface temp
              #('', 'rsu'), # outgoing shortwave
              ('flwout', 'rlu'), # outgoing longwave
              ('fsens', 'wthv', 'hfss_ec', 'hfss_bulk'), # upward sensible heat
              ('flat', 'wqv', 'hfls_bulk'), # upward latent heat
             )
# Plot comparison
f, axs = plt.subplots(len(flux_names), 1, sharex=True, figsize=(10, 3*len(flux_names)))

for i in range(len(flux_names)):
    ds[flux_names[i][0]].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', 
                                                                        ax=axs[i],
                                                                        label='Icepack')
    if len(flux_names[i]) == 2:
        (-1*ds_forc[flux_names[i][1]]).to_pandas().reset_index().plot.line(x='time01', ax=axs[i],
                                                                    label='MET', linestyle=':')
    elif len(flux_names[i]) > 2:
        for j in range(len(flux_names[i])-1):
            (-1*ds_forc[flux_names[i][j+1]]).to_pandas().reset_index().plot.line(x='time10', ax=axs[i],
                                                                        label='MET', linestyle=':',
                                                                        alpha=0.5)
        
    axs[i].set_ylabel(flux_names[i][0])
    axs[i].grid()
    axs[i].legend(flux_names[i])

axs[0].set_xlim([ds.time.values[0], ds.time.values[-1]])
#axs[0].set_xlim([ds.time.values[0], ds.time.values[23]])

plt.show()

# Compare averages
for i in range(len(flux_names)):
    print(flux_names[i][0])
    print(ds[flux_names[i][0]].sel(ni=ni).mean().values)
    for j in range(len(flux_names[i])-1):
        if flux_names[i][j+1] == 'rlu':
            print(flux_names[i][j+1])
            print(-1*ds_forc[flux_names[i][j+1]].sel(time01=slice(ds.time.values[0], 
                                                            ds.time.values[-1])).mean().values)
        else:
            print(flux_names[i][j+1])
            print(-1*ds_forc[flux_names[i][j+1]].sel(time10=slice(ds.time.values[0], 
                                                            ds.time.values[-1])).mean().values)

# Compare oceanic values
# Add variables for sst above freezing
ds['sst_above_frz'] = ds['sst'] - ds['Tf']
ds_ocn['sst_above_frz'] = ds_ocn['tos'] - ds_ocn['tosf']
flux_names = (
              ('fhocn',), # Oceanic heat flux
              ('sst', 'tos'), # sst
              ('sss', 'so'), # sss
              ('Tf', 'tosf'), # ocean freezing temperature
              ('sst_above_frz', 'sst_above_frz'), # t above freezing point
             )
# Plot comparison
f, axs = plt.subplots(len(flux_names), 1, sharex=True, figsize=(10, 3*len(flux_names)))

for i in range(len(flux_names)):
    df_plot = ds[flux_names[i][0]].sel(ni=ni).to_pandas()
    df_plot.name = flux_names[i][0]
    axs[i].plot(df_plot)
    #ds[flux_names[i][0]].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', 
    #                                                                    ax=axs[i],
    #                                                                    label='Icepack')
    if len(flux_names[i]) == 2:
        df_plot = ds_ocn[flux_names[i][1]].to_pandas()
        df_plot.name = flux_names[i][1]
        axs[i].plot(df_plot, linestyle=':')
        
    axs[i].set_ylabel(flux_names[i][0])
    axs[i].grid()

axs[0].set_xlim([ds.time.values[0], ds.time.values[-1]])
#axs[0].set_xlim([ds.time.values[0], ds.time.values[23]])

plt.show()

print(ds['fhocn'].sel(ni=ni).mean().values)

# Compare mass balance terms
flux_names = (
              ('frazil',), # 
              ('congel',), # 
              ('snoice',), # 
              ('meltb',), # 
              ('meltt',),
              ('meltl',),
              ('evap',),
             )
# Plot comparison
f, axs = plt.subplots(len(flux_names), 1, sharex=True, figsize=(10, 3*len(flux_names)))

for i in range(len(flux_names)):
    ds[flux_names[i][0]].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', 
                                                                        ax=axs[i],
                                                                        label='Icepack')
    if len(flux_names[i]) == 2:
        (-1*ds_forc[flux_names[i][1]]).to_pandas().reset_index().plot.line(x='time01', ax=axs[i],
                                                                    label='MET', linestyle=':')
    elif len(flux_names[i]) > 2:
        for j in range(len(flux_names[i])-1):
            (-1*ds_forc[flux_names[i][j+1]]).to_pandas().reset_index().plot.line(x='time10', ax=axs[i],
                                                                        label='MET', linestyle=':',
                                                                        alpha=0.3)
        
    axs[i].set_ylabel(flux_names[i][0])
    axs[i].grid()


axs[0].set_xlim([ds.time.values[0], ds.time.values[-1]])
#axs[0].set_xlim([ds.time.values[0], ds.time.values[23]])

plt.show()

# Compare sums
for i in range(len(flux_names)):
    print(flux_names[i][0])
    print(ds[flux_names[i][0]].sel(ni=ni).sum().values)


# Compare ice state
if ni == 1:
    stakes_site = 'Ridge Ranch/dart_stakes_clu_6' #"Stakes 3/dart_stakes_clu_3"
else:
    stakes_site = "Stakes 3/dart_stakes_clu_3" #'MET Stakes/dart_stakes_clu_5' #'Ridge Ranch/dart_stakes_clu_6' #

stat_names = (('vice', ('Ice thickness (calculated) (cm)',
                        'Drilled ice thickness (cm)')), # ice thickness
              ('vsno', ('Snow depth (calculated) (cm)',)), # snow depth
             )

# What if we limit to just the stakes that we drilled at the end?
last_stakes = df_stak[(df_stak['Measurement date']=='2020-05-01') & 
                      (~df_stak['Drilled ice thickness (cm)'].isna())]['Stake ID'].values
use_last_stakes = False

# Plot comparison
f, axs = plt.subplots(len(stat_names), 1, sharex=True, figsize=(10, 3*len(stat_names)))

for i in range(len(stat_names)):
    ds[stat_names[i][0]].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', 
                                                                        ax=axs[i],
                                                                        label='Icepack')
    for stat_name in stat_names[i][1]:
        if use_last_stakes:
            cond = (df_stak["Site name"]==stakes_site) & df_stak['Stake ID'].isin(last_stakes)
        else:
            cond = df_stak["Site name"]==stakes_site
        (df_stak[cond].groupby('Measurement date'
            )[stat_name].agg(['mean', 'sem'])/100).plot.line(y='mean',ax=axs[i], linestyle=':',
                                                    #marker='o', 
                                                    yerr='sem', capsize=5,
                                                    label=stat_name)

    axs[i].set_ylabel(stat_names[i][0])
    axs[i].grid()
    axs[i].legend()

axs[0].set_xlim([ds.time.values[0], ds.time.values[-1]])
#axs[0].set_xlim([ds.time.values[0], ds.time.values[23]])

plt.show()

# Look at pond variables
# trcr dict contains the trcr indices from log file. 
# Check fortran indexing!
trcr_dict = {'apnd': 19,
             'hpnd': 20,
             'ipnd': 21}

# Plot tracers
f, axs = plt.subplots(len(trcr_dict), 1, sharex=True, figsize=(10, 3*len(trcr_dict)))

for i, key in enumerate(trcr_dict):
    ds['trcr'].sel(ni=ni, ntrcr=trcr_dict[key]).to_pandas().reset_index().plot.line(x='time', 
                                                                        ax=axs[i],
                                                                        label='Icepack')

    axs[i].set_ylabel(key)
    axs[i].grid()
    axs[i].legend()

axs[0].set_xlim([ds.time.values[0], ds.time.values[-1]])
#axs[0].set_xlim([ds.time.values[0], ds.time.values[23]])

plt.show()

print(ds['vice'].sel(ni=ni).isel({'time':-1}).values)
print(ds['aicen'].sel(ni=ni).isel({'time':-1}).values)
print(ds['vicen'].sel(ni=ni).isel({'time':-1}).values/ds['aicen'].sel(ni=ni).isel({'time':-1}).values)

df_ocn = ds[['fhocn','congel','frazil']].sel(ni=ni).to_pandas()
df_ocn['growth'] = df_ocn['congel'] + df_ocn['frazil']

f, ax = plt.subplots(1, 1)

df_ocn.plot.scatter(x='fhocn', y='congel', ax=ax, label='congel', c='b')
df_ocn.plot.scatter(x='fhocn', y='frazil', ax=ax, label='frazil', c='g')
df_ocn.plot.scatter(x='fhocn', y='growth', ax=ax, label='total', c='r')
ax.legend()
plt.show()