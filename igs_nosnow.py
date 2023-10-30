import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"

figure_path = "/home/dcsewall/figures/igs_2023"

plt.rcParams.update({'font.size': 18})
lfont = 14

# Load base case Icepack output
# This run has ktherm=2, init_ice=0.67, init_sno=0.08, qdp=-1.0, ustar_min =0.005
run_name = "igs_nosnow_extice" #"igs_nosnow"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
#filename = os.listdir(data_path)[0]
hist_filename = 'icepack.h.20191129.nc'
ds = xr.open_dataset(os.path.join(hist_path, hist_filename))

# Load met data
forcing_path = os.path.join(icepack_dirs_path, "input", "Icepack_data",
                            "forcing", "MOSAiC")
forc_filename = "MOSAiC_MODF_20191011-20201001_v0.2.nc"
ds_forc = xr.open_dataset(os.path.join(forcing_path, forc_filename))

# Load stakes data
data_path = "/home/dcsewall/data/mass_balance_data"
data_filename = "ablationStakes_hotwireThicknessGauges_MOSAiC.csv"
df_stak = pd.read_csv(os.path.join(data_path, data_filename), parse_dates=[3,4,7])

# Compare with the slab initial condition
ni = 2

# Compare atmospheric fluxes
flux_names = (#('', 'tsns'), # surface temp
              #('', 'rsu'), # outgoing shortwave
              ('flwout', 'rlu'), # outgoing longwave
              ('fsens', 'wthv', 'hfss_ec', 'hfss_bulk'), # upward sensible heat
              ('flat', 'wqv', 'hfls_bulk'), # upward latent heat
             )

legend_names = (
                ('Model', 'Obs.'),
                ('Model', 'Obs. Virt. Pot. T', 'Obs. Eddy Cov.', 'Obs. Bulk'),
                ('Model', 'Obs. Vapor Content', 'Obs. Bulk')
                )

ylabels = ('Longwave\nOut (W m$^{-2}$)', 'Turbulent Sensible\n(W m$^{-2}$)',
           'Turbulent\nLatent (W m$^{-2}$)')

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
        
    axs[i].set_ylabel(ylabels[i])
    axs[i].grid()
    axs[i].legend(legend_names[i], fontsize=lfont)

axs[0].set_xlim(["2019-11-28", ds.time.values[-1]])
axs[-1].set_xlabel('')

# Compare ice state
stakes_site = "Stakes 3/dart_stakes_clu_3"

stat_names = (('vice', 'Ice thickness (calculated) (cm)'), # ice thickness
              ('vsno', 'Snow depth (calculated) (cm)'), # snow depth
             )

ylabels = ('Ice Thickness (m)', 'Snow Depth (m)')

# Plot comparison
f, axs = plt.subplots(len(stat_names), 1, sharex=True, figsize=(10, 3*len(stat_names)))

for i in range(len(stat_names)):
    ds[stat_names[i][0]].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', 
                                                                        ax=axs[i],
                                                                        label='Icepack')
    (df_stak[df_stak["Site name"]==stakes_site].groupby('Measurement date'
        )[stat_names[i][1]].agg(['mean', 'sem'])/100).plot.line(y='mean',ax=axs[i], linestyle=':',
                                                  #marker='o', 
                                                  yerr='sem', capsize=5)

    axs[i].set_ylabel(ylabels[i])
    axs[i].grid()
    axs[i].legend(['Model', 'Obs.'], fontsize=lfont)

axs[0].set_xlim(["2019-11-28", ds.time.values[-1]])
#axs[0].set_xlim([ds.time.values[0], ds.time.values[23]])
axs[-1].set_xlabel('')

# Sensitivity to initial conditions
# 10 cm too thin
run_name = "igs_nosnow_extice_thin" #"igs_nosnow_thinice"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_thin = xr.open_dataset(os.path.join(hist_path, hist_filename))
# 10 cm too thick
run_name = "igs_nosnow_extice_thick" # "igs_nosnow_thickice"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_thick = xr.open_dataset(os.path.join(hist_path, hist_filename))

# Plot comparison
f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6))

ds['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)
(df_stak[df_stak["Site name"]==stakes_site].groupby('Measurement date'
    )['Ice thickness (calculated) (cm)'].agg(['mean', 'sem'])/100).plot.line(y='mean',ax=axs, 
                                                                             linestyle=':',
                                                                             yerr='sem', capsize=5)
ds_thin['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)
ds_thick['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)

axs.set_ylabel('Ice Thickness (m)')
axs.grid()
axs.legend(['Model (init. 0.67 m)', 'Model (init. 0.57 m)', 'Model (init. 0.77 m)', 'Obs.'], fontsize=lfont)

axs.set_xlim(["2019-11-28", ds.time.values[-1]])
axs.set_ylim([1.55, 2.65])
axs.set_xlabel('')

# no ice vs. 10 cm case
# 10 cm too thin
run_name = "igs_snow_10cm" #"igs_nosnow_thinice"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_10cm = xr.open_dataset(os.path.join(hist_path, hist_filename))

# Plot comparison
f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6))

ds_10cm['vice'].sel(ni=1).to_pandas().reset_index().plot.line(x='time', ax=axs)
(df_stak[df_stak["Site name"]==stakes_site].groupby('Measurement date'
    )['Ice thickness (calculated) (cm)'].agg(['mean', 'sem'])/100).plot.line(y='mean',ax=axs, 
                                                                             linestyle=':',
                                                                             yerr='sem', capsize=5)
ds_10cm['vice'].sel(ni=2).to_pandas().reset_index().plot.line(x='time', ax=axs)

axs.set_ylabel('Ice Thickness (m)')
axs.grid()
axs.legend(['Model (init. 0.0 m)', 'Model (init. 0.1 m)', 'Obs.'], fontsize=lfont)

axs.set_xlim(["2019-11-28", ds.time.values[-1]])
#axs.set_ylim([0, 1.2])
axs.set_xlabel('')

print(ds_10cm['vice'].sel(ni=1).values[-1])
print(ds_10cm['vice'].sel(ni=2).values[-1])