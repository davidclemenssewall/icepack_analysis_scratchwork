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
run_name = "test_igs"
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

f.savefig(os.path.join(figure_path, 'atm_flux_igs.png'), bbox_inches='tight')

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

f.savefig(os.path.join(figure_path, 'base_ice_igs.png'), bbox_inches='tight')

# Sensitivity to initial conditions
# 10 cm too thin
run_name = "igs_too_thin"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_thin = xr.open_dataset(os.path.join(hist_path, hist_filename))
# 10 cm too thick
run_name = "igs_too_thick"
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
axs.set_ylim([0.55, 1.35])
axs.set_xlabel('')

f.savefig(os.path.join(figure_path, 'init_cond_igs.png'), bbox_inches='tight')

# Sensitivity to oceanic forcing
# default value of 1 W/m2 representative of thermocline
# qdp 0.1 representative of halocline
run_name = "igs_small_qdp"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_small_qdp = xr.open_dataset(os.path.join(hist_path, hist_filename))
# qdp 3 commonly used value
run_name = "igs_large_qdp"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_large_qdp = xr.open_dataset(os.path.join(hist_path, hist_filename))

# Plot comparison
f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6))

ds['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)
(df_stak[df_stak["Site name"]==stakes_site].groupby('Measurement date'
    )['Ice thickness (calculated) (cm)'].agg(['mean', 'sem'])/100).plot.line(y='mean',ax=axs, 
                                                                             linestyle=':',
                                                                             yerr='sem', capsize=5)
ds_small_qdp['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)
ds_large_qdp['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)

axs.set_ylabel('Ice Thickness (m)')
axs.grid()
axs.legend(['Model (qdp = 1 W m$^{-2}$)', 'Model (qdp = 0.1 W m$^{-2}$)', 
            'Model (qdp = 3 W m$^{-2}$)', 'Obs.'], fontsize=lfont)

axs.set_xlim(["2019-11-28", ds.time.values[-1]])
axs.set_ylim([0.55, 1.35])
axs.set_xlabel('')

f.savefig(os.path.join(figure_path, 'qdp_igs.png'), bbox_inches='tight')

# Sensitivity to ksno
# low is ksno 0.2
run_name = "igs_ksno_low"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_ksno_low = xr.open_dataset(os.path.join(hist_path, hist_filename))
# high is ksno 0.4
run_name = "igs_ksno_high"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_ksno_high = xr.open_dataset(os.path.join(hist_path, hist_filename))

# Plot comparison
f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6))

ds['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)
(df_stak[df_stak["Site name"]==stakes_site].groupby('Measurement date'
    )['Ice thickness (calculated) (cm)'].agg(['mean', 'sem'])/100).plot.line(y='mean',ax=axs, 
                                                                             linestyle=':',
                                                                             yerr='sem', capsize=5)
ds_ksno_low['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)
ds_ksno_high['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)

axs.set_ylabel('Ice Thickness (m)')
axs.grid()
axs.legend(['Model (ksno = 0.3 W m$^{-1}$ K$^{-1}$)', 'Model (ksno = 0.2 W m$^{-1}$ K$^{-1}$)', 
            'Model (ksno = 0.4 W m$^{-1}$ K$^{-1}$)', 'Obs.'], fontsize=lfont)

axs.set_xlim(["2019-11-28", ds.time.values[-1]])
axs.set_ylim([0.55, 1.35])
axs.set_xlabel('')

f.savefig(os.path.join(figure_path, 'ksno_igs.png'), bbox_inches='tight')

# Sensitivity to parameterization
# Bitz and Lipscomb thermodynamics
run_name = "igs_ktherm_1"
hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
hist_filename = 'icepack.h.20191129.nc'
ds_ktherm = xr.open_dataset(os.path.join(hist_path, hist_filename))

# Plot comparison
f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6))

ds['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)
(df_stak[df_stak["Site name"]==stakes_site].groupby('Measurement date'
    )['Ice thickness (calculated) (cm)'].agg(['mean', 'sem'])/100).plot.line(y='mean',ax=axs, 
                                                                             linestyle=':',
                                                                             yerr='sem', capsize=5)
ds_ktherm['vice'].sel(ni=ni).to_pandas().reset_index().plot.line(x='time', ax=axs)

axs.set_ylabel('Ice Thickness (m)')
axs.grid()
axs.legend(['Model (mushy thermo)', 'Model (BL99 thermo)', 
            'Obs.'], fontsize=lfont)

axs.set_xlim(["2019-11-28", ds.time.values[-1]])
axs.set_ylim([0.55, 1.35])
axs.set_xlabel('')

f.savefig(os.path.join(figure_path, 'ktherm_igs.png'), bbox_inches='tight')

# Offsetting errors
# Show how congelation and frazil combine
df = ds[['vice', 'congel', 'frazil']].sel(ni=ni).to_pandas().reset_index()
init_hi = 0.67
df['cum_congel'] = init_hi + np.cumsum(df.congel)

f, axs = plt.subplots(1, 1, figsize=(10,6))

axs.fill_between('time', init_hi, 'cum_congel', data=df, color='tab:blue')
axs.fill_between('time', 'cum_congel', 'vice', data=df, color='tab:green')
(df_stak[df_stak["Site name"]==stakes_site].groupby('Measurement date'
    )['Ice thickness (calculated) (cm)'].agg(['mean', 'sem'])/100).plot.line(y='mean',ax=axs, 
                                                                             linestyle=':',
                                                                             yerr='sem', capsize=5,
                                                                             c='tab:orange',
                                                                             linewidth=5)

axs.set_ylabel('Ice Thickness (m)')
axs.grid()
axs.legend(['Model congelation growth', 'Model frazil growth', 
            'Obs.'], fontsize=lfont)

axs.set_xlim(["2019-11-28", ds.time.values[-1]])
axs.set_ylim([0.55, 1.35])
axs.set_xlabel('')

f.savefig(os.path.join(figure_path, 'congel_frazil_igs.png'), bbox_inches='tight')

