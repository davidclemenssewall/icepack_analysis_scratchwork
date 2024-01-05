"""
Figure for contrasts proposal

"""

import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import datetime

# Function for reading history
def load_icepack_hist(run_name, icepack_dirs_path, hist_filename=None,
                      sst_above_frz=True, trcr_dict=None, trcrn_dict=None):
    """
    Load Icepack history output
    
    Parameters
    ----------
    run_name : str
        Name of the icepack run (directory name in RUN_DIR)
    icepack_dirs_path : str
        Path to root of icepack directory.
    hist_filename : str or None, optional
        Name of specific history file to load. If None load the first file 
        in history directory. Default is None.
    sst_above_frz : bool, optional
        Whether or not to compute the difference between mixed layer freezing 
        point and temperature. Default is True.
    trcr_dict : dict, optional
        Dict for which tracers to convert to data variables. Keys are tracer
        indices and values are names. Default is None.
    trcrn_dict : dict, optional
        Dict for which category tracers to convert to data variables. Keys 
        are tracer indices and values are names. Default is None.

    Returns
    -------
    xarray dataset with Icepack history output

    """

    # Open netCDF
    hist_path = os.path.join(icepack_dirs_path, "runs", run_name, "history")
    if hist_filename is None:
        hist_filename = os.listdir(hist_path)[0]
    ds = xr.open_dataset(os.path.join(hist_path, hist_filename))

    # Create mixed layer freezing point difference
    if sst_above_frz:
        ds['sst_above_frz'] = ds['sst'] - ds['Tf']
    
    # Copy trcr and trcrn data variables
    if trcr_dict is not None:
        for key, value in trcr_dict.items():
            da = ds['trcr'].sel(ntrcr=key)
            da.name = value
            ds[value] = da
    if trcrn_dict is not None:
        for key, value in trcrn_dict.items():
            da = ds['trcrn'].sel(ntrcr=key)
            da.name = value
            ds[value] = da

    # Add the run name as an attribute
    ds.attrs.update({'run_name': run_name})

    return ds

# Function for plotting single Icepack output
def plot_hist_var(ds, var_name, ni, ax, resample=None,
                  cumulative=False, dt=None):
    """
    Plot a single variable from history DS on the given axis

    Parameters
    ----------
    ds : xarray DataSet
    var_name : str
    ni : int
        Which cell of the Icepack output to plot
    ax : matplotlib.pyplot.axes
        Axis object to plot on
    resample : str, optional
        If provided, frequency string for DataFrame.resample(). If None do not
        resample. The default is None.
    cumulative : bool, optional
        Whether the variable should be cumulative, useful for fluxes
    
    Returns
    -------
    handle for matplotlib plot object

    """

    # Get variable as Pandas DataFrame with time as a column
    df = ds[var_name].sel(ni=ni).to_pandas()
    if resample is not None:
        df = df.resample(resample).mean()
    df.reset_index(inplace=True)
    if df.shape[1] == 2:
        df.rename(columns={0: var_name}, inplace=True)
        label = ds.run_name + ' (' + str(ni) + ')'
        # Plot
        h = ax.plot('time', var_name, data=df, label=label)
    else:
        for col_name in df.columns:
            if col_name == 'time':
                continue
            label = ds.run_name + ' (' + str(ni) + ', ' + str(col_name) + ')'
            # Plot
            h = ax.plot('time', col_name, data=df, label=label)

    return h

def plot_forc_var(ds_forc, var_name, ax):
    """
    Plot a single forcing variable

    Parameters
    ----------
    ds : xarray DataSet
        MDF formatted
    var_name : str
    ax : matplotlib.pyplot.axes
        Axis object to plot on
    
    Returns
    -------
    handle for matplotlib plot object

    """

    # Get variable as Pandas DataFrame with time as a column
    df = ds_forc[var_name].to_pandas().reset_index()
    time_name = df.columns[0]
    df.rename(columns={time_name: 'time', 0: var_name}, inplace=True)
    # Plot
    h = ax.plot('time', var_name, data=df, linestyle=':', alpha=0.5)

    return h

def plot_ice_var(df_ice, var_name, site, ax, mean_only=False):
    """
    Plots the requested ice variable from an individual site
    
    Parameters
    ----------
    df_ice : Pandas dataframe
    var_name : str
    site : str
    ax :


    """

    # Get dataframe with columns, time, mean, sem
    df = df_ice.loc[site, var_name].reset_index()
    label = site + " " + var_name
    # plot
    if mean_only or df['sem'].isna().all():
        h = ax.plot('time', 'mean', data=df, linestyle=':', label=label,
                    marker='o')
    else:
        h = ax.errorbar('time', 'mean', yerr='sem', data=df, 
                        capsize=5, linestyle=':', label=label)

    return h

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mosaic_raphael_syi": 'icepack.h.20191129.nc',
            "mosaic_raphael_fyi": 'icepack.h.20191129.nc',
            "mosaic_raphael_fyi_rfracaice": None,
            "mosaic_raphael_fyi_nopndfbd": None,
            "mosaic_raphael_fyi_nopndfbd_pndaspect4": None,
            }
trcr_dict = {19: 'apnd',
             20: 'hpnd',
             21: 'ipnd'}
trcrn_dict = {19: 'apndn',
              20: 'hpndn',
              21: 'ipndn'}

hist_dict = {}
for key, value in run_dict.items():
    hist_dict[key] = load_icepack_hist(run_name=key, 
                                       icepack_dirs_path=ip_dirs_path, 
                                       hist_filename=value, 
                                       trcr_dict=trcr_dict, 
                                       trcrn_dict=trcrn_dict)

# Create ice thickness change
for value in hist_dict.values():
    value['vice_chg'] = value['vice'] - value['vice'].isel(time=0)
    value['pndaspect'] = value['hpnd'] / value['apnd']

# Load ice properties
data_path = "/home/dcsewall/data/mass_balance_data"
data_filename = "ablationStakes_hotwireThicknessGauges_MOSAiC.csv"
df_stak = pd.read_csv(os.path.join(data_path, data_filename), parse_dates=[3,4,7])

# We want to create a dataframe where the row index is site and date and the
# column index is stat_name, and then mean and standard error of mean.
# rename columns to make our life easier
df_renamed = df_stak.rename(columns={'Site name': 'site', 
                                     'Measurement date': 'time',
                                     'Ice thickness (calculated) (cm)': 'hi_hotwire',
                                     'Drilled ice thickness (cm)': 'hi_drill',
                                     'Snow depth (calculated) (cm)': 'hs_gauge',
                                     'Pond depth (cm)': 'hpnd',
                                     'Pond flag': 'apnd'})
# cm to m
for var_name in ['hi_hotwire', 'hi_drill', 'hs_gauge', 'hpnd']:
    df_renamed[var_name] = df_renamed[var_name]/100
stat_names = ['hi_hotwire','hi_drill','hs_gauge','hpnd','apnd']
df_ice = df_renamed.groupby(['site', 'time'])[stat_names].agg(
    ['mean', 'sem'])

ice_var_map = {'vice': ['hi_hotwire', 'hi_drill'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }

# Compare with changing pond parameterizations
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "mosaic_raphael_fyi_nopndfbd": [2],
                 "mosaic_raphael_fyi_nopndfbd_pndaspect4": [2],
                 "mosaic_raphael_fyi_rfracaice": [2],
                 }
var_names = ['vice', 'apnd', 'hpnd']
ice_var_map = {'vice': ['hi_hotwire'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }
ice_sites = ['Reunion Stakes/dart_stakes_clu_12']

xlim = [datetime.datetime.fromisoformat('2020-05-24'),
                  datetime.datetime.fromisoformat('2020-07-15')]
axfont = 14
lfont = 8
dfmt = mdates.DateFormatter("%m-%d")

f, axs = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)

# Loop through each variable
for var_name, ax in zip(var_names, axs):
    # And through each run
    for run_name, nis in run_plot_dict.items():
        # and the desired cell(s) in each run
        for ni in nis:
            h = plot_hist_var(hist_dict[run_name], var_name, ni, ax)
            if run_name != "mosaic_raphael_fyi":
                h[0].set(linestyle=':')
    # Plot site variables
    for site in ice_sites:
        if var_name in ice_var_map:
            for ice_var_name in ice_var_map[var_name]:
                _ = plot_ice_var(df_ice, ice_var_name, site, ax, 
                                    mean_only=True)

for ax in axs:
    ax.set_xlim(xlim)
    ax.grid(alpha=0.5, lw=1)
    ax.xaxis.set_major_formatter(dfmt)
    ax.set_xticks(ax.get_xticks()[::2])

axs[0].set_ylabel('Ice Thickness (m)', fontsize=axfont)
axs[1].set_ylabel('Pond Area Fraction', fontsize=axfont)
axs[2].set_ylabel('Pond Depth (m)', fontsize=axfont)

axs[0].set_ylim([0.7, 1.85])
axs[0].legend(["Icepack Control", "Freeboard", 
               "Depth-to-Area Ratio & Freeboard", "Runoff", 
               "Field Observations"],
               fontsize=lfont)
plt.show()

