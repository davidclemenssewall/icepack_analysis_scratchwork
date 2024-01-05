"""
Tools for comparing Icepack runs with each other and observations

"""

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
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
def plot_hist_var(ds, var_name, ni, ax):
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
    
    Returns
    -------
    handle for matplotlib plot object

    """

    # Get variable as Pandas DataFrame with time as a column
    df = ds[var_name].sel(ni=ni).to_pandas().reset_index()
    if df.shape[1] == 2:
        df.rename(columns={0: var_name}, inplace=True)
        label = ds.run_name + ' (' + str(ni) + ')'
        # Plot
        h = ax.plot('time', var_name, data=df, label=label)
    else:
        for col_name in df.columns:
            if col_name == 'time':
                continue
            df.rename(columns={col_name: var_name + str(col_name)}, 
                      inplace=True)
            label = ds.run_name + ' (' + str(ni) + ', ' + str(col_name) + ')'
            # Plot
            h = ax.plot('time', var_name + str(col_name), data=df, label=label)

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

def plot_ice_var(df_ice, var_name, site, ax):
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
    if df['sem'].isna().all():
        h = ax.plot('time', 'mean', data=df, linestyle=':', label=label,
                    marker='o')
    else:
        h = ax.errorbar('time', 'mean', yerr='sem', data=df, 
                        capsize=5, linestyle=':', label=label)

    return h

def plot_handler(run_plot_dict, var_names, hist_dict, forc_var_map={},
                 ds_forc=None, ice_var_map={}, ice_sites=[], df_ice=None,
                 figsize=None, ax_font=14, lfont=10, xlim=None):
    """
    Handler function for plotting different runs and variables

    Parameters
    ----------
    run_plot_dict : dict
        Dictionary where the keys are the names of the runs to plot and value
        is a list of the cells (ni) to plot
    var_names : iterable
        Variable names to plot
    hist_dict : dict
        Dict containing the Icepack output, keyed on run_name
    forc_var_map : dict, optional
        Dict keyed on variable name and values are lists of forcing variable
        names to compare with. Default is {}.
    ds_forc : xarray Dataset, optional
        MDF formatted forcing dataset. Default is None.
    
    Returns
    -------
    Matplotlib figure object

    """

    # Create figsize
    if figsize is None:
        figsize = (10, 3*len(var_names))
    # Create figure and axes objects
    f, axs = plt.subplots(len(var_names), 1, sharex=True, figsize=figsize)
    if len(var_names) == 1:
        axs = [axs]

    # Loop through each variable
    for var_name, ax in zip(var_names, axs):
        # And through each run
        for run_name, nis in run_plot_dict.items():
            # and the desired cell(s) in each run
            for ni in nis:
                _ = plot_hist_var(hist_dict[run_name], var_name, ni, ax)
        
        # Plot forcing
        if var_name in forc_var_map:
            for forc_var_name in forc_var_map[var_name]:
                _ = plot_forc_var(ds_forc, forc_var_name, ax)
        
        # Plot ice variables
        for site in ice_sites:
            if var_name in ice_var_map:
                for ice_var_name in ice_var_map[var_name]:
                    _ = plot_ice_var(df_ice, ice_var_name, site, ax)
            
        # Axis labels
        ax.set_ylabel(var_name, fontsize=ax_font)
        ax.grid()
        # Legend
        ax.legend(fontsize=lfont, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    
    # xlimits on last plot
    if xlim is not None:
        axs[-1].set_xlim(xlim)

    #plt.show()
    return f, axs

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mosaic_snow2_acc": None,
            "mosaic_snow2_pre": None,
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
    value['hin'] = value['vicen'] / value['aicen']

# Load forcing
forcing_path = os.path.join(ip_dirs_path, "input", "Icepack_data",
                            "forcing", "MOSAiC")
atm_filename = "MOSAiC_MODF_20191011-20201001_v0.2.nc"
ds_atm = xr.open_dataset(os.path.join(forcing_path, atm_filename))
# Surface fluxes have opposite sign in Icepack output
for var_name in ['rlu', 'wthv', 'hfss_ec', 'hfss_bulk', 'wqv', 'hfls_bulk']:
    ds_atm[var_name] = -1*ds_atm[var_name]
ds_atm['hus'] = ds_atm['hus']/1000 # fix specific humidity error
ocn_filename = 'MOSAiC_ocn_MDF_20191007-20200920.nc'
ds_ocn = xr.open_dataset(os.path.join(forcing_path, ocn_filename))
ds_forc = xr.merge([ds_atm, ds_ocn])
# Create temp above frz
ds_forc['sst_above_frz'] = ds_forc['tos'] - ds_forc['tosf']

# mapping from Icepack names to forcing names
forc_var_map = {'flw': ['rld'],
                'fsw': ['rsd'],
                'Tair': ['tas'],
                'Qa': ['hus'],
                'flwout': ['rlu'],
                'fsens': ['wthv', 'hfss_ec', 'hfss_bulk'],
                'flat': ['wqv', 'hfls_bulk'],
                'sst': ['tos'],
                'sss': ['so'],
                'Tf': ['tosf'],
                'sst_above_frz': ['sst_above_frz'],
                }

# Ice data hardcoded from Transects
df_ice = pd.DataFrame(data={'mean': 
                            np.array([0.73, 0.85, 1.27, 1.53, 1.68]),
                            'time': ['2019-11-06', '2019-12-06', 
                                     '2020-02-06', '2020-03-26',
                                     '2020-04-13']})
df_ice['time'] = pd.to_datetime(df_ice.time)
df_ice['site'] = 'NLoop_level'
df_ice['sem'] = np.NaN
df_ice.set_index(['site', 'time'], inplace=True)
df_ice.columns = pd.MultiIndex.from_product([['hi_gems'], df_ice.columns])



ice_var_map = {'vice': ['hi_gems'],
               }

# Explore plotting variables
run_plot_dict = {"mosaic_snow2_acc": [2],
                 "mosaic_snow2_pre": [2],
                }
var_names = ['aice', 'vice', 'vice_chg', 'vsno', 'apnd', 'hpnd', 'ipnd']
site_names = ['NLoop_level']

f = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice)

for k, value in hist_dict.items():
    print(k)
    print(value.sel(ni=2, time="2020-04-30")['vice_chg'].values[0])

# Figures for talk
run_plot_dict = {"mosaic_snow2_acc": [2],
                 "mosaic_snow2_pre": [2],
                }
var_names = ['vice']
site_names = ['NLoop_level']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice)
axs[0].set_ylabel("Ice thickness (m)")
axs[0].set_xlim([datetime.datetime.fromisoformat('2019-11-06'),
            datetime.datetime.fromisoformat('2020-04-30')])
axs[0].set_ylim([0.7, 1.8])
axs[0].legend(labels=['Accumulation', 'Precipitation', 'Observations'])
plt.show()
f.savefig('./figures/Snow2_growth.png', bbox_inches='tight')