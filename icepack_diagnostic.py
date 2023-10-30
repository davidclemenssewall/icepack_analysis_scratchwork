"""
Tools for comparing Icepack runs with each other and observations

"""

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd

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
    df.rename(columns={0: var_name}, inplace=True)
    label = ds.run_name + ' (' + str(ni) + ')'
    # Plot
    h = ax.plot('time', var_name, data=df, label=label)

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

def plot_handler(run_plot_dict, var_names, hist_dict, forc_var_map={},
                 ds_forc=None,
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
        
        # Axis labels
        ax.set_ylabel(var_name, fontsize=ax_font)
        ax.grid()
    
    # Legend on last plot
    axs[-1].legend(fontsize=lfont)
    # xlimits on last plot
    if xlim is not None:
        axs[-1].set_xlim(xlim)

    plt.show()
    return f
    

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mosaic_raphael_syi": 'icepack.h.20191129.nc',
            "mosaic_raphael_fyi": 'icepack.h.20191129.nc'
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


# Explore plotting variables
run_plot_dict = {"mosaic_raphael_syi": [1, 2],
                 "mosaic_raphael_fyi": [1, 2],
                }
var_names = ['vice',
             'vsno',
             'apnd',
             'hpnd',
             'ipnd',
             ]

f = plot_handler(run_plot_dict, var_names, hist_dict)

# Explore plotting with forcing
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['flwout', 'fsens', 'flat', 'sst', 'Tf', 'sst_above_frz']
f = plot_handler(run_plot_dict, var_names, hist_dict, 
                 forc_var_map=forc_var_map, ds_forc=ds_forc)



