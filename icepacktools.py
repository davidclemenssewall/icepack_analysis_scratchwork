"""Tools for examining Icepack output in python

@author: David Clemens-Sewall, NSF NCAR
"""

import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd

# Function for reading history
def load_icepack_hist(run_name, icepack_dirs_path, hist_filename=None,
                      sst_above_frz=True, volp=False,
                      trcr_dict=None, trcrn_dict=None):
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
    volp : bool, optional
        Whether or not to compute the pond volume per grid cell area (units m).
        requires alvl and alvln tracers. Default is False.
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

    # Add pond volume per unit area
    if volp:
        if ('alvl' in ds.data_vars) and ('alvln' in ds.data_vars) and (
            'apnd' in ds.data_vars) and ('apndn' in ds.data_vars) and (
            'hpnd' in ds.data_vars) and ('hpndn' in ds.data_vars):
            ds['volp'] = ds['aice']*ds['alvl']*ds['apnd']*ds['hpnd']
            ds['volpn'] = ds['aicen']*ds['alvln']*ds['apndn']*ds['hpndn']
        else:
            raise RuntimeError("missing data variables needed for volp(n)")

    # Add the run name as an attribute
    ds.attrs.update({'run_name': run_name})

    return ds

# Function for plotting single Icepack output
def plot_hist_var(ds, var_name, ni, ax, resample=None, cumulative=False,
                  mult=1):
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
        Whether the variable should be cumulative, useful for fluxes. The
        default is False.
    mult : float, optional
        Multiplier for values. Useful with cumulative to get the flux into
        correct units. The default is 1.

    Returns
    -------
    handle for matplotlib plot object

    """

    # Get variable as Pandas DataFrame with time as a column
    df = ds[var_name].sel(ni=ni).to_pandas()
    if resample:
        df = df.resample(resample).mean()
    if cumulative:
        df = df.cumsum()
    df *= mult
    df = df.reset_index()

    # Display
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
            h = ax.plot(df['time'], df[col_name], label=label)

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
        A dataframe where the row index is site and date and the
        column index is var_name, and then mean and standard error of mean
    var_name : str
    site : str
    ax :
    mean_only : bool, optional
        Whether to just display mean of ice variable and not st. dev. The 
        defaults is False.

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

def plot_handler(run_plot_dict, var_names, hist_dict, forc_var_map={},
                 ds_forc=None, ice_var_map={}, ice_sites=[], df_ice=None,
                 figsize=None, ax_font=14, lfont=10, xlim=None,
                 mean_only=False, resample=None, cumulative=False,
                 mult=1):
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
    ice_var_map : dict, optional
        Dict keyed on variable name and values are lists of ice df variable
        names to compare with. Default is {}.
    df_ice : Pandas dataframe, optional
        A dataframe where the row index is site and date and the
        column index is var_name, and then mean and standard error of mean
    mean_only : bool, optional
        Whether to just display mean of ice variable and not st. dev. The 
        defaults is False.
    resample : str, optional
        If provided, frequency string for DataFrame.resample(). If None do not
        resample. The default is None.
    cumulative : bool, optional
        Whether the variable should be cumulative, useful for fluxes. The
        default is False.
    mult : float, optional
        Multiplier for values. Useful with cumulative to get the flux into
        correct units. The default is 1.
    
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
                _ = plot_hist_var(hist_dict[run_name], var_name, ni, ax,
                                  resample=resample, cumulative=cumulative,
                                  mult=mult)
        
        # Plot forcing
        if var_name in forc_var_map:
            for forc_var_name in forc_var_map[var_name]:
                _ = plot_forc_var(ds_forc, forc_var_name, ax)
        
        # Plot ice variables
        for site in ice_sites:
            if var_name in ice_var_map:
                for ice_var_name in ice_var_map[var_name]:
                    _ = plot_ice_var(df_ice, ice_var_name, site, ax, 
                                     mean_only=mean_only)
            
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