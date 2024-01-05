"""
Tools for comparing Icepack runs with each other and observations

"""

import xarray as xr
import os
import matplotlib.pyplot as plt
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

# Load history output
ip_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_dict = {"mosaic_raphael_syi": 'icepack.h.20191129.nc',
            "mosaic_raphael_fyi": 'icepack.h.20191129.nc',
            "sheba_raphael_fyi": None,
            "sheba_raphael_myi": None,
            "mosaic_raphael_L2": None,
            "mosaic_raphael_L2_pndaspect": None,
            "sheba_raphael_myi_nov28": None,
            "mosaic_raphael_fyi_rfracaice": None,
            "mosaic_raphael_fyi_nopndfbd": None,
            "mosaic_raphael_fyi_nopndfbd_pndaspect4": None,
            "mosaic_raphael_fyi_topo": None,
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

# Explore plotting variables
run_plot_dict = {"mosaic_raphael_syi": [2],
                 "mosaic_raphael_fyi": [2],
                 "sheba_raphael_fyi": [2],
                 "sheba_raphael_myi": [2],
                 "mosaic_raphael_L2": [2],
                 "sheba_raphael_myi_nov28": [2]
                }
var_names = ['aice', 'vice', 'vice_chg', 'vsno']

f = plot_handler(run_plot_dict, var_names, hist_dict)

var_names = ['apnd', 'hpnd', 'ipnd']

f = plot_handler(run_plot_dict, var_names, hist_dict)

# Test fluxes
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['vice', 'vice_chg']
f = plot_handler(run_plot_dict, var_names, hist_dict)
var_names = ['congel', 'frazil', 'snoice', 'meltt', 'meltl', 'meltb']
f = plot_handler(run_plot_dict, var_names, hist_dict, 
                 cumulative=True, resample='D', mult=24)



# Compare with changing pond parameterizations
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "mosaic_raphael_fyi_nopndfbd": [2],
                 "mosaic_raphael_fyi_nopndfbd_pndaspect4": [2],
                 "mosaic_raphael_fyi_rfracaice": [2],
                 "mosaic_raphael_fyi_topo": [2],
                 }
var_names = ['vice', 'apnd', 'hpnd', 'ipnd']
ice_var_map = {'vice': ['hi_hotwire'],
               'vsno': ['hs_gauge'],
               'hpnd': ['hpnd'],
               'apnd': ['apnd'],
               }
site_names = ['Reunion Stakes/dart_stakes_clu_12']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice, mean_only=True)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-07-28')])
axs[0].set_ylabel('Ice Thickness (m)')
axs[1].set_ylabel('Pond Area Fraction')
axs[2].set_ylabel('Pond Depth (m)')
plt.show()

# Look at per category variables
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "mosaic_raphael_fyi_topo": [2],
                 }
var_names = ['vicen','aicen','apndn','hpndn','ipndn']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-07-28')])
plt.show()

# Look at per category variables
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "mosaic_raphael_fyi_nopndfbd_pndaspect4": [2],
                 }
var_names = ['vicen','aicen','apndn','hpndn','ipndn']
f, axs = plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-07-28')])
plt.show()


# Pond aspect
run_plot_dict = {"mosaic_raphael_L2": [2],
                 "mosaic_raphael_L2_pndaspect": [2],
                }
var_names = ['aice', 'vice', 'vsno', 'apnd', 'hpnd', 'ipnd', 'pndaspect']

f, axs = plot_handler(run_plot_dict, var_names, hist_dict)
axs[-1].set_ylim([0, 1.0])
axs[-1].set_xlim([datetime.datetime.fromisoformat('2020-05-20'),
                  datetime.datetime.fromisoformat('2020-08-07')])
plt.show()

# Explore plotting with forcing
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['flwout', 'fsens', 'flat']
f = plot_handler(run_plot_dict, var_names, hist_dict, 
                 forc_var_map=forc_var_map, ds_forc=ds_forc)
# Explore plotting with forcing
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['sst', 'Tf', 'sst_above_frz']
f = plot_handler(run_plot_dict, var_names, hist_dict, 
                 forc_var_map=forc_var_map, ds_forc=ds_forc)

# Explore plotting fyi with ice state
run_plot_dict = {"mosaic_raphael_fyi": [2]}
var_names = ['vice', 'vsno', 'apnd', 'hpnd']
site_names = ['Stakes 1/dart_stakes_clu_4',
              'Ridge Ranch/dart_stakes_clu_6',
              'Runaway Stakes/dart_stakes_clu_7',
              'Drone Bones/dart_stakes_clu_11',
              'Reunion Stakes/dart_stakes_clu_12']
f = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice)

# Mosaic SYI with ice state
run_plot_dict = {"mosaic_raphael_syi": [2]}
var_names = ['vice', 'vsno', 'apnd', 'hpnd']
site_names = ['Bow Stakes/dart_stakes_clu_1',
              'Stakes 3/dart_stakes_clu_3',
              'MET Stakes/dart_stakes_clu_5',
              'Beanpole Stakes/dart_stakes_clu_13']
f = plot_handler(run_plot_dict, var_names, hist_dict, ice_var_map=ice_var_map,
                 ice_sites=site_names, df_ice=df_ice)


f, ax = plt.subplots(1,1)
plot_ice_var(df_ice, 'hpnd', 'Ridge Ranch/dart_stakes_clu_6', ax)

for k, value in hist_dict.items():
    print(k)
    print(value.sel(ni=2, time="2020-05-12")['vice_chg'].values[0])

# Explore plotting variables
run_plot_dict = {"mosaic_raphael_fyi": [2],
                 "sheba_raphael_myi": [2],
                 }
var_names = ['aice', 'vice', 'apnd', 'hpnd', 'apndn', 'hpndn']

f = plot_handler(run_plot_dict, var_names, hist_dict)
