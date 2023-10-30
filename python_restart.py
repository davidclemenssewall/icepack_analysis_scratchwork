# python_restart.py
# Figure out how to create a netCDF restart file with an arbitrary snow and ice thickness distribution

import xarray as xr
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Start by loading and inspecting a restart file
icepack_dirs_path = "/home/dcsewall/code/docker_icepack_interactive/icepack-dirs"
run_name = "test_simplest_restart"
rest_path = os.path.join(icepack_dirs_path, "runs", run_name, "restart")
hist_filename = 'iced.2020-01-29-21600.nc'
ds_rest = xr.open_dataset(os.path.join(rest_path, hist_filename))

# Constants for grid
ni = 4
kcatbound = 1
nilyr = 7
nslyr = 1

if kcatbound == 1:
    ncat = 5
    cat_bounds = [0.0, 0.60, 1.40, 2.40, 3.60]
else:
    raise RuntimeError()

# For now we can just initialize a single column
# Create aicen, vicen, vsnon given thicknesses
aicen_data = np.array([0.1, 0.5, 0.25, 0.15, 0.1])
hicen_data = np.array([0.3, 1.0, 2.0, 3.0, 5.0])
hsnon_data = np.array([0.06, 0.1, 0.2, 0.3, 0.5])

aicen_in = np.zeros((ncat, ni))
aicen_in[:, 0] = aicen_data
vicen_in = np.zeros((ncat, ni))
vicen_in[:, 0] = aicen_data * hicen_data
vsnon_in = np.zeros((ncat, ni))
vsnon_in[:, 0] = aicen_data * hsnon_data

# Salinity (derive from cores?)
# Need this to set enthalpies
# for now start with Icepack default?
# BL99 prognostic salinity profile
# icepack_therm_shared.F90 icepack_init_thermo
saltmax = 3.2
min_salin = 0.1
nsal = 0.407
msal = 0.573
zn = (np.arange(1, nilyr+1) - 0.5) / nilyr # normalized depths
sprofile =(saltmax / 2) * (1-np.cos(np.pi * zn 
                                    **(nsal /(msal + zn))))
sprofile = np.fmax(sprofile, min_salin)
# that's probably completely unreasonable for MOSAiC but
# it's a start
sice_in = np.zeros((ncat, ni, nilyr))
for nc in np.arange(ncat):
    if aicen_data[nc] > 0:
        sice_in[nc, 0, :] = sprofile

# Air-snow and ice-ocean interface Temperatures
Tair_data = -20.0
sst_data = -1.8
Tsfcn_in = Tair_data * np.ones((ncat, ni))
sst_in = sst_data * np.ones(ni)

# Enthalpies (eventually get T from buoys)
# For now estimate steady-state profile from assumed conductivity
# q = k * dT/dz
# q = ks/hs * (Tair - Tint) = ki/hi * (Tint  - sst)
# ks/hs * Tair + ki/hi * sst = ks/hs * Tint + ki/hi * Tint
# = (ks/hs + ki/hi) * Tint
ksno = 0.3
kice = 2.03
Tint = (ksno * Tair_data / hsnon_data 
        + kice * sst_data / hicen_data
        ) / (ksno / hsnon_data + kice / hicen_data)
zn = (np.arange(1, nslyr+1) - 0.5) / nslyr
Tsno_data = Tair_data + ((Tint - Tair_data) * zn[:,np.newaxis]
                         * np.ones((1,ncat)))
zn = (np.arange(1, nilyr+1) - 0.5) / nilyr
Tice_data = Tint + ((sst_data - Tint) * zn[:,np.newaxis]
                         * np.ones((1,ncat)))
# Snow enthalpy
rhos = 330 # kg/m3 snow density
Lsub = 2.835e6
Lvap  = 2.501e6
Lfresh = Lsub - Lvap # J/Kg
cp_ice = 2106
qsnon_data = -rhos * (Lfresh - cp_ice * Tsno_data)

qsno_in = np.zeros((ncat, ni, nslyr))
for nc in np.arange(ncat):
    if aicen_data[nc] > 0:
        qsno_in[nc, 0, :] = qsnon_data

# Ice enthalpy
cp_ocn = 4218.0
rhow = 1026.0
rhoi = 917.0

def liquidus_brine_salinity_mush(zTin):
     
    bz1_liq = 0.0
    az1_liq = -18.48
    bz2_liq = 62.4
    az2_liq = -10.3085
    Tb_liq = -7.6362968855167352

     # temperature to brine salinity
    J1_liq = bz1_liq / az1_liq
    K1_liq = 0.001
    L1_liq = (1 + bz1_liq/1000) / az1_liq
    J2_liq = bz2_liq  / az2_liq
    K2_liq = 0.001
    L2_liq = (1 + bz2_liq/1000) / az2_liq

    t_high   = (zTin > Tb_liq)
    lsubzero = (zTin <= 0)

    Sbr = (((zTin + J1_liq) / (K1_liq * zTin + L1_liq)) * t_high +
        ((zTin + J2_liq) / (K2_liq * zTin + L2_liq)) * (1 - t_high))

    Sbr = Sbr * lsubzero

    return Sbr

SBr_data = liquidus_brine_salinity_mush(Tice_data)

phi_data = sprofile[:,np.newaxis] / np.fmax(SBr_data, 
                                            sprofile[:,np.newaxis])

zqin = (phi_data * (cp_ocn * rhow - cp_ice * rhoi) * Tice_data + 
        rhoi * cp_ice * Tice_data - (1.0 - phi_data) * rhoi * Lfresh)