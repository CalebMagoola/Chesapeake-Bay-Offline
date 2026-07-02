import numpy as np
import netCDF4
from datetime import datetime
import os
import sys
import glob

"""
This file helps prepare online outputs to be used as offline inputs (history files)
Creates a frc file, hydrodynamic forcing for offline simulation, requires a set of history files
Origination: Julia Crespin
Modified: Caleb Magoola, July 2026
"""

def make_frc(rootdir=os.getcwd(), frc_name='cp_frc_lmd.nc'):
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)

    path= './input/*his*.nc'
    hisfiles= sorted(glob.glob(path))

    his = netCDF4.MFDataset(hisfiles)
    frc = netCDF4.Dataset(os.path.join(rootdir, frc_name), 'w', format='NETCDF4')
    frc.Description = 'Forcing for LMD offline simulation'
    frc.Author = 'cmagoo'
    frc.Created = datetime.now().isoformat()
    frc.type = 'ROMS FRC file'

    vars = ['swrad','shflux', 'sustr', 'svstr']
    dims= ['eta_rho','xi_rho','eta_u','xi_u','eta_v','xi_v']

    for dim in dims[:]:
        frc.createDimension(dim, (np.size(his.dimensions[dim])))

    nt = len(his.dimensions['ocean_time'])
    frc.createDimension('ocean_time', nt)
    frc.createVariable('ocean_time', 'd', ('ocean_time',))
    frc['ocean_time'][:]= his['ocean_time'][:]
    frc['ocean_time'].long_name= his['ocean_time'].long_name
    frc['ocean_time'].units= his['ocean_time'].units
    frc['ocean_time'].calendar= his['ocean_time'].calendar
    frc['ocean_time'].field= "time, scalar, series"

    for var in vars[:]:
        frc.createVariable(var, 'f8', (his[var].dimensions))
        frc[var][:] = his[var][:]
        frc[var].long_name = his[var].long_name
        frc[var].units = his[var].units
        frc[var].field = his[var].field
        frc[var].time = his[var].time

    frc['shflux'].negative_value = "upward flux, cooling"
    frc['shflux'].positive_value = "downward flux, heating"
    frc['swrad'].negative_value = "upward flux, cooling"
    frc['swrad'].positive_value = "downward flux, heating"

    frc.close()
    his.close()

if __name__ == '__main__':

    make_frc(rootdir=os.getcwd())
