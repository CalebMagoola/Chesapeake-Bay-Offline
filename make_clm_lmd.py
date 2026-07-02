import numpy as np
import netCDF4
from datetime import datetime
import os
import sys
import glob

"""
This file helps prepare online outputs to be used as offline inputs (history files)
Creates a clm file, climatology for offline simulation, requires a set of history files
Origination: Julia Crespin
Modified: Caleb Magoola, July 2026
"""

def make_clm(rootdir=os.getcwd(), clm_name='cp_clm_lmd.nc'):
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)

    path= './input/*his*.nc'
    hisfiles= sorted(glob.glob(path))

    his = netCDF4.MFDataset(hisfiles)

    keep = ['zeta', 'ubar', 'vbar', 'u', 'v', 'omega', 'temp', 'salt']
    tnames = ['zeta_time', 'u2d_time', 'v2d_time', 'u3d_time', 'v3d_time',
                'ocean_time', 'temp_time', 'salt_time']

    clm = netCDF4.Dataset(os.path.join(rootdir, clm_name), 'w', format='NETCDF4')
    clm.Description = 'Climatology for offline simulation'
    clm.Author = 'cmagoo'
    clm.Created = datetime.now().isoformat()
    clm.type = 'ROMS CLM file'

    dims= ['xi_rho','xi_u','xi_v','eta_rho','eta_u','eta_v','s_rho','s_w']

    for dim in dims:
        clm.createDimension(dim, np.size(his.dimensions[dim]))

    ########################################

    nt = len(his.dimensions['ocean_time'])
    clm.createDimension('ocean_time', nt)
    clm.createVariable('ocean_time', 'd', ('ocean_time',))
    clm['ocean_time'][:] = his['ocean_time'][:]
    clm['ocean_time'].long_name = his['ocean_time'].long_name
    clm['ocean_time'].units = his['ocean_time'].units
    clm['ocean_time'].calendar =his['ocean_time'].calendar
    clm['ocean_time'].field = "time, scalar, series"

    for var, tname in zip(keep, tnames):
        if tname not in clm.dimensions.keys():
            clm.createDimension(tname, nt)
        if tname != 'ocean_time':
            clm.createVariable(tname, 'd', (tname,))
            clm[tname][:] = his['ocean_time'][:]
            clm[tname].long_name = "time for "+his[var].long_name
            clm[tname].units = his['ocean_time'].units
            clm[tname].calendar =his['ocean_time'].calendar
            clm[tname].field = tname +", scalar, series"

        clm.createVariable(var, 'f8', (his[var].dimensions))
        clm[var][:] = his[var][:]
        clm[var].long_name = his[var].long_name
        clm[var].time = tname
        try:
            clm[var].units = his[var].units
            clm[var].field = his[var].field
        except:
            pass
    """
    ak= ['AKt','AKv','AKs','tke','gls']
    for var in ak[:]:
        clm.createVariable(var, 'f8', (his[var].dimensions))
        clm[var][:] = his[var][:]
        clm[var].long_name = his[var].long_name
        clm[var].time = 'ocean_time'
        clm[var].location = 'face'
        clm[var].units = his[var].units
        clm[var].field = his[var].field
    """
    clm.close()
    his.close()

if __name__ == '__main__':

    make_clm(rootdir=os.getcwd())
