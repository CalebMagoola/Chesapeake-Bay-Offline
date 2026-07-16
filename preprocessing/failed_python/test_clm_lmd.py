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
    assert len(hisfiles) > 0, "ERROR: No history files found in ./input/"

    # 1. Use the first file to initialize dimensions and structure
    with netCDF4.Dataset(hisfiles[0], 'r') as first_his:
        clm = netCDF4.Dataset(os.path.join(rootdir, clm_name), 'w', format='NETCDF4')
        clm.Description = 'Climatology for offline simulation'
        clm.Author = 'cmagoo'
        clm.Created = datetime.now().isoformat()
        clm.type = 'ROMS CLM file'

        dims= ['xi_rho','xi_u','xi_v','eta_rho','eta_u','eta_v','s_rho','s_w']

        for dim in dims:
            clm.createDimension(dim, np.size(first_his.dimensions[dim]))

        # Calculate total ocean_time steps across ALL files
        total_nt = 0
        for f in hisfiles:
            with netCDF4.Dataset(f, 'r') as temp_his:
                total_nt += len(temp_his.dimensions['ocean_time'])

        # Create the standard ocean_time dimension
        clm.createDimension('ocean_time', total_nt)
        clm.createVariable('ocean_time', 'd', ('ocean_time',))
        clm['ocean_time'].long_name = first_his['ocean_time'].long_name
        clm['ocean_time'].units = first_his['ocean_time'].units
        clm['ocean_time'].calendar = first_his['ocean_time'].calendar
        clm['ocean_time'].field = "time, scalar, series"

        # Track the list of variables and their custom time names
        keep = ['zeta', 'ubar', 'vbar', 'u', 'v', 'omega', 'temp', 'salt']
        tnames = ['zeta_time', 'u2d_time', 'v2d_time', 'u3d_time', 'v3d_time',
                  'ocean_time', 'temp_time', 'salt_time']

        # Setup all dimensions, time coordinates, and variables metadata
        for var, tname in zip(keep, tnames):
            if tname not in clm.dimensions.keys():
                clm.createDimension(tname, total_nt)

            if tname != 'ocean_time':
                clm.createVariable(tname, 'd', (tname,))
                clm[tname].long_name = "time for " + first_his[var].long_name
                clm[tname].units = first_his['ocean_time'].units
                clm[tname].calendar = first_his['ocean_time'].calendar
                clm[tname].field = tname + ", scalar, series"

            clm.createVariable(var, 'f8', (first_his[var].dimensions))
            clm[var].long_name = first_his[var].long_name
            clm[var].time = tname
            try:
                clm[var].units = first_his[var].units
                clm[var].field = first_his[var].field
            except:
                pass

    # 2. Loop through all files sequentially to write data blocks
    time_idx = 0
    for i, f in enumerate(hisfiles):
        print(f"Processing climatology file {i+1}/{len(hisfiles)}: {os.path.basename(f)}")
        with netCDF4.Dataset(f, 'r') as his:
            n_times = len(his.dimensions['ocean_time'])
            end_idx = time_idx + n_times

            # Write times to all unique time dimensions
            clm['ocean_time'][time_idx:end_idx] = his['ocean_time'][:]
            for tname in tnames:
                if tname != 'ocean_time':
                    clm[tname][time_idx:end_idx] = his['ocean_time'][:]

            # Stream variable arrays dynamically based on shape (2D/3D/4D)
            for var in keep:
                if len(his[var].dimensions) == 4:    # e.g., (time, s_rho, eta, xi)
                    clm[var][time_idx:end_idx, :, :, :] = his[var][:]
                elif len(his[var].dimensions) == 3:  # e.g., (time, eta, xi)
                    clm[var][time_idx:end_idx, :, :] = his[var][:]

            time_idx = end_idx

    clm.close()
    print(f"\nSUCCESS: Generated climatology file -> {clm_name}")

if __name__ == '__main__':
    make_clm(rootdir=os.getcwd())
