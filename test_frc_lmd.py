import os
import glob
import numpy as np
import netCDF4
from datetime import datetime

def make_frc(rootdir=os.getcwd(), frc_name='cp_frc_lmd.nc'):
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)

    path = './input/*his*.nc'
    hisfiles = sorted(glob.glob(path))
    assert len(hisfiles) > 0, "ERROR: No history files found in"

    # 1. Use the first file to initialize dimensions, variables, and global attributes
    with netCDF4.Dataset(hisfiles[0], 'r') as first_his:
        frc = netCDF4.Dataset(os.path.join(rootdir, frc_name), 'w', format='NETCDF4')
        frc.Description = 'Forcing for LMD offline simulation'
        frc.Author = 'cmagoo'
        frc.Created = datetime.now().isoformat()
        frc.type = 'ROMS FRC file'

        vars_list = ['swrad', 'shflux', 'sustr', 'svstr']
        dims = ['eta_rho', 'xi_rho', 'eta_u', 'xi_u', 'eta_v', 'xi_v']

        # Clone spatial dimensions
        for dim in dims:
            frc.createDimension(dim, len(first_his.dimensions[dim]))

        # Calculate total ocean_time steps across ALL 300 files
        total_nt = 0
        for f in hisfiles:
            with netCDF4.Dataset(f, 'r') as temp_his:
                total_nt += len(temp_his.dimensions['ocean_time'])

        # Create time dimension and variable metadata
        frc.createDimension('ocean_time', total_nt)
        frc.createVariable('ocean_time', 'd', ('ocean_time',))
        frc['ocean_time'].long_name = first_his['ocean_time'].long_name
        frc['ocean_time'].units = first_his['ocean_time'].units
        frc['ocean_time'].calendar = first_his['ocean_time'].calendar
        frc['ocean_time'].field = "time, scalar, series"

        # Initialize the forcing variables with metadata from the first file
        for var in vars_list:
            frc.createVariable(var, 'f8', first_his[var].dimensions)
            frc[var].long_name = first_his[var].long_name
            frc[var].units = first_his[var].units
            frc[var].field = first_his[var].field
            frc[var].time = first_his[var].time

    # 2. Loop through all files sequentially to stream and append data blocks
    time_idx = 0
    for i, f in enumerate(hisfiles):
        print(f"Processing forcing file {i+1}/{len(hisfiles)}: {os.path.basename(f)}")
        with netCDF4.Dataset(f, 'r') as his:
            n_times = len(his.dimensions['ocean_time'])
            end_idx = time_idx + n_times

            # Append time slice
            frc['ocean_time'][time_idx:end_idx] = his['ocean_time'][:]

            # Append data slices for forcing fields
            for var in vars_list:
                frc[var][time_idx:end_idx, :, :] = his[var][:]

            time_idx = end_idx

    # 3. Finalize metadata attributes and close the file safely
    frc['shflux'].negative_value = "upward flux, cooling"
    frc['shflux'].positive_value = "downward flux, heating"
    frc['swrad'].negative_value = "upward flux, cooling"
    frc['swrad'].positive_value = "downward flux, heating"

    frc.close()
    print(f"\nSUCCESS: Generated forcing file -> {frc_name}")

if __name__ == '__main__':
    make_frc(rootdir=os.getcwd())
