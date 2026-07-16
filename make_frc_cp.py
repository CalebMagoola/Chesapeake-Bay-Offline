import numpy as np
import netCDF4
from datetime import datetime
import os
import glob

"""
This file helps prepare online outputs to be used as offline inputs (history files)
Creates a frc file, hydrodynamic forcing for offline simulation.
Uses direct sequential file processing with netCDF4 to bypass MFDataset constraints and memory bloat.
Origination: Julia Crespin
Modified: Caleb Magoola, July 2026
"""

def make_frc(rootdir=os.getcwd(), frc_name='cp_frc_lmd.nc'):
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)

    # Corrected path to point to the actual history folder
    path = './history/*his*.nc'
    hisfiles = sorted(glob.glob(path))
    
    if not hisfiles:
        raise FileNotFoundError(f"No history files found matching pattern: {path}")

    num_files = len(hisfiles)
    print(f"Found {num_files} history files to process for forcing data.")
    
    # Setup the output folder path
    output_dir = os.path.join(rootdir, 'input')
    os.makedirs(output_dir, exist_ok=True) 
    output_path = os.path.join(output_dir, frc_name)

    vars_to_keep = ['swrad', 'shflux', 'sustr', 'svstr']

    # Open the first file to clone the dimensions and structure
    print("Reading dimensions from the first history file...")
    with netCDF4.Dataset(hisfiles[0], 'r') as first_nc:
        
        # Calculate total steps dynamically across all files
        total_time_steps = 0
        for f in hisfiles:
            with netCDF4.Dataset(f, 'r') as temp_nc:
                total_time_steps += len(temp_nc.dimensions['ocean_time'])
        print(f"Total forcing time steps to concatenate: {total_time_steps}")

        frc = netCDF4.Dataset(output_path, 'w', format='NETCDF4')
        frc.Description = 'Forcing for LMD offline simulation'
        frc.Author = 'cmagoo & Gemini'
        frc.Created = datetime.now().isoformat()
        frc.type = 'ROMS FRC file'

        # Copy spatial grid dimensions
        dims = ['eta_rho', 'xi_rho', 'eta_u', 'xi_u', 'eta_v', 'xi_v']
        for dim in dims:
            if dim in first_nc.dimensions:
                frc.createDimension(dim, len(first_nc.dimensions[dim]))

        # Define unlimited time dimension
        frc.createDimension('ocean_time', None)
        frc.createVariable('ocean_time', 'f8', ('ocean_time',))
        
        # Clone time metadata attributes
        frc['ocean_time'].long_name = first_nc['ocean_time'].long_name
        frc['ocean_time'].units = first_nc['ocean_time'].units
        if 'calendar' in first_nc['ocean_time'].ncattrs():
            frc['ocean_time'].calendar = first_nc['ocean_time'].calendar
        frc['ocean_time'].field = "time, scalar, series"

        # Build definitions for the forcing variables
        for var in vars_to_keep:
            if var not in first_nc.variables:
                continue
            frc.createVariable(var, 'f8', first_nc[var].dimensions)
            frc[var].long_name = first_nc[var].long_name
            frc[var].units = first_nc[var].units
            frc[var].field = first_nc[var].field
            frc[var].time = first_nc[var].time

        # Set specific heat flux attributes
        frc['shflux'].negative_value = "upward flux, cooling"
        frc['shflux'].positive_value = "downward flux, heating"
        frc['swrad'].negative_value = "upward flux, cooling"
        frc['swrad'].positive_value = "downward flux, heating"

    # Sequentially read and append chunk arrays from each history file
    current_index = 0
    for idx, f in enumerate(hisfiles):
        print(f"[{idx+1}/{num_files}] Extracting forcing: {os.path.basename(f)}...", flush=True)
        
        with netCDF4.Dataset(f, 'r') as src:
            n_steps = len(src.dimensions['ocean_time'])
            
            # Append time coordinates
            frc['ocean_time'][current_index:current_index + n_steps] = src['ocean_time'][:]
            
            # Append 2D variable slices in bulk
            for var in vars_to_keep:
                if var not in src.variables:
                    continue
                frc[var][current_index:current_index + n_steps, ...] = src[var][:]
                
            current_index += n_steps

    frc.close()
    print(f"\nForcing file preprocessing successfully completed! Saved to: {output_path}")

if __name__ == '__main__':
    make_frc(rootdir=os.getcwd())
