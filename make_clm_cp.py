import numpy as np
import netCDF4
from datetime import datetime
import os
import glob

"""
This file prepares online ROMS outputs to be used as offline inputs (climatology files).
Uses direct sequential file processing with netCDF4 to bypass slow multi-file indexing.
"""

def make_clm(rootdir=os.getcwd(), clm_name='ches_clm_modern.nc'):
    assert os.path.exists(rootdir), ('%s does not exist.' % rootdir)

    # 1. Look for history files in the "history" folder
    path = './history/*his*.nc'
    hisfiles = sorted(glob.glob(path))
    
    if not hisfiles:
        raise FileNotFoundError(f"No history files found matching pattern: {path}")

    num_files = len(hisfiles)
    print(f"Found {num_files} history files in './history/' to process.")
    
    # Setup the output folder and dataset
    output_dir = os.path.join(rootdir, 'input')
    os.makedirs(output_dir, exist_ok=True) 
    output_path = os.path.join(output_dir, clm_name)

    # Core physical and BGC variables to keep
    keep = ['zeta', 'ubar', 'vbar', 'u', 'v', 'omega', 'temp', 'salt', 'NO3', 'NH4', 'chlorophyll', 'PO4', 'oxygen']
    tnames = ['zeta_time', 'u2d_time', 'v2d_time', 'u3d_time', 'v3d_time', 'ocean_time', 'temp_time', 'salt_time', 
              'ocean_time', 'ocean_time', 'ocean_time', 'ocean_time', 'ocean_time']

    # Map our input-to-output variable names (handling 'oxygen' -> 'oxyg')
    var_map = {v: ('oxyg' if v == 'oxygen' else v) for v in keep}

    # Open the first file to copy the structure and dimensions
    print("Reading structure from the first history file...")
    with netCDF4.Dataset(hisfiles[0], 'r') as first_nc:
        
        # Determine total ocean_time steps across all files
        total_time_steps = 0
        for f in hisfiles:
            with netCDF4.Dataset(f, 'r') as temp_nc:
                total_time_steps += len(temp_nc.dimensions['ocean_time'])
        print(f"Total time steps to concatenate across all files: {total_time_steps}")

        clm = netCDF4.Dataset(output_path, 'w', format='NETCDF4')
        clm.Description = 'Climatology for modern offline simulation'
        clm.Author = 'Gemini & Collab'
        clm.Created = datetime.now().isoformat()
        clm.type = 'ROMS CLM file'

        # Copy grid dimensions
        dims = ['xi_rho', 'xi_u', 'xi_v', 'eta_rho', 'eta_u', 'eta_v', 's_rho', 's_w']
        for dim in dims:
            if dim in first_nc.dimensions:
                clm.createDimension(dim, len(first_nc.dimensions[dim]))

        # Define time dimensions (unlimited) and coordinate variables
        clm.createDimension('ocean_time', None)
        clm.createVariable('ocean_time', 'f8', ('ocean_time',))
        
        # Get raw calendar and units
        raw_units = first_nc['ocean_time'].units
        calendar = getattr(first_nc['ocean_time'], 'calendar', 'proleptic_gregorian')
        
        # Convert seconds units to days
        if "seconds" in raw_units:
            new_units = raw_units.replace("seconds", "days")
            convert_to_days = True
        else:
            new_units = raw_units
            convert_to_days = False

        clm['ocean_time'].long_name = "time since initialization"
        clm['ocean_time'].units = new_units
        clm['ocean_time'].calendar = calendar
        clm['ocean_time'].field = "time, scalar, series"

        # Create variable definitions and their custom time coordinates
        for var, tname in zip(keep, tnames):
            out_var_name = var_map[var]
            if var not in first_nc.variables:
                continue

            # Create the sub-time dimension and tracking coordinate variable if unique
            if tname not in clm.dimensions:
                clm.createDimension(tname, None)
                clm.createVariable(tname, 'f8', (tname,))
                clm[tname].long_name = f"time for {first_nc[var].long_name}"
                clm[tname].units = new_units
                clm[tname].calendar = calendar
                clm[tname].field = f"{tname}, scalar, series"

            # Re-map variable dimensions to use the designated tname
            in_dims = first_nc[var].dimensions
            out_dims = [tname if d == 'ocean_time' else d for d in in_dims]

            clm.createVariable(out_var_name, 'f8', tuple(out_dims))
            clm[out_var_name].long_name = first_nc[var].long_name
            clm[out_var_name].time = tname
            
            if 'units' in first_nc[var].ncattrs():
                clm[out_var_name].units = first_nc[var].units
            if 'field' in first_nc[var].ncattrs():
                clm[out_var_name].field = first_nc[var].field

    # 2. Sequentially read each file and append the data slices
    current_index = 0
    for idx, f in enumerate(hisfiles):
        print(f"[{idx+1}/{num_files}] Reading: {os.path.basename(f)}...")
        
        with netCDF4.Dataset(f, 'r') as src:
            n_steps = len(src.dimensions['ocean_time'])
            
            # Read and scale time
            src_time = src['ocean_time'][:]
            if convert_to_days:
                src_time = src_time / 86400.0
            
            # Append times
            clm['ocean_time'][current_index:current_index + n_steps] = src_time
            
            # Append variable data
            for var in keep:
                out_var_name = var_map[var]
                if var not in src.variables:
                    continue
                
                tname = clm[out_var_name].time
                # Write to the specific segment of the unlimited time dimension
                clm[tname][current_index:current_index + n_steps] = src_time
                clm[out_var_name][current_index:current_index + n_steps, ...] = src[var][:]
                
            current_index += n_steps

    clm.close()
    print(f"\nPre-processing successfully finished! Saved output to: {output_path}")

if __name__ == '__main__':
    make_clm(rootdir=os.getcwd())
