#!/bin/bash -l

#SBATCH -J jobMODEL_OFFLINE
#SBATCH --mail-user=calebmagoola@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -o ./offline_run.log
#SBATCH --partition=hercules
#SBATCH --chdir=.
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --mem=45G
#SBATCH --time=4:00:00

# Clear out login environment configurations 
module purge

# Load the exact high-performance libraries matching your online run
module load intel-oneapi-compilers/2022.2.1
module load intel-oneapi-mpi/2021.7.1
module load hdf5/1.12.2
module load netcdf-c/4.9.0
module load netcdf-fortran/4.6.0
module load zlib/1.2.13

# Load and activate your verified Python environment
# module load miniconda3/24.3.0
# source $(conda info --base)/etc/profile.d/conda.sh
# conda activate base

# Prevent memory array allocation stack overflows (Critical for 3D grids)
ulimit -s unlimited

# Execute using srun across your 24 parallel processors pointing to the correct binary location
srun -n $SLURM_NTASKS ./romsM external/cp_offline_lmd1.in
