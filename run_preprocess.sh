#!/bin/bash -l

#SBATCH -J roms_preprocess
#SBATCH --account=aoml-oap-oa2o
#SBATCH --mail-user=caleb.magoola@noaa.gov
#SBATCH --mail-type=FAIL             # Changed to FAIL so it only bugs you if it crashes
#SBATCH -o ./preprocess_%j.log       # Use %j so multiple test runs don't overwrite the same log
#SBATCH --partition=hercules         # Sticks to the Hercules compute partition
#SBATCH --chdir=.
#SBATCH --nodes=1
#SBATCH --ntasks=1                   # Changed to 1: Python loops run on a single core by default
#SBATCH --ntasks-per-node=1          # Changed to 1: No need to reserve 24 cores for a single Python thread
#SBATCH --mem=30G                    # 20G is plenty of headroom to hold the largest 3D arrays in memory
#SBATCH --time=03:00:00 


# Load the exact working environment we found earlier
module load miniconda3/24.3.0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate base

# Run the heavy climatology script safely on a compute node
echo "Starting Climatology Generation..."
python3 -u test_clm_lmd.py
echo "Preprocessing Finished!"
