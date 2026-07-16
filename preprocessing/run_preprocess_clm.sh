#!/bin/bash
#SBATCH -J run_preprocess_clm
#SBATCH --mail-user=calebmagoola@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -o ./preprocessing/preprocess_clm.log
#SBATCH --partition=hercules
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# Move up to the root directory where the python script and history folder live
cd $SLURM_SUBMIT_DIR

# Run the Python script
python3 make_clm_cp.py
