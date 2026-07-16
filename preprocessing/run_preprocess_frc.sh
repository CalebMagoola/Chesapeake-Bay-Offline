#!/bin/bash
#SBATCH -J run_preprocess_frc
#SBATCH --mail-user=caleb.magoola@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -o ./preprocessing/preprocess_frc.log
#SBATCH --partition=hercules
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00            # Safe 2 hours allocation

# Go to the folder you submitted from (Chesapeake-Bay-Offline)
cd $SLURM_SUBMIT_DIR
python3 make_frc_cp.py
