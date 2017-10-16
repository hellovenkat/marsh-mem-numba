#!/bin/sh

#SBATCH --job-name=numba_12cores
#SBATCH --output=job%j.out
#SBATCH --error job%j.err
#SBATCH -p all 
###Number of Cores Max 20
#SBATCH -n 12
###Number of Nodes
#SBATCH -N 1

source /share/apps/Modules/3.2.10/init/sh
module load python/anaconda4.3

#echo "test"
#date
#hostname
#echo "done"
#echo "Here u go"
python MEM_testing_numba_DEM5.py

