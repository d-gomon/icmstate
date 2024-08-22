#!/bin/bash
#SBATCH --output=%x_%j.out
#SBATCH --mail-user="<d.gomon@math.leidenuniv.nl>"
#SBATCH --mail-type="ALL"

#
# loading the R module
#
module load ALICE/default
module load R/4.4.0-gfbf-2023a

# get the current working directory
CWD=$(pwd)
echo "I am currently in $CWD"
# get the current time and date
DATE=$(date)
echo "It is now $DATE"

# Run the file
Rscript General_Simulation.R $scenario $n $n_obs $N $method $RNG

echo "#### Finished script"