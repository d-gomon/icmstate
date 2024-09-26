#!/bin/bash


#Scenario 4

## MSM (Time homogeneous) - Too slow to run exact observations...

sbatch --job-name="sc4n100obs6N1000msm" --export=scenario=4,n=100,n_obs=6,N=1000,method="msm",RNG=1 --time=2-1:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-long MSMsimulation.slurm

sbatch --job-name="sc4n300obs6N1000msm" --export=scenario=4,n=300,n_obs=6,N=1000,method="msm",RNG=1 --time=2-1:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-long MSMsimulation.slurm

sbatch --job-name="sc4n500obs6N500msm" --export=scenario=4,n=500,n_obs=6,N=500,method="msm",RNG=1 --time=2-1:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-long MSMsimulation.slurm

sbatch --job-name="sc4n500obs6N500msm2" --export=scenario=4,n=500,n_obs=6,N=500,method="msm",RNG=2 --time=2-1:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-long MSMsimulation.slurm


## multinomial EM

sbatch --job-name="sc4n100obs6N1000bin" --export=scenario=4,n=100,n_obs=6,N=1000,method="multinomial",RNG=1 --time=2:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm

sbatch --job-name="sc4n300obs6N1000bin" --export=scenario=4,n=300,n_obs=6,N=1000,method="multinomial",RNG=1 --time=20:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

sbatch --job-name="sc4n500obs6N500bin" --export=scenario=4,n=500,n_obs=6,N=500,method="multinomial",RNG=1 --time=23:50:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

sbatch --job-name="sc4n500obs6N500bin2" --export=scenario=4,n=500,n_obs=6,N=500,method="multinomial",RNG=2 --time=23:50:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

## Poisson EM - DOES NOT WORK WITH EXACT OBSERVATIONS

#sbatch --job-name="sc4n100obs6N1000pois" --export=scenario=4,n=100,n_obs=6,N=1000,method="poisson" --time=2:00:00 \
#--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm

#sbatch --job-name="sc4n300obs6N1000pois" --export=scenario=4,n=300,n_obs=6,N=1000,method="poisson" --time=20:00:00 \
#--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

#sbatch --job-name="sc4n500obs6N500pois" --export=scenario=4,n=500,n_obs=6,N=500,method="poisson" --time=23:00:00 \
#--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm




