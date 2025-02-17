#!/bin/bash


#Scenario 2

## MSM (Time homogeneous)

sbatch --job-name="sc5n100obs6N1000msm" --export=scenario=5,n=100,n_obs=6,N=1000,method="msm",RNG=1 --time=1:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm

sbatch --job-name="sc5n300obs6N1000msm" --export=scenario=5,n=300,n_obs=6,N=1000,method="msm",RNG=1 --time=1:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm

sbatch --job-name="sc5n500obs6N500msm" --export=scenario=5,n=500,n_obs=6,N=500,method="msm",RNG=1 --time=2:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm

sbatch --job-name="sc5n500obs6N500msm2" --export=scenario=5,n=500,n_obs=6,N=500,method="msm",RNG=2 --time=2:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm


## multinomial EM

sbatch --job-name="sc5n100obs6N1000bin" --export=scenario=5,n=100,n_obs=6,N=1000,method="multinomial",RNG=1 --time=2:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm

sbatch --job-name="sc5n300obs6N1000bin" --export=scenario=5,n=300,n_obs=6,N=1000,method="multinomial",RNG=1 --time=20:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

sbatch --job-name="sc5n500obs6N500bin" --export=scenario=5,n=500,n_obs=6,N=500,method="multinomial",RNG=1 --time=23:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

sbatch --job-name="sc5n500obs6N500bin2" --export=scenario=5,n=500,n_obs=6,N=500,method="multinomial",RNG=2 --time=23:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

## Poisson EM

sbatch --job-name="sc5n100obs6N1000pois" --export=scenario=5,n=100,n_obs=6,N=1000,method="poisson",RNG=1 --time=2:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm

sbatch --job-name="sc5n300obs6N1000pois" --export=scenario=5,n=300,n_obs=6,N=1000,method="poisson",RNG=1 --time=20:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

sbatch --job-name="sc5n500obs6N500pois" --export=scenario=5,n=500,n_obs=6,N=500,method="poisson",RNG=1 --time=23:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm

sbatch --job-name="sc5n500obs6N500pois2" --export=scenario=5,n=500,n_obs=6,N=500,method="poisson",RNG=2 --time=23:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-medium MSMsimulation.slurm




