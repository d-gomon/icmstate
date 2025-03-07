#!/bin/bash


#Scenario 2

## MSM (Time homogeneous)

sbatch --job-name="sc6n500obs4N1000msm" --export=scenario=6,n=500,n_obs=4,N=1000,method="msm",RNG=1 --time=2:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm

sbatch --job-name="sc6n500obs10N1000msm" --export=scenario=6,n=500,n_obs=10,N=1000,method="msm",RNG=1 --time=2:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-short MSMsimulation.slurm



## multinomial EM

sbatch --job-name="sc6n500obs4N1000mult" --export=scenario=6,n=500,n_obs=4,N=1000,method="multinomial",RNG=1 --time=2-15:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-long MSMsimulation.slurm

sbatch --job-name="sc6n500obs10N1000mult" --export=scenario=6,n=500,n_obs=10,N=1000,method="multinomial",RNG=1 --time=2-15:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-long MSMsimulation.slurm

## Poisson EM


sbatch --job-name="sc6n500obs4N1000pois" --export=scenario=6,n=500,n_obs=4,N=1000,method="poisson",RNG=1 --time=5-15:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-long MSMsimulation.slurm

sbatch --job-name="sc6n500obs10N1000pois" --export=scenario=6,n=500,n_obs=10,N=1000,method="poisson",RNG=1 --time=5-15:00:00 \
--ntasks=1 --cpus-per-task=12 --mem-per-cpu=10000M --partition=cpu-long MSMsimulation.slurm




