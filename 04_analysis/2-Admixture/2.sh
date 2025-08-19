#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=m
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.ad09


for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv plink_miss09.bed $K | tee 09.log${K}.out; done
