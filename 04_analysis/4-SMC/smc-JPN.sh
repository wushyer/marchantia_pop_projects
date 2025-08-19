#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=c
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.asian

module load gsl/2.6-gcc-8.3.0
mkdir Japanese_smc
while read CHROMOSOME; do smc++ vcf2smc Japanese.recode.vcf.gz Japanese_smc/Japanese.${CHROMOSOME}.smc.gz ${CHROMOSOME} Japanese:Japan10,Japan12,Japan13,Japan16,Japan17,Japan18,Japan19,Japan1,Japan20,Japan21,Japan22,Japan26,Japan27,Japan28,Japan29,Japan2,Japan30,Japan4,Japan5,Japan6,Japan7,Japan9,Tak1F7,Tak2F7; done < chromosomes.list
smc++ estimate --timepoints 0 1000000000 -o Japanese_result/ 2.5e-9 Japanese_smc/Japanese.*.smc.gz	
