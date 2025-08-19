#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=m
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.Europe

module load gsl/2.6-gcc-8.3.0
#mkdir European_smc
#while read CHROMOSOME; do smc++ vcf2smc European.recode.vcf.gz European_smc/European.${CHROMOSOME}.smc.gz ${CHROMOSOME} European:bgld1,Bonn2,Bonn3,Bonn6,Cam1,Cam2,Croatia,Czech1,Denmark1,Ehrenhausen3,Ehrenhausen6,France1,Germany1,Grafenegg1,Grafenegg2,Ireland10,Ireland12,Ireland13,Ireland18,Ireland19,Ireland1,Ireland21,Ireland22,Ireland2,Ireland36,Ireland37,Ireland38,Ireland3,Ireland5,Ireland6,Ireland7,MaG1,MiBa1,MiBa2,Norwich1,Norwich2,Norwich4,Oxford12,Oxford15,Oxford4,Oxford6,Schubert13,Schubert2,Schubert4,Schubert5,Schubert6,Schubert8,Schubert9,Sopron1,Zurich14,Zurich6,Zurich9; done < chromosomes.list
smc++ estimate --timepoints 0 10000000000 -o European_result/ 2.5e-9 European_smc/European.*.smc.gz
	
