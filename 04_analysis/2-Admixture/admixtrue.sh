#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=m
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.ad

#vcftools --vcf merged.nuclear_variants.final.snps.vcf --plink --out plink
#plink --noweb --file plink --geno 0.1 --maf 0.05 --make-bed --out plink

#vcftools --vcf merged.nuclear_variants.miss09.final.snps.vcf --plink --out plink_miss09
#plink --noweb --file plink_miss09 --geno 0.1 --maf 0.05 --make-bed --out plink_miss09

for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv plink.bed $K | tee log${K}.out; done
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv plink_miss09.bed $K | tee log${K}.out; done
