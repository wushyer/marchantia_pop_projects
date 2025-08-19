#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=c
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.pixy3

module load tabix/0.2.6-gcccore-7.3.0
pixy --stats pi fst dxy --vcf merge.pure.pixy.vcf.gz --populations sample.info --window_size 20000  --n_cores 20 --output_prefix all_pixy_20K  --bypass_invariant_check yes
pixy --stats pi fst dxy --vcf merge.pure.pixy.vcf.gz --populations sample.info --window_size 500  --n_cores 20 --output_prefix all_pixy_500  --bypass_invariant_check yes
