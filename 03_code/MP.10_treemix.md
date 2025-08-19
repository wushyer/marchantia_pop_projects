# Treemix


- https://ppp.readthedocs.io/en/latest/PPP_pages/examples.html
- https://speciationgenomics.github.io/Treemix/

## Contents
- Treemix
- plotting treemix data
- Estimating the optimal number of migration edges
- calculate the variance explained by the data

```bash
conda activate py-popgen

export LD_LIBRARY_PATH=/lustre/scratch118/infgen/team133/sd21/software/anaconda2/envs/py-popgen/lib/

mkdir /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/05_ANALYSIS/TREEMIX
cd /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/05_ANALYSIS/TREEMIX

vcftools \
--gzvcf Trichuris_trichiura.cohort.nuclear_variants.final.recode.vcf.gz \
--max-missing 1 \
--keep ../../04_VARIANTS/GATK_HC_MERGED/nuclear_3x_animal.list \
--bed ../../01_REF/trichuris_trichiura.autosomeLG.bed \
--recode --stdout | gzip > treemix.vcf.gz

./ldPruning.sh treemix.vcf.gz

# before pruning
vcftools --gzvcf treemix.vcf.gz
#> After filtering, kept 180443 out of a possible 356541 Sites

# after pruning
vcftools --gzvcf treemix.LDpruned.vcf.gz
#> After filtering, kept 53671 out of a possible 109195 Sites

bcftools query -l treemix.vcf.gz | awk '{split($1,pop,"."); print $1"\t"$1"\t"pop[2]}' > data.clust

# this "data.clust" needs to be edited to put the population in, eg.

"
AN_DNK_COG_EN_0012	AN_DNK_COG_EN_0012	ANCIENT
AN_NLD_KAM_EN_0034	AN_NLD_KAM_EN_0034	ANCIENT
MN_CHN_GUA_HS_001	MN_CHN_GUA_HS_001	CHN
MN_CHN_GUA_HS_002	MN_CHN_GUA_HS_002	CHN
MN_CHN_GUA_HS_003	MN_CHN_GUA_HS_003	CHN
"


cat data.clust | cut -f3 | sort | uniq > populations.list

./vcf2treemix.sh treemix.LDpruned.vcf.gz data.clust


# run treemix, across a range of migration edges, and across a range of seeds
#--- the range of seeds is used by optM below to estimate the optimal number of edges.

for i in {0..5}; do
     for j in {0..10}; do
          treemix -i treemix.LDpruned.treemix.frq.gz -seed $j -m $i -o treemix.m_$i.s_$j -root COLOBUS -bootstrap -k 500  > treemix_${i}_log &
     done;
done

```

## plotting treemix data

```R

# load libraries
library(RColorBrewer)
library(R.utils)

source("plotting_funcs.R")
prefix="treemix"

# plot trees across range of migration edges
par(mfrow=c(1,2))

pdf("treemix_edges_2-6.pdf")
for(edge in 0:5){
     plot_tree(cex=0.8,paste0(prefix,".m_",edge,".s_4"))
     title(paste(edge,"edges"))
}
dev.off()

for(edge in 0:5){
     pdf(paste0("plot_treemix_tree_m-",edge,".pdf"))
     plot_tree(cex=0.8,paste0(prefix,".m_",edge,".s_3"))
     plot_resid(stem=paste0(prefix,".m_",edge,".s_3"),pop_order="populations.list")
     title(paste(edge,"edges"))
     dev.off()
}


# plot residuals across range of migration edges - positive values suggest admixture
for(edge in 0:5){
     pdf(paste0("plot_treemix_residuals_m-",edge,".pdf"))
     plot_resid(stem=paste0(prefix,".m_",edge,".s_4"),pop_order="populations.list")
     title(paste(edge,"edges")
     dev.off()
}
```


## Estimating the optimal number of migration edges
- : https://rdrr.io/cran/OptM/#vignettes
```R
# run in the folder with the treemix output files
#install.packages("OptM")

library(OptM)

optM(folder="./")

#> The maximum value for delta m was 10.7624 at m = 3 edges.

# remake the plots, using 2 migration edges
prefix="treemix"

par(mfrow=c(1,1))

plot_tree(cex=0.8,paste0(prefix,".m_3.s_4"))
title(paste(1,"edges"))

plot_resid(stem=paste0(prefix,".",2),pop_order="populations.list")
title(paste(1,"edges"))

```

## calculate the variance explained by the data
- https://github.com/darencard/RADpipe/blob/master/treemixVarianceExplained.R

```bash

Rscript treemixVarianceExplained.R treemix.m_1.s_2

#> Standard error for all entries in the covariance matrix estimated from the data	0.000566063693877551
#> Variance of relatedness between populations explained by the model	0.997250751108934

```
