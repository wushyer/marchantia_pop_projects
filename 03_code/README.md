# Population genomics of ancient and modern *Trichuris trichiura*

## Overview of workflow and code used in the analysis
The code used in this project has been divided into workbooks based on the analysis performed, which are listed below. Each workbook can be accessed using the hyperlinks.

### Workbook:[Genome analyses](MP.00_genome_analyses.md)
- Reference genome
- Genome completeness using BUSCO
- Transfer of gene models using LiftOff
- Running interproscan on liftoff annotation

### Workbook:[Mapping](MP.01_mapping.md)
- Project setup
- Prepare reference genome
- Raw sequence data
- metadata
- trimming of the raw reads
- FIX: merge duplicate ancient read sets
- Mapping

### Workbook:[Kraken](MP.02_kraken_contmaination.md)
- kraken - modern samples
- kraken - ancient samples

### Workbook:[DNA damage](MP.03_dna_damage.md)
- pmdtools
- plot deamination frequencies
- trim bases from reads in bam to remove deamination

### Workbook:[Genome coverage](MP.04_genome_coverage.md)
- Genome wide coverage
- Generate quantitative stats on coverage for supplementary tables
- generate some coverage plots
- Genome-wide coverage to determine worm sex

### Workbook:[Variant calling and filtering](MP.05_variant_calling_and_filtering.md)
- Genome scope to estimate heterozyosity
- GATK
- Filter the VCF - SNPable
- Filter the VCF - hardfilter
- Querying SNP and INDEL QC profiles to determine thresholds for filters
- Applying filters to the variants
- merge VCFs
- Filter genotypes based on depth per genotype
- Sample missingness
- Generate an ALL SITES variant set for running pixy

### Workbook:[Sampling sites](MP.06_sampling_site_maps_and_data.md)
- World map
- Sampling timepoints

### Workbook:[PCAs](MP.07_PCAs.md)
- PCA of mitochondrial variants
- PCA of nuclear variants

### Workbook:[ANGSD](MP.08_ANGSD.md)
- IBS and Coviance matrices
- Single allele analyses of mitochondrial and nuclear variants

### Workbook:[NGSadmix](MP.09_NGSadmix.md)
- Admixture plots
- Clumpak to determine optimal K

### Workbook:[Treemix](MP.10_treemix.md)
- Treemix
- plotting treemix data
- Estimating the optimal number of migration edges
- calculate the variance explained by the data

### Workbook:[Admixtools](MP.11_admixtools.md)
- Prepared data and run admixtools
- plotting admix data

### Workbook:[Population history using SMC++](MP.12_smc++.md)
- run smc++
- plot all smcpp plots together

### Workbook:[Genome-wide genetic variation](MP.13_genomewide_genetic_variation.md)
- Running pixy to calculate nucleotide diversity, dXY and Fst between groups
- analyses of nucleotide diversity (Pi)
- dXY and Fst
- extracting top X% of Fst values for each comparison
- Private and shared variation between populations
- Relatedness and kinship between samples in a population
- genome-wide analyses presented in the original draft of the manuscript

