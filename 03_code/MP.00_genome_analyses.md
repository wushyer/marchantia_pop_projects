# Genome analyses


## Contents:
- Reference genome
- Genome completeness using BUSCO
- Transfer of gene models using LiftOff
- Running interproscan on liftoff annotation



## Reference genome
- the reference is an unpublished Trichuris trichiura assembly generated in the Parasite Genomics team at Sanger
- this is a new PacBio assembly, rather than an improvement from the original assembly published by Foth et al (2014)
- this assembly is from worms extracted from Peter Nejsum that originated in Uganda, whereas the Foth assembly was from Ecuadorian worms.
- this assembly didn't have an annotation when the project started

### Genome stats
- genome stats of the new assembly
     - sum = 80573711, n = 113, ave = 713041.69, largest = 29164577
     - N50 = 11299416, n = 2
     - N60 = 9167782, n = 3
     - N70 = 5100676, n = 5
     - N80 = 2017282, n = 7
     - N90 = 643749, n = 14
     - N100 = 1808, n = 113
     - N_count = 250000
     - Gaps = 25




## Genome completeness using BUSCO
```bash

export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/pkgs/augustus-3.1-0/bin:$PATH"
export PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/pkgs/augustus-3.1-0/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/nfs/users/nfs_s/sd21/lustre118_link/software/anaconda2/pkgs/augustus-3.1-0/config"

bsub.py --threads 20 --queue long 20 busco_tt_new "busco -i trichuris_trichiura.fa -l /nfs/users/nfs_s/sd21/lustre118_link/databases/busco/metazoa_odb9 -o new_assembly_metazoa --mode genome --long -sp caenorhabditis -f --cpu 20"

bsub.py --threads 20 --queue long 20 busco_tt_old "busco -i trichuris_trichiura.PRJEB535.WBPS12.genomic.fa -l /nfs/users/nfs_s/sd21/lustre118_link/databases/busco/metazoa_odb9 -o old_assembly_metazoa --mode genome --long -sp caenorhabditis -f --cpu 20"

```

- Summarized benchmarking in BUSCO notation for file trichuris_trichiura.PRJEB535.WBPS12.genomic.fa
- BUSCO was run in mode: genome
     - C:81.9%[S:79.8%,D:2.1%],F:1.8%,M:16.3%,n:978
          - 801	Complete BUSCOs (C)
          - 780	Complete and single-copy BUSCOs (S)
          - 21	Complete and duplicated BUSCOs (D)
          - 18	Fragmented BUSCOs (F)
          - 159	Missing BUSCOs (M)
          - 978	Total BUSCO groups searched


- Summarized benchmarking in BUSCO notation for file trichuris_trichiura.fa
- BUSCO was run in mode: genome
     - C:81.4%[S:79.3%,D:2.1%],F:1.7%,M:16.9%,n:978
          - 797	Complete BUSCOs (C)
          - 776	Complete and single-copy BUSCOs (S)
          - 21	Complete and duplicated BUSCOs (D)
          - 17	Fragmented BUSCOs (F)
          - 164	Missing BUSCOs (M)
          - 978	Total BUSCO groups searched

- stats are slightly worse for the new assembly. Givne the new assembly incorporated PacBio data, I wonder to what degree it was polished....





## Transfer of gene models using LiftOff
- need to identify genes in new genome assembly, which does not yet have a stable annotations
- propose to simply lift over gene models from original annotation, simply for the purpose of asking, what genes are present
- doesnt have to be perfect
- using / trying the new "liftoff" tool - https://github.com/agshumate/Liftoff; https://doi.org/10.1093/bioinformatics/btaa1016


```bash

conda activate py37
pip3 install liftoff --user
conda install -c bioconda minimap2

wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/trichuris_trichiura/PRJEB535/trichuris_trichiura.PRJEB535.WBPS15.annotations.gff3.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/trichuris_trichiura/PRJEB535/trichuris_trichiura.PRJEB535.WBPS15.genomic.fa.gz

liftoff -g trichuris_trichiura.PRJEB535.WBPS15.annotations.gff3.gz trichuris_trichiura.fa trichuris_trichiura.PRJEB535.WBPS15.genomic.fa -o newversion.gff -u unmapped_features.txt


wc -l unmapped_features.txt
#> 1199 unmapped_features.txt

awk $3="gene" '{print}' newversion.gff
#> 8451 genes


# to be cross-compatible with Apollo, which has the same genome version but with updated scaffold names (changed by Faye after I had started the analysis), had to run the following on the genome and annotation

perl rename_scaffolds.pl trichuris_trichiura.fa > trichuris_trichiura.renamed.fa

perl rename_scaffolds.pl newversion.gff > newversion.renamed.gff

# where "rename_scaffolds.pl" is:

#!/usr/bin/env perl
use strict;
use warnings;
while (<>){
    s/Trichuris_trichiura/TTRE/;
    s/_00_/_unplaced_/;
    s/_1_/_chr1_/;
    s/_2_001/_chr2/;
    s/_3_/_chr3_/;
    s/_[0]+/_scaffold/;
    print $_;
}

# found that, becasue the original version has been annotated interpro descriptions etc, the liftover has misformatted some of the lines (in the same way that apollo dumps of the Haemonchus annotations were misformatted). To fix, ran the following:

cat newversion.renamed.gff | while read LINE; do
     if [[ $LINE =~ (^TTRE*|\#) ]]; then
          echo -ne "\n${LINE} ";
          else
          echo -n "${LINE} ";
     fi;
done > newversion.renamed.v2.gff3


cat newversion.gff | while read LINE; do
     if [[ $LINE =~ (^Trichuris*|\#) ]]; then
          echo -ne "\n${LINE} ";
          else
          echo -n "${LINE} ";
     fi;
done > newversion.v2.gff3

echo -e "##gff-version 3" > UPDATED_annotation.gff3
cat newversion.renamed.v2.gff3 | sort -k1,1 -k4,4n | sed '/^$/d' >> UPDATED_annotation.gff3

```

- ran really quickly, maybe a few minutes at most
- shows 1199 features in the "unmapped_features.txt" - have not idea how "bad" this is, but possible ok(?) coming from a draft assembly to something that is more complete. Very likely some at least are errors, but possibly not a big deal for my use.
- 8451 genes

- will attempt to map protein sequences back to genome with exonerate

```bash
# extract protein sequences
gffread <(zcat trichuris_trichiura.PRJEB535.WBPS15.annotations.gff3.gz) -g trichuris_trichiura.PRJEB535.WBPS15.genomic.fa -y old_proteins.fa
gffread <(zcat trichuris_trichiura.PRJEB535.WBPS15.annotations.gff3.gz) -g trichuris_trichiura.PRJEB535.WBPS15.genomic.fa -x old_cds.fa

# fix line lengths
fastaq to_fasta -l0 old_proteins.fa tmp; mv tmp old_proteins.fa
fastaq to_fasta -l0 old_cds.fa tmp; mv tmp old_cds.fa

# extract protein sequences for the unmapped genes
cat unmapped_features.txt | sed 's/'gene:'//g' | while read -r GENE; do
     grep -A1 ${GENE} old_proteins.fa >> unmapped_proteins.fa;
done

cat unmapped_features.txt | sed 's/'gene:'//g' | while read -r GENE; do
     grep -A1 ${GENE} old_cds.fa >> unmapped_cds.fa;
done

~sd21/bash_scripts/run_exonerate_splitter trichuris_trichiura.fa unmapped_proteins.fa

cat split_exonerate*out | Exonerate_to_evm_gff3.pl - > merged_exonerate.output

rm split_exonerate*out
rm x*
rm run_split*

```

### run diamond to query species
```bash
# load diamond
module load diamond/0.9.24--ha888412_1

diamond blastp \
--db /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/uniprot_ref_proteomes.dmnd \
--query unmapped_proteins.fa \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
--taxonnodes /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/nodes.dmp \
--tmpdir /dev/shm \
--taxonmap /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/prot.accession2taxid \
--max-target-seqs 1 > unmapped_proteins.diamond.out

diamond blastx \
--db /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/uniprot_ref_proteomes.dmnd \
--query unmapped_cds.fa \
--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
--taxonnodes /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/nodes.dmp \
--tmpdir /dev/shm \
--taxonmap /nfs/users/nfs_s/sd21/lustre118_link/databases/diamond/prot.accession2taxid \
--max-target-seqs 1 > unmapped_cds.diamond.out


cut -f4 unmapped_proteins.diamond.out | sort | uniq -c | sort -k1nr

```

- 796 of 1000 hits classified as contaminants
- 204 not in uniprot refseq database

- 588 83333 Escherichia coli K-12
- 204  No hit
- 161 224308     Bacillus subtilis
- 13 6239   Caenorhabditis elegans
-  3 10090  Mus musculus
-  3 224911 Bradyrhizobium diazoefficiens
-  3 243274 Thermotoga maritima
-  2 100226 Streptomyces coelicolor
-  2 1111708     Synechocystis sp
-  2 122586 Neisseria meningitidis
-  2 190304 Fusobacterium nucleatum
-  2 208964 Pseudomonas aeruginosa
-  2 224324 Aquifex aeolicus
-  2 7955   Danio rerio
-  1 243090 Rhodopirellula baltica
-  1 243230 Deinococcus radiodurans
-  1 243232 Methanocaldococcus jannaschii
-  1 273057 Saccharolobus solfataricus
-  1 44689 Dictyostelium discoideum
-  1 515635 Dictyoglomus turgidum
-  1 7227 Drosophila melanogaster
-  1 83332 Mycobacterium tuberculosis
-  1 8364 Xenopus tropicalis
-  1 85962 Helicobacter pylori
-  1 9913 Bos taurus

```bash
 cut -f13 unmapped_cds.diamond.out | sort | uniq -c | sort -k1nr
```

-      598 83333 Escherichia coli K-12
-      206
-      161 224308  Bacillus subtilis
-       13 6239
-        3 10090
-        3 224911
-        3 243274
-        2 100226
-        2 1111708
-        2 122586
-        2 190304
-        2 208964
-        2 224324
-        2 7955
-        1 243090
-        1 243230
-        1 243232
-        1 273057
-        1 44689
-        1 515635
-        1 7227
-        1 83332
-        1 8364
-        1 85962
-        1 9913




## Running interproscan on liftoff annotation

```bash
# make a proteins file from annotation
gffread -y PROTEINS.fa -g trichuris_trichiura.fa liftover_annotation.gff3

# remove stop codons from protein sequences - IPS doesnt like it
sed -e 's/\.//g' PROTEINS.fa > tmp; mv tmp PROTEINS.fa

# load module
module load interproscan/5.39-77.0-W01

# run IPS
farm_interproscan -a PROTEINS.fa -o IPS.output.gff

# lift over GO terms from interproscan to GFF
extract_interproscan_go_terms -i IPS.output.gff -e liftover_annotation.gff3



# extract GO terms from mRNAs in GFF
awk '$3=="mRNA" {print $0}' OFS="\t" liftover_annotation.gff3.go.gff | grep "Ontology_term" | cut -f 9 | cut -f3,4 -d ";" | sed -e 's/;/\t/g' -e 's/Name=//g' -e 's/"//g' -e 's/Ontology_term=//g' -e 's/-.*\t/\t/g' -e 's/,/\t/g' -e 's/\..*\t/\t/g' > annotation_GO_per_gene.txt

# work through columns to split multiple go terms per gene into one term per gene, repeating the gene name if multiple GO terms present
awk '{for(i=2; i<=NF; i++) {print $1,$i}}' OFS="\t" annotation_GO_per_gene.txt > annotation_GO_per_gene_split.txt

```
