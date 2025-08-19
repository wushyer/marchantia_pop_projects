# Mapping


## Contents
- Project setup
- Prepare reference genome
- Raw sequence data
- metadata
- trimming of the raw reads
- FIX: merge duplicate ancient read sets
- Mapping


## Project setup

```bash
# working directory
cd /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura
WORKING_DIR=/nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura

mkdir 00_SCRIPTS 01_REF 02_RAW 03_MAPPING 04_VARIANTS 05_ANALYSIS
```

## Prepare reference genome
```bash
cd /nfs/users/nfs_s/sd21/lustre118_link/trichuris_trichiura/01_REF

# make a generic bwa index for mapping later on
bwa index trichuris_trichiura.fa

# make a generic dict file for SNP calling later on
samtools dict trichuris_trichiura.fa > trichuris_trichiura.dict

```


## Raw sequencing data
- Sequencing was performed at the Centre for Geogenetics (CGG), University of Copenhagen
- I received the raw sequencing data from Peter on a hard drive, and copied all to the lustre environment as is into the directory 02_RAW
- want to collect all of the FASTQ files into one place, and then begin processing them

```bash
cd ${WORKING_DIR}/02_RAW

for i in *gz; do
     echo ${i%_???.fastq.gz};
done | sort | uniq | while read NAME; do
          cat ${NAME}* > ${NAME}.merged.fastq.gz;
     done &

# sample used in the Foth 2014 genome paper
pf data --type lane --id 6929_8 -l ./ --filetype fastq

```

- this generated a total of 130 R1 files and 60 R2 files
- the difference in number reflects the fact that the "modern" samples have been sequenced using PE reads whereas the "ancient" samples have been sequenced using SE reads. Will need to treat these a little differently later on.
- there is also a difference in the number of samples with useable data (based on some sample informaiton spreadsheets from MS / PJ), and the total number
     - some samples are controls
     - some samples presumably didnt work, and so we excluded from the preliminary analyses
     - some samples are duplicates, and so will need to be merged at some point too.
---

## metadata
- I have begun to curate the metadata, to get a better idea of what data is actually available, and to ensure that I have consistent informaiton for all samples.
- The existing naming scheme is not very informative to me at least, and so will rename everything with a simple format
- GoogleSheet:  https://docs.google.com/spreadsheets/d/1PiapiaZZw0g0i3lN0feXxqEVupvaAuOm0IViXnZPFqk/edit?usp=sharing

---
## trimming of the raw reads
- trimming using AdapterRemoval, which seems to be used for a few different ancient DNA projects. I think it is because it trims and merges, which generally improves the mapping scores off some poor qual end of reads.
- Tool: https://buildmedia.readthedocs.org/media/pdf/adapterremoval/latest/adapterremoval.pdf

```bash
cd ${WORKING_DIR}/02_RAW

### run the trimming

# main samples
while read OLD_NAME NEW_NAME; do
     bsub.py --threads 4 20 adapter_remove_modern "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_PE.sh ${OLD_NAME} ${NEW_NAME}";
done < ${WORKING_DIR}/modern.sample_list

while read OLD_NAME NEW_NAME; do
     bsub.py --threads 4 20 adapter_remove_ancient "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_SE.sh ${OLD_NAME} ${NEW_NAME}";
done < ${WORKING_DIR}/ancient.sample_list

# other samples
while read OLD_NAME NEW_NAME; do
     bsub.py --threads 4 20 adapter_remove_others_PE "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_PE.sh ${OLD_NAME} ${NEW_NAME}";
done < ${WORKING_DIR}/others_PE.sample_list

while read OLD_NAME NEW_NAME; do
     bsub.py --threads 4 20 adapter_remove_others_SE "${WORKING_DIR}/00_SCRIPTS/run_adapter_remove_SE.sh ${OLD_NAME} ${NEW_NAME}";
done < ${WORKING_DIR}/others_SE.sample_list

```

where "run_adapter_remove_PE.sh" is:

```bash
#!/bin/bash
# adaptor remove PE - modern samples
OLD_NAME=${1}
NEW_NAME=${2}

/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/adapterremoval/bin/AdapterRemoval \
--file1 ${OLD_NAME}_R1.merged.fastq.gz \
--file2 ${OLD_NAME}_R2.merged.fastq.gz \
--basename ${NEW_NAME}_PE \
--trimns --trimqualities --collapse --threads 4

```

and where "run_adapter_remove_SE.sh" is:

```bash
#!/bin/bash
# single end - ancient samples
OLD_NAME=${1}
NEW_NAME=${2}

/nfs/users/nfs_s/sd21/lustre118_link/software/ANCIENT/adapterremoval/bin/AdapterRemoval \
--file1 ${OLD_NAME}_R1.merged.fastq.gz \
--basename ${NEW_NAME}_SE \
--trimns --trimqualities --threads 4

```



## FIX: merge duplicate ancient read sets
These seem to have been a single sample, extracted twice (or in two ways, perhaps two washes of a column) and sequenced individually. Given the low coverage, these should be merged (I think).

```bash
# merge the duplicated set, renaming with a unique name that refers to both the originals, and then remove the originals
cat AN_DNK_COG_EN_001_SE.truncated AN_DNK_COG_EN_002_SE.truncated > AN_DNK_COG_EN_0012_SE.truncated; rm AN_DNK_COG_EN_001_SE.truncated AN_DNK_COG_EN_002_SE.truncated
cat AN_DNK_COG_EN_003_SE.truncated AN_DNK_COG_EN_004_SE.truncated > AN_DNK_COG_EN_0034_SE.truncated; rm AN_DNK_COG_EN_003_SE.truncated AN_DNK_COG_EN_004_SE.truncated
cat AN_DNK_COG_EN_005_SE.truncated AN_DNK_COG_EN_006_SE.truncated > AN_DNK_COG_EN_0056_SE.truncated; rm AN_DNK_COG_EN_005_SE.truncated AN_DNK_COG_EN_006_SE.truncated
cat AN_DNK_COG_EN_007_SE.truncated AN_DNK_COG_EN_008_SE.truncated > AN_DNK_COG_EN_0078_SE.truncated; rm AN_DNK_COG_EN_007_SE.truncated AN_DNK_COG_EN_008_SE.truncated
cat AN_DNK_COK_EN_001_SE.truncated AN_DNK_COK_EN_002_SE.truncated > AN_DNK_COK_EN_0012_SE.truncated; rm AN_DNK_COK_EN_001_SE.truncated AN_DNK_COK_EN_002_SE.truncated
cat AN_DNK_COK_EN_003_SE.truncated AN_DNK_COK_EN_004_SE.truncated > AN_DNK_COK_EN_0034_SE.truncated; rm AN_DNK_COK_EN_003_SE.truncated AN_DNK_COK_EN_004_SE.truncated
cat AN_DNK_OBM_EN_001_SE.truncated AN_DNK_OBM_EN_002_SE.truncated > AN_DNK_OBM_EN_0012_SE.truncated; rm AN_DNK_OBM_EN_001_SE.truncated AN_DNK_OBM_EN_002_SE.truncated
cat AN_DNK_OBM_EN_003_SE.truncated AN_DNK_OBM_EN_004_SE.truncated > AN_DNK_OBM_EN_0034_SE.truncated; rm AN_DNK_OBM_EN_003_SE.truncated AN_DNK_OBM_EN_004_SE.truncated
cat AN_DNK_OBM_EN_005_SE.truncated AN_DNK_OBM_EN_006_SE.truncated > AN_DNK_OBM_EN_0056_SE.truncated; rm AN_DNK_OBM_EN_005_SE.truncated AN_DNK_OBM_EN_006_SE.truncated
cat AN_DNK_OBM_EN_007_SE.truncated AN_DNK_OBM_EN_008_SE.truncated > AN_DNK_OBM_EN_0078_SE.truncated; rm AN_DNK_OBM_EN_007_SE.truncated AN_DNK_OBM_EN_008_SE.truncated
cat AN_DNK_OBM_EN_009_SE.truncated AN_DNK_OBM_EN_010_SE.truncated > AN_DNK_OBM_EN_0910_SE.truncated; rm AN_DNK_OBM_EN_009_SE.truncated AN_DNK_OBM_EN_010_SE.truncated
cat AN_LTU_VIL_EN_001_SE.truncated AN_LTU_VIL_EN_002_SE.truncated > AN_LTU_VIL_EN_0012_SE.truncated; rm AN_LTU_VIL_EN_001_SE.truncated AN_LTU_VIL_EN_002_SE.truncated
cat AN_NLD_KAM_EN_001_SE.truncated AN_NLD_KAM_EN_002_SE.truncated > AN_NLD_KAM_EN_0012_SE.truncated; rm AN_NLD_KAM_EN_001_SE.truncated AN_NLD_KAM_EN_002_SE.truncated
cat AN_NLD_KAM_EN_003_SE.truncated AN_NLD_KAM_EN_004_SE.truncated > AN_NLD_KAM_EN_0034_SE.truncated; rm AN_NLD_KAM_EN_003_SE.truncated AN_NLD_KAM_EN_004_SE.truncated
cat AN_NLD_ZWO_EN_001_SE.truncated AN_NLD_ZWO_NA_002_SE.truncated > AN_NLD_ZWO_EN_0012_SE.truncated; rm AN_NLD_ZWO_EN_001_SE.truncated AN_NLD_ZWO_NA_002_SE.truncated

# make a new sample list to work from using new names
ls -1 AN* | cut -c-18 > ../ancient.sample_list_v2

```

### FIX: additional samples to include

```bash
# found some additional samples I had not initially included - they were misclassified in the "other", so correcting them and add to the "ancient.sample_list_v2"
mv OTHER_MSOE43_027_SE.truncated AN_DNK_VIB_EN_001_SE.truncated
mv OTHER_MSOE44_028_SE.truncated AN_DNK_VIB_EN_002_SE.truncated
mv OTHER_MSOE45_029_SE.truncated AN_DNK_VIB_EN_003_SE.truncated
mv OTHER_MSOE46B_030_SE.truncated AN_DNK_VIB_EN_004_SE.truncated
mv OTHER_MSOE46C_031_SE.truncated AN_DNK_VIB_EN_005_SE.truncated

cat AN_DNK_VIB_EN_001_SE.truncated AN_DNK_VIB_EN_002_SE.truncated > AN_DNK_VIB_EN_0012_SE.truncated; rm AN_DNK_VIB_EN_001_SE.truncated AN_DNK_VIB_EN_002_SE.truncated
cat AN_DNK_VIB_EN_003_SE.truncated AN_DNK_VIB_EN_004_SE.truncated AN_DNK_VIB_EN_005_SE.truncated > AN_DNK_VIB_EN_345_SE.truncated; rm AN_DNK_VIB_EN_003_SE.truncated AN_DNK_VIB_EN_004_SE.truncated AN_DNK_VIB_EN_005_SE.truncated

ls -1 AN_DNK_VIB_EN_0012_SE.truncated AN_DNK_VIB_EN_345_SE.truncated | cut -c-18 > ../ancient.sample_list_v2


# clean up
rm *discarded *settings

```


## Mapping
- Need to map ancient and modern samples a little differently, based on the fact that ancient samples only have SE reads, whereas modern samples are PE.
- Below are two mapping scripts for each approach.

### script for mapping the ancient samples - these are all single end (SE) reads

```bash
#!/bin/bash
# map SE reads
OLD_NAME=${1}
NEW_NAME=${2}

# map the SE reads
bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M ${WORKING_DIR}/01_REF/trichuris_trichiura.fa ${WORKING_DIR}/02_RAW/${NEW_NAME}_SE.truncated |\
     samtools view --threads 4 -b - |\
     samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.SE.bam - ;

# mark duplicate reads
java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=${NEW_NAME}.tmp.sort.SE.bam OUTPUT=${NEW_NAME}.tmp2.bam METRICS_FILE=${NEW_NAME}.tmp.metrics TMP_DIR=$PWD/${NEW_NAME}.tmp;

# generate mapping stats from unfilted bam
samtools flagstat ${NEW_NAME}.tmp2.bam > ${NEW_NAME}.flagstat;

# filter bam to only contain mapped reads
samtools view --threads 4 -F 4 -b -o ${NEW_NAME}.bam ${NEW_NAME}.tmp2.bam;

# index the filtered bam
samtools index -b ${NEW_NAME}.bam;

# clean up unnecessary tmp files
rm -r ${NEW_NAME}.*tmp*

```


### script for mapping modern samples - these are all paired-end (PE) reads

```bash
#!/bin/bash
# map PE reads
OLD_NAME=${1}
NEW_NAME=${2}

# map the PE reads
bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M ${WORKING_DIR}/01_REF/trichuris_trichiura.fa ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair1.truncated ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair2.truncated |\
     samtools view --threads 4 -b - |\
     samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.PE.bam - ;

#merge SE reads that are either merged PE reads, or singletons after adapter filtering
cat ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.singleton.truncated ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.collapsed.truncated ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.collapsed > ${NEW_NAME}_SE.tmp.fastq ;

# map SE reads
bwa mem -t 4 -R $(echo "@RG\tRG:${NEW_NAME}\tID:${NEW_NAME}\tSM:${NEW_NAME}") -Y -M ${WORKING_DIR}/01_REF/trichuris_trichiura.fa ${NEW_NAME}_SE.tmp.fastq | \
     samtools view --threads 4 -b - |\
     samtools sort --threads 4 -o ${NEW_NAME}.tmp.sort.SE.bam - ;

# merge PE and SE bams
samtools merge ${NEW_NAME}.tmp.bam ${NEW_NAME}.tmp.sort.PE.bam ${NEW_NAME}.tmp.sort.SE.bam;

# mark duplicate reads
java -Xmx20g -jar /nfs/users/nfs_s/sd21/lustre118_link/software/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=${NEW_NAME}.tmp.bam OUTPUT=${NEW_NAME}.tmp2.bam METRICS_FILE=${NEW_NAME}.tmp.metrics TMP_DIR=$PWD/${NEW_NAME}.tmp;

# generate mapping stats from unfilted bam
samtools flagstat ${NEW_NAME}.tmp2.bam > ${NEW_NAME}.flagstat;

# filter bam to only contain mapped reads
samtools view --threads 4 -F 4 -b -o ${NEW_NAME}.bam ${NEW_NAME}.tmp2.bam;

# index the filtered bam
samtools index -b ${NEW_NAME}.bam;

# clean up unnecessary tmp files
rm -r ${NEW_NAME}.*tmp*

```


```bash

cd ${WORKING_DIR}/03_MAPPING

# run the mapping jobs for modern and ancient samples
while read OLD_NAME NEW_NAME; do
     bsub.py --threads 4 20 mapping_modern "${WORKING_DIR}/00_SCRIPTS/run_map_modern_PE.sh ${OLD_NAME} ${NEW_NAME}" ;
done < ${WORKING_DIR}/modern.sample_list

while read OLD_NAME NEW_NAME; do
     bsub.py --threads 4 20 mapping_otherPE "${WORKING_DIR}/00_SCRIPTS/run_map_modern_PE.sh ${OLD_NAME} ${NEW_NAME}" ;
done < ${WORKING_DIR}/others_PE.sample_list

# run the mapping jobs for the control and other samples
#while read OLD_NAME NEW_NAME; do bsub.py --threads 4 20 mapping_ancient "${WORKING_DIR}/00_SCRIPTS/run_map_ancient_SE.sh ${OLD_NAME} ${NEW_NAME}" ; done < ${WORKING_DIR}/ancient.sample_list_v2
while read OLD_NAME NEW_NAME; do
     bsub.py --threads 4 20 mapping_otherSE "${WORKING_DIR}/00_SCRIPTS/run_map_ancient_SE.sh ${OLD_NAME} ${NEW_NAME}" ;
done < ${WORKING_DIR}/others_SE.sample_list

mkdir MAPPED_OTHER
mv OTHER_M* MAPPED_OTHER

mkdir MAPPED_CONTROL
mv CONTROL_* MAPPED_CONTROL


# rerun of ancient samples due to name change
while read NEW_NAME; do
     bsub.py --threads 4 20 mapping_ancient "${WORKING_DIR}/00_SCRIPTS/run_map_ancient_SE.sh NULL ${NEW_NAME}" ;
done < ${WORKING_DIR}/ancient.sample_list_v2



multiqc *flagstat --title mapping

```

[Mapping MultiQC report](../04_analysis/mapping_multiqc_report.html)
