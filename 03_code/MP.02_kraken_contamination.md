# Kraken


## Contents
- kraken - modern samples
- kraken - ancient samples

of trimmed reads post mapping
- The mapping shows that there is variable mapping rates, and that for some samples there is very poor mapping.
- This is particularly the case for the ancient samples, which is to be expected to a degree, given they are both old and collected from the environment.
- Kraken might give some insight into this, given they might be heavily contaminated with bacteria etc.

```bash

# load kraken 2
module load kraken2/2.0.8_beta=pl526h6bb024c_0-c1

# run kraken on the modern PE trimmed reads
while read OLD_NAME NEW_NAME; do
     bsub.py 10 kraken2 "kraken2 --db /lustre/scratch118/infgen/pathogen/pathpipe/kraken/minikraken_20190423/minikraken2_v1_8GB \
     --report ${WORKING_DIR}/02_RAW/${NEW_NAME}.kraken2report \
     --paired ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair1.truncated \
     ${WORKING_DIR}/02_RAW/${NEW_NAME}_PE.pair2.truncated";
done < ${WORKING_DIR}/modern.sample_list

# run kraken on the ancient SE trimmed reads
# while read OLD_NAME NEW_NAME; do
#      bsub.py 10 kraken2_SE "kraken2 --db /lustre/scratch118/infgen/pathogen/pathpipe/kraken/minikraken_20190423/minikraken2_v1_8GB
#      --report ${WORKING_DIR}/02_RAW/${NEW_NAME}.kraken2report
#      ${WORKING_DIR}/02_RAW/${NEW_NAME}_SE.truncated";
# done < ${WORKING_DIR}/ancient.sample_list

while read NEW_NAME; do
     bsub.py 5 kraken2_SE "kraken2 --db /lustre/scratch118/infgen/pathogen/pathpipe/kraken/minikraken_20190423/minikraken2_v1_8GB \
     --report ${WORKING_DIR}/02_RAW/${NEW_NAME}.kraken2report \
     ${WORKING_DIR}/02_RAW/${NEW_NAME}_SE.truncated";
done < ${WORKING_DIR}/ancient.sample_list_v2

# once the kraken runs have completed, run multiqc .
multiqc *kraken2report --title kraken

```

[Kraken MultiQC report](../04_analysis/kraken/kraken_multiqc_report.html)


- the output shows most samples have a small degree of contamination based on hits in the kraken database
- non have a lot of contamination, which is slightly surprising
- this alone doesnt explain the mismapping, althoguh, it simply may mean that the putative contaminant is not present in the kraken databased
- could try
     - blasting some sequences
     - mapping to other reference - are there other nematodes in these environmental ancient samples?
