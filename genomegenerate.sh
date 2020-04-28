#!/bin/bash

set -e

###############################
## Indexing reference genome ##
###############################

# Relationship between input parameters and the ones used here
outn=$1; fastqfolder=$2; fastaf=$3; annotf=$4; threads=$5

printf "\nIndexing the reference genome...\n"
# Remove directory if exists and create a new one with all permissions
rm -rf genome_${outn} && \
  mkdir genome_${outn} && \
  chmod +xwr genome_${outn}

# Get one of the FASTQ files in order to calculate the read lenght later
seqfilename=$(ls ${fastqfolder} | head -n1)
seqfile="${fastqfolder}/${seqfilename}"

# Check if file is compressed and calculate length of the sequencing reads
if file --mime-type ${seqfile} | grep -q gzip; then
  num_bases=$(zcat < ${seqfile} | head | awk 'NR==2' | wc -c)
else
  num_bases=$(head ${seqfile} | awk 'NR==2' | wc -c)
fi
# Compute needed parameter according to STAR manual
num_overhang=$(echo ${num_bases} - 1 | bc -l)

# Index genome
STAR --runThreadN ${threads} \
--runMode genomeGenerate \
--genomeSAindexNbases 8 \
--outFileNamePrefix genome_${outn}/ \
--outTmpDir temporal \
--genomeDir genome_${outn} \
--genomeFastaFiles ${fastaf} \
--sjdbGTFfile ${annotf} \
--sjdbGTFfeatureExon CDS \
--sjdbOverhang $num_overhang

# Remove temporal directory if exists
rm -rf temporal*
echo "Reference genome indexed sucessfully."
