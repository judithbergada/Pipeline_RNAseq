#!/bin/bash

set -e

###################################################
## Mapping reads to the indexed reference genome ##
###################################################

# Relationship between input parameters and the ones used here
outn=$1; fastqfolder=$2; annotf=$3; threads=$4; revfastqfolder=$5

printf "\nMapping...\n"
# Remove temporal directory if exists
rm -rf temporal
rm -rf temporal*

# Remove directories if existing and create new ones with all permissions
rm -rf intermediate_${outn} && \
  mkdir intermediate_${outn} && \
  chmod +xwr intermediate_${outn}
rm -rf outputs_${outn} && \
  mkdir outputs_${outn} && \
  chmod +xwr outputs_${outn}

n=0
for i in ${fastqfolder}/* ; do
  # Get the needed reverse file matching the forward one
  rvffiles=${revfastqfolder}/*
  j=$(echo "${rvffiles[$n]}")
  # Take name of fastqfile ignoring zip, fastq and common part
  namef=$(echo $i | sed 's/\.gz//g' | sed 's/\.fastq//g' | sed 's/\.fq//g')
  # Remove directory from the fastqfile name and keep only the name of file
  namef=$(echo ${namef} | sed 's/.*\///g')

  # Check if file is compressed
  if file --mime-type $i | grep -q gzip; then
    readflag='--readFilesCommand "gunzip -c"'
  fi

  # Perform alignment
  STAR --runThreadN ${threads} \
    --runMode alignReads \
    --genomeDir genome_${outn} \
    --readFilesIn $i $j \
    --twopassMode Basic \
    --outSAMtype None \
    --quantMode GeneCounts \
    --sjdbGTFfile ${annotf} \
    --sjdbGTFfeatureExon CDS \
    --outTmpDir temporal \
    --outFileNamePrefix "${namef}_" \
    ${readflag}

  # Remove temporal directory if exists
  rm -rf temporal
  rm -rf temporal*

  # Get only needed columns of outputs, which contain un-stranded gene counts
  cat ${namef}_ReadsPerGene.out.tab | \
  cut -f1,2 > "intermediate_${outn}/reads_${namef}.tab"

  # Save useful files considering user parameters
  mv "${namef}_Log.final.out" \
  intermediate_${outn}/Statistics_alignment_${namef}.txt

  # Remove files that can lead to problems in future iterations
  rm ${namef}_Log* ${namef}_SJ.out.tab
  rm -r *STAR*
  n=$n+1
done

##################################################
## Prepare final table with the counts per gene ##
##################################################

# Write first column of the table: Gene names
commandpaste="<(sort -V intermediate_${outn}/reads_${namef}.tab | cut -f1)"

# Write header
header="ID"

# Write the rest of the columns: number of counts per SampleName
for i in intermediate_${outn}/reads*; do
  commandpaste="${commandpaste} <(sort -V $i | cut -f2)"
  name=$(echo $i | sed 's/.*intermediate_.*\/reads_/reads_/g' | sed 's/.tab//g')
  header="${header}\t${name}"
done

# Paste everything
eval paste ${commandpaste} > "intermediate_${outn}/counts.txt"
echo -e ${header} | \
cat - intermediate_${outn}/counts.txt > outputs_${outn}/table.counts.txt

# Remove not needed files
rm *ReadsPerGene*
echo "Mapping finished succesfully."
