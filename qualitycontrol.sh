#!/bin/bash

set -e

#############################################
## Quality control of the sequencing reads ##
#############################################

# Relationship between input parameters and the ones used here
outn=$1; fastqfolder=$2; threads=$3

printf "\nPerforming Quality Control of sequencing reads...\n"

# Remove directory if exists and create a new one with all permissions
rm -rf qc_${outn} && mkdir qc_${outn} && chmod +xwr qc_${outn}

fastqc -q -t ${threads} -o qc_${outn} ${fastqfolder}/*
echo "QC finished successfully."
