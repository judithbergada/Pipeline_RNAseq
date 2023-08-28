#!/bin/bash

set -e
###################
## Assign inputs ##
###################

# Define usage of the script
function print_usage {
  printf """
Usage: RNAseq     [-h or --help]
                  [-f or --fastqfolder]
                  [-a or --annotfolder]
                  [-c or --condfiles]
                  [-r or --revfastqfolder]
                  [-o or --outname]
                  [-t or --threads]
"""
}
# Describe usage of the tool and provide help
function print_help {
  print_usage
  printf """
Optional arguments:
    -h, --help:
                Show this help message and exit.
    -o, --outname:
                Name of your analysis.
                It will be used to name the output files.
                Default: rnaseq.
    -t, --threads:
                Number of threads that will be used.
                It must be an integer.
                Default: 8.
    -r, --revfastqfolder:
                If you are using paired-end reads, this should be the
                path to the folder that contains ALL your REVERSE reads.
                Only FASTQ files should be placed in it.
                Files might (or might not) be compressed.
                Note: ONLY REVERSE reads must be found in this folder, and
                we recommend to name the files exactly as the forward reads but
                adding _2 or _reverse at the end (followed by .fastq or .fq.gz).

Required arguments:
    -f, --fastqfolder:
                Path to the folder that contains ALL your FASTQ files.
                Only FASTQ files should be placed in it.
                Files might (or might not) be compressed.
                Note: if you are using paired-end reads, this should be the
                folder that contains ONLY your FORWARD reads.
    -a, --annotfolder:
                Path to the folder that contains the annotation files.
                Only 2 files should be placed in this folder:
                1) Whole genome FASTA file, 2) GTF file.
    -c, --condfiles:
                Path to a tab-delimited text file describing the
                experimental design or comparisons. 2 columns are needed:
                Column 1: name of each FASTQ file, without extension.
                Column 2: condition of each file (experiment or control).
                You can specify as many files as you want, but
                they must be writen within quotes.
                E.g.: --condfiles 'comparison1.txt comparison2.txt'
"""
}

# Define inputs
for ARGS in "$@"; do
  shift
        case "$ARGS" in
                "--fastqfolder") set -- "$@" "-f" ;;
                "--annotfolder") set -- "$@" "-a" ;;
                "--condfiles") set -- "$@" "-c" ;;
                "--revfastqfolder") set -- "$@" "-r" ;;
                "--outname") set -- "$@" "-o" ;;
                "--threads") set -- "$@" "-t" ;;
                "--help") set -- "$@" "-h" ;;
                *) set - "$@" "$ARGS"
        esac
done

# Define defaults
outn="rnaseq"; threads=8; revfastqfolder=""

# Define all parameters
while getopts 'f:a:c:r::o::t::h' flag; do
        case "${flag}" in
                f) fastqfolder=${OPTARG} ;;
                a) annotfolder=${OPTARG} ;;
                c) condfiles=${OPTARG} ;;
                r) revfastqfolder=${OPTARG} ;;
                o) outn=${OPTARG} ;;
                t) threads=${OPTARG} ;;
                h) print_help
                   exit 1;;
                *) print_usage
                    exit 1;;
        esac
done


##############################
## Identify Software Errors ##
##############################

printf "\nChecking if required software is installed...\n"
# Check installation of the required software.
# If something is missing, show how to install it.
if ! [ -x "$(command -v fastqc)" ]; then
  echo "Missing: fastqc not found"
  echo "Information on the installation:"
  echo "http://www.bioinformatics.babraham.ac.uk/projects/download.html"
  exit 127
fi
if ! [ -x "$(command -v STAR)" ]; then
  echo "Missing: STAR not found"
  echo "Information on the installation:"
  echo "https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
  exit 127
fi
echo "Required software is properly installed."


########################################
## Identify Errors in Inputs Required ##
########################################

printf "\nChecking if required inputs are correct...\n"
# Check if directories providing FASTQ files and ANNOTATION files exist
if [ ! -d ${fastqfolder} ]; then
  echo "Error: --fastqfolder doesn't exist."
  echo "Solution: check if the path to this directory is correct."
  exit 1
fi
if [ ! -d ${annotfolder} ]; then
  echo "Error: --annotfolder doesn't exist."
  echo "Solution: check if the path to this directory is correct."
  exit 1
fi

# Check if only FASTQ files are provided in the fastqfolder
for seqfile in ${fastqfolder}/*; do
  # Get the extension of the file, if file is compressed
  if file --mime-type ${seqfile} | grep -q gzip; then
    filename=${seqfile%.*}
    extension=${filename##*.}
  else # Get the extension of the file, if file is NOT compressed
    extension=${seqfile##*.}
  fi
  # Check if extension is fastq or fq
  extension=$(tr "[:upper:]" "[:lower:]" <<< ${extension})
  if [[ ${extension} != "fq" && ${extension} != "fastq" ]]; then
    echo "Error: --fastqfolder should only contain FASTQ files."
    echo "Solution: remove any other file from this directory."
    exit 1
  fi
done

# Check if FASTA file & GTF file are provided in the annotfolder
i=1
for annotfile in ${annotfolder}/*; do
  extension=${annotfile##*.}
  # Check if extension is fasta, fa or gtf
  extension=$(tr "[:upper:]" "[:lower:]" <<< $extension)
  if [[ ${extension} == "gtf" ]]; then
    annotf=${annotfile}
  elif [[ ${extension} == "fa" || ${extension} == "fasta" ]]; then
    fastaf=${annotfile}
  else
    echo "Error: --annotfolder should only contain a FASTA and a GTF file."
    echo "Solution: remove any other file from this directory."
    exit 1
  fi
  # Check if there are at most 2 files in the folder
  if [[ $i -gt 2 ]]; then
    echo "Error: --annotfolder should only contain 2 files."
    echo "Solution: remove the other files from this directory."
    exit 1
  fi
  i=$((i+1))
done
# Check that both FASTA and GTF are provided
if [[ ${annotf} == "" || ${fastaf} == "" ]]; then
  echo "Error: --annotfolder does not contain a FASTA and a GTF file."
  echo "Solution: add FASTA and GTF file to the directory."
  exit 1
fi

# Check if input file showing the experimental design exists
for comparison in $condfiles; do
  if [ ! -f ${comparison} ]; then
    echo "Error: --condfiles ${comparison} does not exist."
    echo "Solution: check if the path to this file is correct."
    exit 1
  fi
done

echo "Required inputs seem correct."


########################################
## Identify Errors in Optional Inputs ##
########################################

printf "\nChecking if optional inputs are correct...\n"
# Check if the number of threads is an integer
if ! [[ ${threads} =~ ^[0-9]+$ ]]; then
  echo "Error: --threads is not an integer."
  echo "Solution: remove this optional parameter or use an integer."
  exit 1
fi
# Check if directory providing reverse FASTQ files exist
if [[ ${revfastqfolder} != "" ]]; then
  if [ ! -d ${revfastqfolder} ]; then
    echo "Error: --revfastqfolder doesn't exist."
    echo "Solution: check if the path to this directory is correct."
    exit 1
  fi
  # Check if only FASTQ files are provided in the revfastqfolder
  for revseqfile in ${revfastqfolder}/*; do
    # Get the extension of the file, if file is compressed
    if file --mime-type ${revseqfile} | grep -q gzip; then
      revfilename=${revseqfile%.*}
      revextension=${revfilename##*.}
    else # Get the extension of the file, if file is NOT compressed
      revextension=${revseqfile##*.}
    fi
    # Check if extension is fastq or fq
    revextension=$(tr "[:upper:]" "[:lower:]" <<< ${revextension})
    if [[ ${revextension} != "fq" && ${revextension} != "fastq" ]]; then
      echo "Error: --revfastqfolder should only contain FASTQ files."
      echo "Solution: remove any other file from this directory."
      exit 1
    fi
  done
fi
echo "Optional inputs seem correct."


###########################
## Start RNAseq analysis ##
###########################

printf "\nRNAseq analysis in progress...\n"
# Get the directory of this script to be able to call the others
currentdir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Perform quality control
bash ${currentdir}/qualitycontrol.sh \
  "${outn}" "${fastqfolder}" "${threads}" "${revfastqfolder}" # Inputs
if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if error

# Index the genome files
bash ${currentdir}/genomegenerate.sh \
  "${outn}" "${fastqfolder}" "${fastaf}" "${annotf}" "${threads}" # Inputs
if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.

if [[ ${revfastqfolder} != "" ]]; then
  bash ${currentdir}/alignment_pairedreads.sh \
    "${outn}" "${fastqfolder}" "${annotf}" "${threads}" "${revfastqfolder}"
  if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.
else
  bash ${currentdir}/alignment.sh \
    "${outn}" "${fastqfolder}" "${annotf}" "${threads}" # Inputs
  if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.
fi

# Perform analysis of differential expression
for comparison in $condfiles; do
  Rscript --vanilla $currentdir/differentialexpression.R \
          ${outn} ${comparison}
  if [[ $(echo $?) != 0 ]]; then exit 1; fi # Exit if there has been an error.
done

echo "Analysis completed sucessfully!"

##########
## DONE ##
##########
