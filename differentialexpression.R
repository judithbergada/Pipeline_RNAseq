#!/usr/bin/env Rscript

#########################################
## Analysis of differential expression ##
#########################################

# Check if DESeq2 package is properly installed and load it
if (!requireNamespace("DESeq2", quietly = T)) {
  print("Installing DESeq2")
  if (!requireNamespace("BiocManager", quietly = T)) {
    install.packages("BiocManager", repos = "https://stat.ethz.ch/CRAN/")
  }
  BiocManager::install("DESeq2", ask=F, update=F)
  print("DESeq2 has been installed")
}
library("DESeq2")

# Relationship between input parameters and the ones used here
inputdata = commandArgs(trailingOnly=TRUE)
outn = inputdata[1]
conditions_txt = inputdata[2]

print("Performing differential expression analysis...")

#################
## Import data ##
#################

# Open counts file obtained from the alignment
info = read.table(paste("outputs_", outn, "/table.counts.txt", sep=""),
  stringsAsFactors = F, header = T, check.names = F)[-c(1:4),]

# Use first column names as rownames
rownames(info) = info$ID
info = info[,-1]

# Open file of conditions provided by the user
conditions = read.table(conditions_txt, stringsAsFactors = F, header = F)

# Add "reads" in front of the file name
conditions$V1 = paste("reads_", conditions$V1, sep = "")

# Transpose table
conditions = as.data.frame(t(conditions))
colnames(conditions) = as.character(unlist(conditions[1,]))
conditions = conditions[-1,]

# Arrange counts file with information provided by user
info = info[colnames(info) %in% colnames(conditions)]
info = rbind(info, conditions)

# Change name of samples to be only condition or control
colnames(info) = info[nrow(info),]
info = info[-nrow(info),]
info=as.data.frame(apply(info, c(1,2), as.integer))


##############
## Analysis ##
##############

# Take only data for which there are some counts
info_data = info[rowSums(info)>0, ]

# Set conditions for comparison
colData<-data.frame(condition=factor(colnames(info_data)))

# Perform analysis
dds <- DESeqDataSetFromMatrix(info_data, colData, formula(~condition))
dds <- DESeq(dds, betaPrior = TRUE)
res <- results(dds)

# Create dataframe with DESEQ2 results
final_results = data.frame(res)

# Remove rows with NAN values after DESeq
# Reason: DESeq produces NANs to optimise multiple testing
final_results = final_results[!is.na(final_results$padj),]
final_results$Gene_names = rownames(final_results)

# Save data
compname = gsub(".*/", "", conditions_txt)
write.table(final_results,
  file = paste("outputs_", outn, "/Diff_expression_", compname, sep=""),
  quote = F, sep = "\t", row.names = F, col.names = T)

print("Differential expression finished successfully.")

##########
## DONE ##
##########
