#!/usr/bin/Rscript

# To make an index for a .vcf file
# Here I use VariantAnnotation

# # module load R/3.1.2
# # R

args <- commandArgs(TRUE) # first argument will be bamfile name
# args <- "Benlim.chrM.var.vcf"

library(VariantAnnotation, quietly = TRUE)
vcfname <- args[1]
compressedVcfname <- paste(vcfname, ".bgz", sep = "")
zipped <- bgzip(vcfname, compressedVcfname)
idx <- indexTabix(zipped, "vcf")
