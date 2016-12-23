#!/usr/bin/Rscript

# Need to make an index if selecting variants from a split compressed .var.vcf.gz file
# Here I use VariantAnnotation

# # module load R/3.1.2
# # R

args <- commandArgs(TRUE) # first argument will be bamfile name
# args <- "Benlim.chrM.var.vcf.gz"

library(VariantAnnotation, quietly = TRUE)
vcfname <- args[1]
compressedVcfname <- sub(".gz$", ".bgz", vcfname)
zipped <- bgzip(vcfname, compressedVcfname)
idx <- indexTabix(zipped, "vcf")
