#!/usr/bin/Rscript

# To make an index for a .vcf file
# Here I use VariantAnnotation

# module load R/3.1.2
# R

args <- commandArgs(TRUE) # first argument will be vcf.gz name (might not have the gz)
# args <- "SculpinNoSculpin.var.vcf.gz"

library(VariantAnnotation, quietly = TRUE)

# File names
vcfname <- args[1]
root <- sub("[.]gz$", "", vcfname) # remove the "gz" if present
compressedVcfname <- paste(root, "bgz", sep = ".") # save to this instead of .gz, which already exists

# compress and save file name
bgzname <- bgzip(vcfname, compressedVcfname)

# Create the index
idx <- indexTabix(bgzname, "vcf")

# optional: rename .bgz to .gz (this will overwrite original gz file)
# file.rename(bgzname, sub("vcf.bgz", "vcf.gz", bgzname))
# file.rename(idx, sub("vcf.bgz.tbi", "vcf.gz.tbi", idx))
