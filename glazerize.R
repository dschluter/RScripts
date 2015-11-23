#!/usr/bin/Rscript

# MAKE THIS A PART OF MAKING VCFRESULTS and GOODINVARIANTS 

# Converts cordinates to Glazer et al assembly.
# Doesn't split yet

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrUn")
# args <- c("BenlimPax22pacMar7", "chrXXI")

project <- args[1]
chrname <- args[2]
chrno 				<- gsub("^chr", "", chrname)
chrNumeric <- chrno
chrNumeric[chrno != "Un"] <- as.numeric( as.roman( chrNumeric[chrno != "Un"] ) )

vcfresultsfile	<- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")

cat("\nLoading vcfresults file\n")
load(file = vcfresultsfile)   # object is "vcfresults"
# lapply(vcfresults, object.size)

# names(vcfresults)
# [1] "groupnames"        "groupcodes"        "control"           "vcf"              
# [5] "altUsedList"       "snpTypeList"       "alleleFreqByGroup"

library(VariantAnnotation)

pos <- start(rowData(vcfresults$vcf))

z <- g$glazerConvertCoordinate(rep(chrNumeric, length(pos)), pos)
vcfresults$newChr <- z[ ,"newChr"]
vcfresults$newPos <- z[ ,"newPos"]

names(vcfresults)

*now SPLIT VCFRESULTS*