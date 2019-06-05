#!/usr/bin/Rscript

# Combines vcfresultsPart files having the same newChr (after Glazerized)
# chrname is the chromosome name whose vcfresultsPart objects are to be joined

# Run in Unix as " Rscript combineVcfResultsParts.R ... "

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

library(VariantAnnotation)

args <- commandArgs(TRUE) # project chrname 
# args <- c("Benlim_99.0_SNP", "chrXXI")

project <- args[1]
chrname <- args[2]

chrNumeric <- g$chrname2numeric(chrname)
# [1] 21

cat("\nchrNumeric is", chrNumeric, "\n")

# Results file
vcfresultsNewName <- paste(project, chrname, "vcfNew", "rdd", sep = ".")
# [1] "Benlim_99.0_SNP.chrXXI.vcfNew.rdd"

# glob to identify all vcfresultsPart files, they will have chrNumeric (e.g., "21") before the ".rdd"
vcfresultsPartNames <- paste(project, "chr*.vcfresultsPart", chrNumeric, "rdd", sep = ".")
# [1] "Benlim_99.0_SNP e.chr*.vcfresultsPart.21.rdd"

# The corresponding file names is whatever order (*watch order in subsequent files*)
z <- list.files( pattern=glob2rx( vcfresultsPartNames ), ignore.case=TRUE )
	# [1] "Benlim_99.0_SNP.chrUn.vcfresultsPart.21.rdd" 
	# [2] "Benlim_99.0_SNP.chrXXI.vcfresultsPart.21.rdd"

# Load the vcfresultsPart files in the order given above, place in a temporary list: "parts"
parts <- vector("list", length(z)) # initiate
oldChrNames <- sapply( strsplit(z, split = "[.]"), function(x){x[grepl("chr",x)]})
# [1] "chrUn"  "chrXXI"

# Reorder them
ord <- order(g$chrname2numeric(oldChrNames))
# [1] 2 1
names(parts) <- oldChrNames[ord]
z <- z[ord]

for(i in 1:length(z)){
	load(z[i]) # object name is vcfPart
	parts[[i]] <- vcfPart
	}

# Combine the vcf parts
vcf <- parts[[1]]
if(length(parts) > 1){	
	for(i in 2:length(parts)){
		vcf <- rbind(vcf, parts[[i]])
		}
	}

save(vcf, file = vcfresultsNewName) # object is vcf
# load(vcfresultsNewName) # object is vcf
