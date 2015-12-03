#!/usr/bin/Rscript

# Combines goodInvPart files having the same newChr (after Glazerized)
# chrname is the chromosome name whose goodInvPart objects are to be joined

# Run in Unix as " Rscript combinegoodInvParts.R ... "

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

args <- commandArgs(TRUE) # project chrname 
# args <- c("BenlimAllMarine", "chrXXI")

project <- args[1]
chrname <- args[2]

chrNumeric <- g$chrname2numeric(chrname)
# [1] 21

# Results file
goodInvNewName <- paste(project, chrname, "goodInvNew", "rdd", sep = ".")
# [1] "BenlimAllMarine.chrXXI.goodInvNew.rdd"

# glob to identify all goodInvPart files, they will have chrNumeric (e.g., "21") before the ".rdd"
goodInvPartNames <- paste(project, "chr*.goodInvPart", chrNumeric, "rdd", sep = ".")
# [1] "BenlimAllMarine.chr*.goodInvPart.21.rdd"

# The corresponding file names
z <- list.files( pattern=glob2rx( goodInvPartNames ), ignore.case=TRUE )
	# [1] "BenlimAllMarine.chrUn.goodInvPart.21.rdd" 
	# [2] "BenlimAllMarine.chrXXI.goodInvPart.21.rdd"

# Load the goodInvPart files, place in a temporary list: "parts"
parts <- vector("list", length(z)) # initiate
oldChrNames <- sapply( strsplit(z, split = "[.]"), function(x){x[2]})
# [1] "chrUn"  "chrXXI"
names(parts) <- oldChrNames

for(i in 1:length(z)){
	load(z[i]) # name of object is goodInvariantsPart
	parts[[i]] <- goodInvariantsPart
	}
	
goodInvariants <- do.call("rbind", parts)

save(goodInvariants, file = goodInvNewName) # object is goodInvariants
# load(goodInvNewName) # object is goodInvariants
