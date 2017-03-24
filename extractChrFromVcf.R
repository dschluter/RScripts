#!/usr/bin/Rscript

# Filters bgzipped Vcf file, pulls out a single CHR 

# Run in Unix as " Rscript extractChrFromVcf.R ... "

# setwd("~/Desktop")
# git("genome.r")

# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

# Expect to read these arguments from args
chrname 	<- NULL
vcfFile 	<- NULL
genome 		<- "gasAcu1pitx1new.fa"
chunksize 	<- 10000

args <- commandArgs(TRUE)
# args <- c("chrname=chrM", "vcfFile= Benlim_99.0_SNP.vcf.gz", "genome=gasAcu1pitx1new.fa")

# Parses the args into a data frame with two columns (V1=left and V2=right of each "=" sign)
# and then assigns V2 to variables whose names are in V1 

x <- read.table(text = args, sep = "=", colClasses = "character", strip.white = TRUE)
for(i in 1:nrow(x)){ assign(x[i,1], x[i,2]) }

if(is.null(vcfFile)) stop("Provide vcfFile= in arguments")
if(is.null(chrname)) stop("Provide chrname= in arguments")
cat("\nchrname is", chrname, "\n")

tbiFile <- paste(vcfFile, "tbi", sep = ".")
if(!file.exists(tbiFile)) stop(paste("tabix file", tbiFile, "missing"))

g$extractChrFromVcfFile(chrname = chrname, bgzfile = vcfFile, genome = genome, chunksize = chunksize)