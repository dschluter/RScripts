#!/usr/bin/Rscript

# Takes file "Benlim.chrM.DP.inv.gz" and drops fish indicated by "drop" argument
#	and writes to "Benlim.chrM.sel.inv.gz"

# Run in Unix as " Rscript selectInvariants.R invFile drop "

# qsub -I -l walltime=03:00:00 -l mem=2gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

# Expect to read these arguments from args

invariantsummaryname <- NULL
dropfish <- NULL

args <- commandArgs(TRUE)
# args <- c("Benlim.chrM.DP.inv.gz", "Marine-Pac-Bamfield-VI17-Sara,Marine-Pac-Oyster-12-Sara,Marine-Pac-Salmon-01-Sara")

invariantsummaryname<- args[1]
dropfish			<- args[2:length(args)]
dropfish			<- unlist(strsplit(dropfish, split=",")) # breaks apart drop if not already a vector

if(is.null(invariantsummaryname)) stop("Provide invariantsummaryname in arguments")
if(is.null(dropfish)) stop("Provide fishnames to drop in arguments")

selfilename    <- sub("DP", "sel", invariantsummaryname)

nLinesAtaTime <- 100000
# nLinesAtaTime <- 10000

INFILE <- file(invariantsummaryname, open = "r")
OUTFILE <- gzfile(selfilename, "w")

x <- readLines(con = INFILE, n = nLinesAtaTime)
nlines <- length(x)

# Break up the lines, store in a list
x1 <- strsplit(x, split = "\t")
headerline <- x1[[1]]

# Figure out which elements belong to those being dropped
dropwhich <- sapply(dropfish, function(x){grep(x, headerline)})

if(length(dropwhich) != length(dropfish)) cat("\nWarning: dropwhich not same length as dropfish vector")

# Loop
while(nlines >0){
	
	# Drop the offending fish
	x2 <- lapply(x1, function(x1){x1[-dropwhich]})
	x3 <- sapply(x2, function(x2){paste(x2, collapse="\t")})
	
	writeLines(x3, OUTFILE)
	
	x <- readLines(con = INFILE, n = nLinesAtaTime)
	nlines <- length(x)
	x1 <- strsplit(x, split = "\t")

	} # end while loop

close(INFILE)
close(OUTFILE)
