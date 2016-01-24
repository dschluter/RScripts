#!/usr/bin/Rscript

# Plots genome scanes for interesting pairs by chromosomes using "g$plotSlidewinInterestingPairsByChr"
# Uses the files named " project.chr.pop1.pop2.blockstats.500.rdd "
# Requires file "glazerFileS4 NewScaffoldOrder.csv" in current directory

# qsub -I -l walltime=04:00:00 -l mem=4gb 
# module load R/3.1.2
# R

# setwd("~/Desktop")
# git("genome.r")

# set ymax to zero to ignore this argument

# these are the defaults anyway in g$plotSlidewinInterestingPairsByChr
orderChr <- TRUE
Glazerize <- TRUE
scafFile <- "glazerFileS4 NewScaffoldOrder.csv"

args <- commandArgs(TRUE) # project chrname method stepsize nsteps.per.window windowNmin ymax interestingPairs[list]
# args <- c( "BenlimAllMarine", "chrXXI", "vara", 500, 5, 100, 0.03, list(c("paxl", "paxb"), c("pril", "prib")) )

interestingPairs <- args[8:length(args)]
args <- unlist(args[1:7])

project <- args[1]
chrname <- args[2]
method <- args[3]
stepsize <- as.integer(args[4])
nsteps.per.window <- as.integer(args[5])
windowNmin <- as.integer(args[6])
ymax <- as.numeric(args[7]) # set ymax to zero to ignore this argument

# Check that all blockstats files are present
z <- sapply(interestingPairs, paste, collapse = ".")
z <- apply(expand.grid(chrname, z), 1, paste, collapse = ".")
z <- paste(project, z, "blockstats", stepsize, "rdd", sep = ".")
if( !all(file.exists(z)) ){
	print( z[!file.exists(z)] )
	stop("The above blockstats files are missing\n")
	}

g$plotSlidewinInterestingPairsByChr(project = project, chrname = chrname, method = method, ymax = ymax,
		interestingPairs = interestingPairs, stepsize = 500, nsteps.per.window = 5, windowNmin = 100, 
		orderChr = orderChr, Glazerize = Glazerize, scafFile = scafFile)
