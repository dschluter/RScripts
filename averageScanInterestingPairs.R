#!/usr/bin/Rscript

# Averages sliding window results of multiple pairwise genome scans and plots the results to a pdf file.
# Separate sliding windows are carried out for each pair, and then an unweighted average is calculated.

# args are project, pairtype, method, chrname, ymax)

# see "scanInterestingPairsByChr.R" for information on arguments

# Analyses files named " project.chr.pop1.pop2.blockstats.stepsize.rdd "

# qsub -I -l walltime=04:00:00 -l mem=4gb 
# module load R/3.1.2
# R

# setwd("~/Desktop")
# git("genome.r")

# set ymax to zero to ignore this argument

args <- commandArgs(TRUE) 
# args <- c( "BenlimAllMarine", "species-pairs", "vara", "all", 0.03)
# args <- c( "BenlimAllMarine", "species-pairs", "vara", "all", 0)
# args <- c( "BenlimAllMarine", "species-pairs", "css", "chrVIIpitx1", 0.04)
# args <- c( "BenlimAllMarine", "species-pairs", "css", "all", 0.04)
# args <- c( "BenlimAllMarine", "species-pairs", "fst", "chrXXI", 1)

project <- args[1]
pairtype <- args[2]
method <- args[3]
chromosomes <- args[4]
ymax <- as.numeric(args[5])

# these are defaults, not arguments
stepsize = 500
nsteps.per.window = 5  	# window size is (nsteps.per.window)*(stepsize), e.g., 5*500 = 2500
windowNmin = 100 		# windowNmin is minimum number of good bases in window
orderChr <- TRUE
Glazerize <- TRUE
scafFile <- "glazerFileS4 NewScaffoldOrder.csv"

if(pairtype == "species-pairs") 
	interestingPairs <- list(
		c("paxl", "paxb"), c("pril", "prib"), c("qryl", "qryb"), c("ensl", "ensb")
		)

if(pairtype == "marinepac-pairs") 
	interestingPairs <- list(
		c("marine-pac", "paxb"), c("marine-pac", "prib"), c("marine-pac", "qryb"), c("marine-pac", "ensb"),
		c("marine-pac", "paxl"), c("marine-pac", "pril"), c("marine-pac", "qryl"), c("marine-pac", "ensl")
		)

if(pairtype == "marinepac-benthic") 
	interestingPairs <- list(
		c("marine-pac", "paxb"), c("marine-pac", "prib"), c("marine-pac", "qryb"), c("marine-pac", "ensb")
		)

if(pairtype == "marinepac-limnetic") 
	interestingPairs <- list(
		c("marine-pac", "paxl"), c("marine-pac", "pril"), c("marine-pac", "qryl"), c("marine-pac", "ensl")
		)

if(pairtype == "solitary-benthic") 
	interestingPairs <- list(
		c("paxb", "solitary"), c("prib", "solitary"), c("qryb", "solitary"), c("ensb", "solitary")
		)

npairs <- length(interestingPairs)

# chrname must be a single chromosome (e.g., "chrXXI") or "all"
# chrname <- c("chrXXI", "chrXX")
if(chromosomes == "all") 
	chrname <- gsub("[.]fa", "", list.files( pattern=glob2rx( "chr*.fa"), ignore.case=TRUE )) else 
	chrname <- chromosomes[1]

# Check that all blockstats files are present
z <- sapply(interestingPairs, paste, collapse = ".")
z <- apply(expand.grid(chrname, z), 1, paste, collapse = ".")
z <- paste(project, z, "blockstats", stepsize, "rdd", sep = ".")
if( !all(file.exists(z)) ){
	print( z[!file.exists(z)] )
	stop("The above blockstats files are missing\n")
	}

if(Glazerize) x <- read.csv(scafFile)

# Order the chromosomes by their numeric values, leaving out the specially named ones
if(orderChr){
	chrSpecial <- intersect(chrname, c("chrVIIpitx1", "chrUn", "chrM"))
	if(length(chrSpecial) > 0) chrname <- chrname[!is.element(chrname, chrSpecial)]
	chrNumeric <- sapply(chrname, g$chrname2numeric)
	chrname <- chrname[ order(chrNumeric) ]
	if(length(chrSpecial) > 0) chrname <- c(chrname, chrSpecial)
	}

# Plot 3 scans per page
plotsPerPage <- 3
pdf( paste(project, pairtype, method, "meanscan", "pdf", sep = ".") )
par(mfrow = c(plotsPerPage, 1), mar = c(2, 2, 1.5, 1) + 0.1)

# If ymax is set, include in header
ymaxHeader <- ""
if(ymax > 0) ymaxHeader <- paste("(ymax ", ymax, ")", sep = "")

library(zoo)
if(tolower(method) == "css") library(gdata)

for(i in chrname){
	# i <- "chrXXI"
	# Grab data needed to include an underline to indicate location of old assembly of this chromosome
	drawOldAssembly <- FALSE
	if( Glazerize & !is.element(i, c("chrVIIpitx1", "chrUn", "chrM")) ){
		chrNumeric <- g$chrname2numeric( i )
		# grab those rows for which both old and new are on the given chrname
		z <- x[x$NewChr == chrNumeric & x$OldChr == chrNumeric,] 
		zstart <- z$NewStart/10^6
		zend <- z$NewEnd/10^6
		drawOldAssembly <- TRUE
		}

	# Read all the blockstats files into a list
	blockstatsList <- vector("list", length=npairs)
	names(blockstatsList) <- sapply(interestingPairs, paste, collapse = ".")
	for(k in 1:npairs){
		# k <- 1
		blockstatsfile 	<- paste(project, i, names(blockstatsList)[k], "blockstats", stepsize, "rdd", sep = ".")
		load(blockstatsfile) # object name is blockstats
		blockstatsList[[k]] <- blockstats
		}
	rm(blockstats)
		
	# Another list for the sliding window results
	slideWinList <- vector("list", length=npairs)
	names(slideWinList) <- sapply(interestingPairs, paste, collapse = ".")

	if(tolower(method) == "vara"){
		for(k in 1:npairs){
			# k <- 1
			FSTwin <- g$slidewin(blockstatsList[[k]], method = c("FST"), 
						nsteps.per.window = nsteps.per.window, windowNmin = windowNmin)
			# Need to convert raw variances to *per-base*
			VarPerBase <- FSTwin$totVARa/FSTwin$nbases
			slideWinList[[k]] <- VarPerBase	
			}
		ibaseMillions <- FSTwin$ibase/10^6
		z <- do.call("cbind.data.frame", slideWinList)
		meanVarPerBase <- apply(z, 1, mean, na.rm = TRUE)
		meanVarPerBase[is.nan(meanVarPerBase)] <- NA
		ylim = range(meanVarPerBase, na.rm=TRUE)
		if(ymax > 0){
			# meanVarPerBase[meanVarPerBase > ymax] <- ymax
			ylim[2] <- ymax
			}
		header <- paste(c(i, "   /    ", ymaxHeader), collapse = " ")
		plot(meanVarPerBase ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
		if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(meanVarPerBase, na.rm=TRUE), col = "blue", lwd = 2)
		} else

	if(tolower(method) == "fst"){
		for(k in 1:npairs){
			# k <- 1
			FSTwin <- g$slidewin(blockstatsList[[k]], method = c("FST"), 
						nsteps.per.window = nsteps.per.window, windowNmin = windowNmin)
			slideWinList[[k]] <- FSTwin$fst
			}
		ibaseMillions <- FSTwin$ibase/10^6
		z <- do.call("cbind.data.frame", slideWinList)
		meanFst <- apply(z, 1, mean, na.rm = TRUE)
		meanFst[is.nan(meanFst)] <- NA
		ylim = range(meanFst, na.rm=TRUE)
		if(ymax > 0){
			ylim[2] <- ymax
			}
		plot(meanFst ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
		if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(FSTwin$fst, na.rm=TRUE), col = "blue", lwd = 2)
		} else

	if(tolower(method) == "css"){
		for(k in 1:npairs){
			# k <- 1
			CSSwin <- g$slidewin(blockstatsList[[k]], method = c("CSS"), 
						nsteps.per.window = nsteps.per.window, windowNmin = windowNmin)
			# Need to convert to *per-base*
			cssPerBase <- CSSwin$CSS/CSSwin$nbases
			slideWinList[[k]] <- VarPerBase	
			}			
		ibaseMillions <- CSSwin$ibase/10^6
		z <- do.call("cbind.data.frame", slideWinList)
		meanCssPerBase <- apply(z, 1, mean, na.rm = TRUE)
		meanCssPerBase[is.nan(meanCssPerBase)] <- NA
		ylim = range(meanCssPerBase, na.rm=TRUE)
		if(ymax > 0){
			# meanCssPerBase[meanCssPerBase > ymax] <- ymax
			ylim[2] <- ymax
			}
		header <- paste(c(i, "   /    ", ymaxHeader), collapse = " ")
		plot(meanCssPerBase ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
		if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(meanVarPerBase, na.rm=TRUE), col = "blue", lwd = 2)
		} else stop("Method must be vara, fst, or css")
	
	}
dev.off()
