#!/usr/bin/Rscript

# Averages sliding window results of multiple pairwise genome scans and plots the results to a pdf file.
# Separate sliding windows are carried out for each pair, and then an unweighted average is calculated.

# Carries out sliding window analysis and plots the results to a pdf file, one chromosome per page.
# All interestingPairs are plotted on the same page, whose length is adjusted according to their number
# Unweighted average of multiple sliding window scans is provided at the bottom of each page.

# Analyses files named " project.chr.pop1.pop2.blockstats.stepsize.rdd "

# qsub -I -l walltime=04:00:00 -l mem=4gb 
# module load R/3.1.2
# R

# setwd("~/Desktop")
# git("genome.r")

# set ymax to zero to ignore this argument

chrname 		<- NULL # "chrXXI" or "all" (must be a single chromosome or all chromosomes)
project 		<- NULL
pairtype		<- NULL # "speciesPairs", "marinepacPairs" or "solitaryBenthic"
method			<- NULL # "vara" is totVarA, variance among, per base; 
						# "fst" is weir-cockerham Fst
						# "css" is my slightly modified version of felicity's CSS score, per base
stepsize 		<- 500
nsteps.per.window <- 5  # window size is (nsteps.per.window)*(stepsize), e.g., 5*500 = 2500
windowNmin 		<- 100 	# windowNmin is minimum number of good bases in window
ymax	  		<-   0  # maximum y-value in plots; 0 means not specified
genomeDir		<- NULL # where genomes and "glazerFileS4 NewScaffoldOrder.csv" are located

# Defaults (can't be overridden yet)
orderChr  <- TRUE
Glazerize <- TRUE
scafFile  <- "glazerFileS4 NewScaffoldOrder.csv"

args <- commandArgs(TRUE) 
# args <- c("chrname=chrM","project=Benlim", "pairtype=speciesPairs", "method=vara", "stepsize=500", "nsteps.per.window=5", "windowNmin=100", "ymax=0", "genomeDir=~/tmp")

# Parses the args into a data frame with two columns and then assigns variables 
x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){assign(x[i,1], x[i,2])}
# x
                 # V1            V2
# 1           chrname          chrM
# 2           project        Benlim
# 3          pairtype speciesPairs
# 4            method          vara
# 5          stepsize           500
# 6 nsteps.per.window             5
# 7        windowNmin           100
# 8              ymax             0
# 9         genomeDir         ~/tmp

if(is.null(chrname)) stop("Provide chrname= in arguments")
if(is.null(project)) stop("Provide project= in arguments")
if(is.null(pairtype)) stop("Provide pairtype= in arguments")
if(is.null(method)) stop("Provide method= in arguments")
if(is.null(stepsize)) stop("Provide stepsize= in arguments")
if(is.null(nsteps.per.window)) stop("Provide nsteps.per.window= in arguments")
if(is.null(windowNmin)) stop("Provide windowNmin= in arguments")
if(is.null(ymax)) stop("Provide ymax= in arguments")
if(is.null(genomeDir)) stop("Provide genomeDir= in arguments")

ymax 				<- as.numeric(ymax)
stepwise 			<- as.numeric(stepsize)
nsteps.per.window 	<- as.numeric(nsteps.per.window)
windowNmin 			<- as.numeric(windowNmin)

# get chrnames if chrname="all"
z <- list.files(pattern=glob2rx("*.vcfresultsNew.rdd"))
if(chrname == "all")  chrname <- read.table(text = z, sep = ".", colClasses = "character")$V2

if(pairtype == "speciesPairs") 
	interestingPairs <- list(
		c("paxl", "paxb"), c("pril", "prib"), c("qryl", "qryb"), c("ensl", "ensb")
		)

if(pairtype == "marinepacPairs") 
	interestingPairs <- list(
		c("marine-pac", "paxb"), c("marine-pac", "prib"), c("marine-pac", "qryb"), c("marine-pac", "ensb"),
		c("marine-pac", "paxl"), c("marine-pac", "pril"), c("marine-pac", "qryl"), c("marine-pac", "ensl")
		)

if(pairtype == "marinepacBenthic") 
	interestingPairs <- list(
		c("marine-pac", "paxb"), c("marine-pac", "prib"), c("marine-pac", "qryb"), c("marine-pac", "ensb")
		)

if(pairtype == "marinepacLimnetic") 
	interestingPairs <- list(
		c("marine-pac", "paxl"), c("marine-pac", "pril"), c("marine-pac", "qryl"), c("marine-pac", "ensl")
		)

if(pairtype == "solitaryBenthic") 
	interestingPairs <- list(
		c("paxb", "solitary"), c("prib", "solitary"), c("qryb", "solitary"), c("ensb", "solitary")
		)

npairs <- length(interestingPairs)


# Check that all blockstats files are present
z <- sapply(interestingPairs, paste, collapse = ".")
z <- apply(expand.grid(chrname, z), 1, paste, collapse = ".")
z <- paste(project, z, "blockstats", stepsize, "rdd", sep = ".")
if( !all(file.exists(z)) ){
	print( z[!file.exists(z)] )
	stop("The above blockstats files are missing\n")
	}

if(Glazerize) x <- read.csv(paste(genomeDir, scafFile, sep = "/"))

# Order the chromosomes by their numeric values, leaving out the specially named ones
if(orderChr){
	chrSpecial <- intersect(chrname, c("chrVIIpitx1", "chrUn", "chrM"))
	if(length(chrSpecial) > 0) chrname <- chrname[!is.element(chrname, chrSpecial)]
	chrNumeric <- sapply(chrname, g$chrname2numeric)
	chrname <- chrname[ order(chrNumeric) ]
	if(length(chrSpecial) > 0) chrname <- c(chrname, chrSpecial)
	}

# Adjust page size according to the number of interestingPairs (plus 1 for the unweighted mean scan)
pdf( paste(project, pairtype, method, "slidewin", "pdf", sep = "."), 
				height = max(11, round((npairs + 1) * 1.5)) )
par(mfrow = c(npairs+1, 1), mar = c(2, 2, 1.5, 1) + 0.1)

# If ymax is set, include in header
ymaxHeader <- ""
if(ymax > 0) ymaxHeader <- paste("(ymax ", ymax, ")", sep = "")

library(zoo)
if(tolower(method) == "css") library(gdata)

for(i in chrname){
	# i <- "chrM"
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
		
	# Store the sliding window results for each pair of species in another list
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
		
		# Plot individual pairs
		for(k in 1:npairs){
			# k <- 1
			ylim = range(slideWinList[[k]], na.rm=TRUE)
			if(ymax > 0) ylim[2] <- ymax
			header <- paste(c(names(slideWinList)[k], "  /  ", i, ymaxHeader), collapse = " ")
			plot(slideWinList[[k]] ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
			if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(slideWinList[[k]], na.rm=TRUE), 
					col = "blue", lwd = 2)
			}
		# plot average scan
		ylim = range(meanVarPerBase, na.rm=TRUE)
		if(ymax > 0) ylim[2] <- ymax
		plot(meanVarPerBase ~ ibaseMillions, type="l", lwd = 0.5, col = "red", main = "Average scan", ylim = ylim)
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

		# Plot individual pairs
		for(k in 1:npairs){
			# k <- 1
			ylim = range(slideWinList[[k]], na.rm=TRUE)
			if(ymax > 0) ylim[2] <- ymax
			header <- paste(c(names(slideWinList)[k], "  /  ", i, ymaxHeader), collapse = " ")
			plot(slideWinList[[k]] ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
			if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(slideWinList[[k]], na.rm=TRUE), 
					col = "blue", lwd = 2)
			}
		# plot average scan
		ylim = range(meanFst, na.rm=TRUE)
		if(ymax > 0) ylim[2] <- ymax
		plot(meanFst ~ ibaseMillions, type="l", lwd = 0.5, col = "red", main = "Average scan", ylim = ylim)
		if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(meanFst, na.rm=TRUE), col = "blue", lwd = 2)
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
		# Plot individual pairs
		for(k in 1:npairs){
			# k <- 1
			ylim = range(slideWinList[[k]], na.rm=TRUE)
			if(ymax > 0) ylim[2] <- ymax
			header <- paste(c(names(slideWinList)[k], "  /  ", i, ymaxHeader), collapse = " ")
			plot(slideWinList[[k]] ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
			if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(slideWinList[[k]], na.rm=TRUE), 
					col = "blue", lwd = 2)
			}
		# plot average scan
		ylim = range(meanCssPerBase, na.rm=TRUE)
		if(ymax > 0) ylim[2] <- ymax
		header <- paste(c(i, "   /    ", ymaxHeader), collapse = " ")
		plot(meanCssPerBase ~ ibaseMillions, type="l", lwd = 0.5, col = "red", main = "Average scan", ylim = ylim)
		if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(meanCssPerBase, na.rm=TRUE), col = "blue", lwd = 2)
		} else stop("Method must be vara, fst, or css")
	
	}
dev.off()
