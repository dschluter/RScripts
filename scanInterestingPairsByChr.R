#!/usr/bin/Rscript

# Carries out sliding window analysis and plots the results to a pdf file, one chromosome per page.
# All interestingPairs are plotted on the same page, whose length is adjusted according to their number

# interestingPairs must be a list, eg interestingPairs <- list(c("paxl", "paxb"), c("pril", "prib"))
# chromosomes must be a list, eg chromosomes <- list(c("chrIV", "chrXXI"))
# params must be a list, eg params <- list( stepsize = 500, nsteps.per.window = 5, windowNmin = 100 )

# methods: 
# 	"vara" is totVarA, variance among, per base
# 	"fst" is weir-cockerham Fst
# 	"css" is my slightly modified version of felicity's CSS score, per base

# stepsize is the size of the block in the corresponding blockstats file
# nsteps.per.window 	# window size is (nsteps.per.window)*(stepsize), e.g., 5*500 = 2500
# windowNmin	 		# minimum number of good bases in window to include, otherwise not included

# Glazerize = TRUE results in line segments being included in plot to locate the old assembly of that chr
# 		"glazerFileS4 NewScaffoldOrder.csv" must be in local directory

# Uses files named " project.chr.pop1.pop2.blockstats.stepsize.rdd "


# qsub -I -l walltime=04:00:00 -l mem=4gb 
# module load R/3.1.2
# R

# setwd("~/Desktop")
# git("genome.r")

# set ymax to zero to ignore this argument

args <- commandArgs(TRUE) 
# args <- c( "BenlimAllMarine", "species-pairs", "vara", "all", 0)

project <- args[1]
pairs <- args[2]
method <- args[3]
chromosomes <- args[4]
ymax <- args[5]

# these are defaults, not arguments
stepsize = 500
nsteps.per.window = 5  	# window size is (nsteps.per.window)*(stepsize), e.g., 5*500 = 2500
windowNmin = 100 		# windowNmin is minimum number of good bases in window
orderChr <- TRUE
Glazerize <- TRUE
scafFile <- "glazerFileS4 NewScaffoldOrder.csv"

if(pairs == "species-pairs") 
	interestingPairs <- list(
		c("paxl", "paxb"), c("pril", "prib"), c("qryl", "qryb"), c("ensl", "ensb")
		)

if(pairs == "marinepac-pairs") 
	interestingPairs <- list(
		c("marine-pac", "paxb"), c("marine-pac", "prib"), c("marine-pac", "qryb"), c("marine-pac", "ensb"),
		c("marine-pac", "paxl"), c("marine-pac", "pril"), c("marine-pac", "qryl"), c("marine-pac", "ensl")
		)

if(pairs == "solitary-benthic") 
	interestingPairs <- list(
		c("solitary", "paxb"), c("solitary", "prib"), c("solitary", "qryb"), c("solitary", "ensb")
		)

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
		
# Adjust page size according to the number of interestingPairs
npairs <- length(interestingPairs)
pdf( paste(project, pairs, method, "slidewin", "pdf", sep = "."), height = max(11, round(npairs * 1.5)) )
par(mfrow = c(npairs, 1), mar = c(2, 2, 1.5, 1) + 0.1)

# If ymax is set, include in header
ymaxHeader <- ""
if(ymax > 0) ymaxHeader <- paste("(ymax ", ymax, ")", sep = "")

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

	lapply(interestingPairs, function(groupnames){
		# groupnames <- interestingPairs[[1]]
		blockstatsfile 	<- paste(project, i, paste(groupnames, collapse = "."), "blockstats", stepsize, "rdd", sep = ".")
		load(blockstatsfile) # object name is blockstats
		header <- paste(c(groupnames, "      /      ", i, ymaxHeader), collapse = " ")

		if(tolower(method) == "vara"){
			FSTwin <- g$slidewin(blockstats, method = c("FST"), nsteps.per.window = nsteps.per.window, 
				windowNmin = windowNmin)
			# Need to convert raw variances to *per-base*
			VarPerBase <- FSTwin$totVARa/FSTwin$nbases
			ibaseMillions <- FSTwin$ibase/10^6
			ylim = range(VarPerBase, na.rm=TRUE)
			if(ymax > 0){
				VarPerBase[VarPerBase > ymax] <- ymax
				ylim[2] <- ymax
				}
			plot(VarPerBase ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
			if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(VarPerBase, na.rm=TRUE), col = "blue", lwd = 2)
			} else

		if(tolower(method) == "fst"){
			FSTwin <- g$slidewin(blockstats, method = c("FST"), nsteps.per.window = nsteps.per.window, 
				windowNmin = windowNmin)
			Fst <- FSTwin$fst
			ibaseMillions <- FSTwin$ibase/10^6
			ylim = range(Fst, na.rm=TRUE)
			if(ymax > 0){
				Fst[Fst > ymax] <- ymax
				ylim[2] <- ymax
				}
			plot(Fst ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
			if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(FSTwin$fst, na.rm=TRUE), col = "blue", lwd = 2)
			} else

		if(tolower(method) == "css"){
			CSSwin <- g$slidewin(blockstats, method = c("CSS"), nsteps.per.window = nsteps.per.window, 
				windowNmin = windowNmin)
			cssPerBase <- CSSwin$CSS/CSSwin$nbases
			ibaseMillions <- CSSwin$ibase/10^6
			ylim = range(cssPerBase, na.rm=TRUE)
			if(ymax > 0){
				cssPerBase[cssPerBase > ymax] <- ymax
				ylim[2] <- ymax
				}
			plot(cssPerBase ~ ibaseMillions, type="l", lwd = 0.5, main = header, ylim = ylim)
			if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(cssPerBase, na.rm=TRUE), col = "blue", lwd = 2)
			} else stop("Method must be vara, fst, or css")
				
		})
	}
dev.off()
