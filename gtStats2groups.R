#!/usr/bin/Rscript

# -------------------------------------------------------------------
# Analysis of differences between two groups (e.g., benthic-limnetic)

# This code assumes that you are working only with ONE chromosome
# Note: At most 3 ALT alleles permitted per snp in this version

# Code analyzes all snps, including indels. 
# Need to add step if want true snps only

# Q: How did Jones et al handle indels in the 21 genomes project?
# Answer: they did not call indels

# Obtains stats on one pair of populations
# qsub -I -l walltime=02:00:00 -l mem=4gb # work interactively; use "exit" to exit
# module load R/3.1.2
# R

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrXXI", "paxl", "paxb")

includePfisher 	<- TRUE
includeFst     	<- TRUE
trueSnpOnly		<- FALSE
getGenotypeFreq <- FALSE

project <- args[1]
chrname <- args[2]
groupnames <- args[3:length(args)]
vcfname			<- paste(project, ".", chrname, ".var.vcf", sep="")
vcfresultsfile	<- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")

gtstatsfile 	<- paste(project, chrname, paste(groupnames, collapse = "-"), "rdd", sep = ".")
gtstats <- list()

if(length(groupnames) > 2 ) stop("Provide names of only two groups")

cat("\nLoading vcfresults file\n")
load(file = vcfresultsfile)

# names(vcfresults)
# [1] "groupcodes"        "vcf"               "GT"               
# [4] "nCalledGenotypes"  "nAlt"              "altUsedList"      
# [7] "nAltUsed"          "snpTypeList"       "alleleFreqByGroup"

library(VariantAnnotation)
# library(GenomicFeatures)

# Group codes for thw two groups of interest here
# 1 and 2 will indicate the two groups of interest; all other fish are assigned a code of 0
fishnames <- samples(header(vcfresults$vcf)) # pulls fish names
groupcodes <- rep(0, length(fishnames))		 # initialize
for(i in 1:length(groupcodes)){
	x <- grep(groupnames[i], fishnames, ignore.case = TRUE)
	groupcodes[x] <- i
	}
cat("\ngroupcodes:\n")
print(groupcodes)
# [1] 0 0 0 0 0 0 0 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

# Grab genotypes of the two groups of interest only
groups <- groupcodes[groupcodes > 0]
# groups
 # [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

# Pull out the genotypes and allele frequencies for just the two groups being analyzed
genotypes <- vcfresults$GT[ , groupcodes > 0] 
alleleFreqByGroup <- lapply(vcfresults$alleleFreqByGroup, function(x){x[groupnames,]})

# Delete rare alleles (< 5%) ? do this in "saveSnpsByGroup.R"

GTminFrac <- 2/3
nInd <- as.vector(table(groups)) 	# n individuals genotyped in each group eg 11 11  7
nMin <- floor(GTminFrac*nInd) 		# minimum number required in each group eg 7 7 4
nMin <- mapply(rep(5, length(groupnames)), nMin, FUN = min) # criterion is 5 or nMin, whichever is smaller, eg 5 5 4

# OK, need to change "saveSnpsByGroup.R" to include all 0, 1, 2, 3 alleles using factor to make sure is all is a matrix

alleles.per.fish <- lapply(vcfresults$alleleFreqByGroup, function(x){
	rowSums(x, na.rm = TRUE)[groupnames] # can't rowsums x[groupnames,] because not always a matrix; cbind too slow
	})
alleles.per.locus <- lapply(vcfresults$alleleFreqByGroup, function(x){
	colSums(x[groupnames,], na.rm = TRUE)
	})

# -----------------
# Prepare genotypes

if(trueSnpOnly){
	# Remove all snp that are not true snp
	}

# Transpose genotype array - lapply on SNPs easier this way
genotypes <- as.data.frame(t(genotypes), stringsAsFactors = FALSE) 
names(genotypes) <- rownames(vcfresults$GT)

#genotypes[1:22,1:5]
                            # chrXXI:6316_C/T chrXXI:10364_T/A chrXXI:10365_G/T chrXXI:16024_C/T chrXXI:16025_C/G
# PaxBen-PxBmale5-GS11                   <NA>             <NA>             <NA>              0/0              0/0
# PaxBen-PxBmale6-GS12                   <NA>             <NA>             <NA>              0/0              0/0
# PaxBen-PxBmale8-GS10                   <NA>              1/1              1/1              0/0              0/0
# PaxBen-PxCL09femaleBF6-GS13            <NA>             <NA>             <NA>             <NA>             <NA>
# PaxBen-RPxCL09maleBM2-GS9              <NA>              1/1              1/1              0/0              0/0
# PaxLim-PxCL09maleLM1-GS14              <NA>             <NA>             <NA>              0/0              0/0
# PaxLim-PxLfemale6-GS18                  1/1              1/1              1/1              0/0              0/0
# PaxLim-PxLmale102-GS16                 <NA>             <NA>             <NA>              0/0              0/0
# PaxLim-PxLmale106-GS15                 <NA>             <NA>             <NA>              0/0              0/0
# PaxLim-PxLmale107-GS17                 <NA>             <NA>             <NA>              0/0              0/0
# paxb04                                 <NA>             <NA>             <NA>              0/0              0/0
# paxb05                                 <NA>             <NA>             <NA>              0/0              0/0
# paxb06                                 <NA>             <NA>             <NA>              0/0              0/0
# paxb07                                 <NA>             <NA>             <NA>              0/0              0/0
# paxb08                                 <NA>             <NA>             <NA>             <NA>             <NA>
# paxb09                                 <NA>              1/1              1/1              0/0              0/0
# paxl01                                 <NA>             <NA>             <NA>              0/0              0/0
# paxl05                                 <NA>             <NA>             <NA>              0/0              0/0
# paxl09                                 <NA>             <NA>             <NA>              0/0              0/0
# paxl10                                 <NA>             <NA>             <NA>             <NA>             <NA>
# paxl13                                  0/1             <NA>             <NA>              0/0              0/0
# paxl14                                 <NA>             <NA>             <NA>              0/0              0/0

# table(unlist(genotypes))
    # 0/0     0/1     0/2     0/3     1/1     1/2     1/3     2/2     2/3     3/3 
# 4206848  697773   16965     479  984865    4817     141   21180     112     700

# Allele frequency table by group
cat("\nCalculating allele frequencies by group at every snp (maximum of 3 ALT alleles assumed)\n")
z <- split(genotypes, groups)
z1 <- lapply(z[[1]], function(x){ unlist(strsplit(x, split = "/")) }) # split genotypes first group
z2 <- lapply(z[[2]], function(x){ unlist(strsplit(x, split = "/")) }) # split genotypes second group
alleleFreqByGroup <- mapply(z1, z2, FUN = function(x,y){
	z <- table( c( rep(1,length(x)), rep(2,length(y)) ), c(x,y))
	})
# $`chrXXI:6316_C/T`
            
             # 0 1
  # marine-pac 2 0
  # paxb       0 0
  # paxl       1 3
  
# Save key results
gtstats$groupnames <- groupnames 	# "paxl" "paxb"
gtstats$groups <- groups 			# 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
gtstats$genotypes <- genotypes		# columns are loci
gtstats$trueSnpOnly <- trueSnpOnly	# Are results based only on the true snps between the two groups?
gtstats$alleleFreqByGroup <- alleleFreqByGroup
rm(alleleFreqByGroup)

# ------------------
# Genotype frequencies
if(getGenotypeFreq){
	# Genotype frequencies	
	z1 <- outer(c(0,1),c(0,1),paste, sep="/")
	z2 <- outer(c(0,1,2),c(0,1,2),paste, sep="/")
	z3 <- outer(c(0,1,2,3),c(0,1,2,3),paste, sep="/")
	gtypes1 <- z1[upper.tri(z1,diag=TRUE)] # "0/0" "0/1" "1/1"
	gtypes2 <- z2[upper.tri(z2,diag=TRUE)] # "0/0" "0/1" "1/1" "0/2" "1/2" "2/2"
	gtypes3 <- z3[upper.tri(z3,diag=TRUE)] # "0/0" "0/1" "1/1" "0/2" "1/2" "2/2" "0/3" "1/3" "2/3" "3/3" 
	
	# Genotype contingency table for every SNP: genotype by population
	# Inside, change gtypes to a factor so that the same columns are present in every table
	
	cat("\nCalculating genotype frequencies at every snp (maximum of 3 ALT alleles assumed)\n")
	genotypeFreqTable <- lapply(genotypes,function(x){
					z <- table(groups, factor(x, levels = gtypes3))
					}) 
					
	# head(genotypeFreqTable) # 0 is REF allele
	# $`1`   
	# pop 0/0 0/1 1/1 0/2 1/2 2/2 0/3 1/3 2/3 3/3
	  # 1   0   1   1   0   0   0   0   0   0   0
	  # 2   0   0   0   0   0   0   0   0   0   0
	# $`2`
	gtstats$genotypeFreqTable <- genotypeFreqTable
	rm(genotypeFreqTable)
	}

# Determine the subset of loci to be analyzed.
# Requires at least 1 allele per population, and  
genoSubset

# ----------
# Check the indexes of the different possible genotypes when genotypes are multi-allelic
# (need to do this just once -- nothing is saved for later use)
# Generate all possible genotypes for multi-allelic SNP 
# Allele 0 always refers to the REF allele, 1-3 to ALT alleles
# So if REF is C and ALT is A, 00 = CC, 01 = CA and 11 = AA  
# 
# z1 <- outer(c(0,1),c(0,1),paste, sep="/")
# z2 <- outer(c(0,1,2),c(0,1,2),paste, sep="/")
# z3 <- outer(c(0,1,2,3),c(0,1,2,3),paste, sep="/")
# gtypes1 <- z1[upper.tri(z1,diag=TRUE)] # "0/0" "0/1" "1/1"
# gtypes2 <- z2[upper.tri(z2,diag=TRUE)] # "0/0" "0/1" "1/1" "0/2" "1/2" "2/2"
# gtypes3 <- z3[upper.tri(z3,diag=TRUE)] # "0/0" "0/1" "1/1" "0/2" "1/2" "2/2" "0/3" "1/3" "2/3" "3/3"


# -------------------------------------------------------------
# Group differences in allele frequencies using Fisher exact test
# Note that one can use logistic regression to get equivalent of odds ratio when there are more than 2x2 categories.

if(includePfisher){ 
	
	pop <- as.character(groups)  # population names
	
	# ---
	  
	# pop 0/0 0/1 1/1 0/2 1/2 2/2 0/3 1/3 2/3 3/3
	  # 1   0   0   1   0   0   0   0   0   0   0
	  # 2   0   0   3   0   0   0   0   0   0   0

	gtstats$genotypeFreqTable <- genotypeFreqTable
	rm(genotypeFreqTable)
	
	# cat("\nExample genotype frequency table (0 is REF allele)\n")
	# print(vcfresults$genotypeFreqTable[1])
	# $`1`  
	# pop 0/0 0/1 1/1 0/2 1/2 2/2 0/3 1/3 2/3 3/3
	  # 1   6   0   0   0   0   0   0   0   0   0
	  # 2   2   3   0   0   0   0   0   0   0   0
		
	# ---
	# Allele count - a list with tables of the number of 0, 1, 2, and 3 alleles
	# First get the matrix needed to convert genotypes to counts
	
	cat("\nCalculating allele frequencies at every snp (maximum of 3 ALT alleles assumed)\n")
	z <- strsplit(gtypes3, split="/")
	xmat <- cbind(sapply(z,function(x){sum(x=="0")}), sapply(z,function(x){sum(x=="1")}), 
				sapply(z,function(x){sum(x=="2")}), sapply(z,function(x){sum(x=="3")}))

	alleleFreqTable <- lapply(gtstats$genotypeFreqTable, function(x){
		z <- as.data.frame(x %*% xmat)
		names(z) <- c("n0","n1","n2","n3")
		return(z)
		})		
	gtstats$alleleFreqTable <- alleleFreqTable
	rm(alleleFreqTable)
	# print(gtstats$alleleFreqTable[6])
	# $`6`
	  # n0 n1 n2 n3
	# 1 17  3  0  0
	# 2 12  6  0  0
	
	# Total number of alleles present in each species at snp (row sums of alleleFreqTable)
	allele.popsums <- lapply(gtstats$alleleFreqTable, rowSums)
	gtstats$allele.popsums <- allele.popsums
	# allele.popsums[6]
	# $`6`
	 # 1  2 
	# 20 18
	
	# Indicate those alleles actually used in this pair of populations, or whether fixed (colSums of alleleFreqTable)
	allele.colsums <- lapply(gtstats$alleleFreqTable, colSums)
	gtstats$allele.colsums <- allele.colsums
	# allele.colsums[6]
	# $`6`
	# n0 n1 n2 n3 
	# 29  9  0  0

	fixed <- sapply(allele.colsums, function(x){length(x[x>0]) < 2})
	gtstats$fixed <- fixed
	
	# Confirming that fixed is always false when nAltUsed > 1 
	# ** not any more, because nAltUsed is based on *all fish* analyzed, whereas fixed is here based on just 2 groups
	# table(gtstats$fixed, nAltUsed, useNA = "always")
	      # nAltUsed
	             # 0      1      2      3   <NA>
	  # FALSE      0 172664  10428    674      0
	  # TRUE   26102  70751    970     18      0
	  # <NA>       0      0      0      0      0	


	# ---
	# Fisher exact tests of allele and genotype frequencies
	# Using ALLELE frequencies
	
	cat("\nCalculating Fisher exact test p-values for group allele freq differences (maximum of 3 ALT alleles assumed)\n")
	pfisher <- mapply(gtstats$alleleFreqTable, gtstats$fixed, gtstats$allele.colsums, FUN = function(x,y,z){
		# i <- 1; x = gtstats$alleleFreqTable[[i]]; y = gtstats$fixed[[i]]; z = gtstats$allele.colsums[[i]]
		if(y){pfisher <- NA}
		else{
			tab <- x[ , z > 0]
			# print(tab)
			# nalleles <- ncol(tab)
			# print(nalleles)
			pfisher <- fisher.test(as.matrix(tab))$p.value
			}
		return(pfisher)
		})

	gtstats$pfisher <- pfisher
	rm(pfisher)
	
	# Write to a file in bedGraph format
	# g$write.bedGraph(paste(chrname,".pfisher.bedGraph", sep=""), chrinfo = chrname, vcfresults$POS, vcfresults$POS+1, round(-log10(pfisher),3),
		# graphname = "-log10Pfisher", col = "41,155,241", altcol = "100,100,100", ylinemark = 4,
		# description = "-log10Pfisher") 

	# Plot where the significant SNPs reside
	
	# pdf(file = paste(vcfdir, project, ".", chrname, ".pfisherPlot", ".pdf", sep=""))	
	# x <- start(vcfresults$rowdata)/(10^6)
	# y <- (-log10(vcfresults$pfisher) + rnorm(length(x), mean=0, sd=.005))
	# colour <- (round(-log10(vcfresults$pfisher)/1.2)+1)
	# plot( y ~ x, pch = ".", ylab="-log10 Pfisher", col = colour, xlab = chrname,
		# main = "Fisher exact test group allele freq differences (max 3 ALT alleles)") 
	# dev.off()
	
	cat("\nTransition-transversion ratio for pfisher < 0.01 (uses first ALT allele if more than one)\n")
	#snptype <- sapply(vcfresults$snptypeList, function(x){x[1]}) # first ALT allele: snp ins or del?

	snp <- unlist(vcfresults$snpTypeList[nAltUsed > 0])
	ref <- rep(as.vector(ref(vcfresults$vcf)), nAltUsed)
	alt <- unlist(vcfresults$altUsedList)

	snp <- snp[!is.na(snp) & snp == "snp"]
	ref <- ref[!is.na(snp) & snp == "snp"]
	alt <- alt[!is.na(snp) & snp == "snp"]
	ref <- substr(ref, 1, 1) # keep the first letter of REF
	alt <- substr(alt, 1, 1) # keep the first letter of ALT
	print(g$tstv(ref[ref != alt], alt[ref != alt]))

	REF <- as.character(vcfresults$rowdata$REF)
	ALT <- sapply(vcfresults$ALTlist, function(x){x[1]})
	print( g$tstv(REF[snptype == "snp" & vcfresults$pfisher <= 0.01], ALT[snptype == "snp" & vcfresults$pfisher <= 0.01]) )
	# $R
	# [1] 1.258801
	# $tstv
	     # alt
	# ref     pur   pyr
	  # pur 60301 47947
	  # pyr 47864 60306
	# $tot
	# [1] 216418
	}

# cat("\nInterim save\n")
# vcfresults <- vcfresults
# save(vcfresults, file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = ""))
# load(file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = ""))
# vcfresults <- vcfresults

# --------------

# library(Geneland)
# # Try Geneland method - this method below assumes that all genotypes are single-digit
# # Note: It was very fast but only gave a total Fis and Fst, not variance components per locus.
# #	so not as useful as hierfstat
# It also required the installation of a gigantic code development app

# genotypes <- as.data.frame(t(vcfresults$GT), stringsAsFactors = FALSE)
# GTsplit <- lapply(genotypes, function(x){
	# vcfresults <- substr(x, 1, 1)
	# x2 <- substr(x, 3, 3)
	# cbind(vcfresults, x2)
	# })
# GTsplit <- do.call("cbind", GTsplit)
# z <- Fstat(GTsplit, 2, groups)
# z
# $Fis
# [1] 0.02649774 0.05939996

# $Fst
          # [,1]      [,2]
# [1,] 0.0000000 0.4434629
# [2,] 0.4434629 0.0000000
	
# --------------

if(includeFst){

	library(hierfstat)
	
	cat("\nCalculating Fst (slowly) between groups\n")
	
	# Per-locus Fst using Weir and Cockrham estimates
	# Returns NA if there isn't at least 8 alleles for each pop
	# Fst is calculated as tsiga/sum(c(tsiga, tsigb, tsigw)), where
	# tsiga <- sum(lsiga)
	# tsigb <- sum(lsigb)
	# tsigw <- sum(lsigw)
	# lsiga, lsigb, and lsigw are the variance components: for each locus among populations,
	#   among individuals within populations, and within individuals, respectively
	# wc command does not like duplicate row names

	# Uses alleleFreqTable as calculated earlier
	# A snp needs at least 8 alleles for each pop to be analyzed
	# A snp needs >= 2 different alleles across all individuals (ie ALT snp not fixed in L and B) or wc will protest

	#whichsnp <- ( sapply(totalleles.lb, function(x){all(x >= nminalleles)}) & !vcfresults$fixed.lb )

	# sum(whichsnp)
	# [1] 127278  # number of TRUE values

	pop <- as.character(groups)  # population names

	# Convert "GT" to hierfstat format:
	# Columns are genotypes indicated as 11 (0/0), 12 (0/1), 22 (1/1), or NA (NA)

	# This next operation is redundant (did it already above) but can't use the duplicate names
	# Redoing it here so that names remain as "V1" "V2" etc

	# genotypes <- as.data.frame( t(vcfresults$GT), stringsAsFactors=FALSE) 

	#genotypes[,1:25]

	geno <- unlist(genotypes)
	
	geno[geno == "0/0"] <- "11"
	geno[geno == "0/1"] <- "12"
	geno[geno == "1/1"] <- "22"
	geno[geno == "0/2"] <- "13"
	geno[geno == "1/2"] <- "23"
	geno[geno == "2/2"] <- "33"
	geno[geno == "0/3"] <- "14"
	geno[geno == "1/3"] <- "24"
	geno[geno == "2/3"] <- "34"
	geno[geno == "3/3"] <- "44"

	geno <- as.data.frame(matrix(as.integer(geno), nrow = length(pop)))
	#names(geno) <- as.character(1:ncol(genotypes)) # name genotypes by xrow to avoid duplicate names

	#head(names(geno)) 
	#[1] "V1" "V2" "V3" "V4" "V5" "V6"
	#tail(names(geno))
	#[1] "V281602" "V281603" "V281604" "V281605" "V281606" "V281607"

	#table(pop)
	# pop
	 # 1  2 
	# 11 11 
	# Must convert pop to a number -- using factor(pop) means that pops will be in alphabetical order
	xpop <- as.numeric(factor(pop))

	cat("Populations represented:\n", xpop, "\n")
	cat("Individuals corresponding to groups:\n", colnames(vcfresults$GT), "\n")

	# Calculate Fst - # took about 40 min
	# z <- wc(cbind(xpop, geno))

	# Use this slightly faster version - took about 20 minutes
	z <- g$wc.revised(cbind(xpop,geno))

	cat("\nFst\n")
	print(z$FST)
	#[1] 0.4434733
	# head(z$sigma.loc)

	# Check, per locus measures: ** this is the one we're keeping
	# tsiga <- sum(z$sigma.loc[,"lsiga"], na.rm=TRUE)
	# tsigb <- sum(z$sigma.loc[,"lsigb"], na.rm=TRUE)
	# tsigw <- sum(z$sigma.loc[,"lsigw"], na.rm=TRUE)
	# tsiga/sum(c(tsiga, tsigb, tsigw))

	# Another check, using per allele measures:
	# tsiga <- sum(z$sigma[,"siga"], na.rm=TRUE)
	# tsigb <- sum(z$sigma[,"sigb"], na.rm=TRUE)
	# tsigw <- sum(z$sigma[,"sigw"], na.rm=TRUE)
	# tsiga/sum(c(tsiga, tsigb, tsigw))

	z1 <- cbind.data.frame(z$sigma.loc, fst = z$per.loc$FST, fis = z$per.loc$FIS)
	# rownames(z1) <- rownames(vcfresults$GT) # error because of duplicate names
	
	# range of variance among groups at single loci - maximum is 1, same as psd
	# range(z1$fst$lsiga)
	# [1] -0.1666667  1.0000000
	
	vcfresults$fst <- z1
	rm(z1)

	pdf(file = paste(vcfdir, project, ".", chrname, ".FstPlot", ".pdf", sep=""))

	# Plot Fst
	# add a tiny random number to reduce discrete bands
	y <- vcfresults$fst$fst + rnorm(length(start(vcfresults$rowdata)), mean = 0, sd = .005)
	x <- start(vcfresults$rowdata)/(10^6)
	colour <- abs(round(abs(vcfresults$fst$fst) * 4.8)) + 1
	plot( y ~ x, pch = ".", ylab="Fst", col = colour, xlab="MB",
		main = "Fst (max 3 ALT alleles)" ) 

	dev.off()
	
	# ----
	# Transition to transversion ratios for significant snp - snp only (not indels) first ALT allele only
	# Transition-transversion ratios - snp only (not indels) first ALT allele only

	cat("\nTransition-transversion ratio for true snps having Fst >= 0.95 snp (uses first ALT allele if more than one)\n")
	snptype <- sapply(vcfresults$snptypeList, function(x){x[1]}) # first ALT allele: snp ins or del?
	REF <- as.character(vcfresults$rowdata$REF)
	ALT <- sapply(vcfresults$ALTlist, function(x){x[1]})
	print(g$tstv(REF[snptype == "snp" & vcfresults$fst$fst >= .95], ALT[snptype == "snp" & vcfresults$fst$fst >= .95]))
	# $R
	# [1] 1.237113
	# $tstv
	     # alt
	# ref   pur pyr
	  # pur 469 396
	  # pyr 380 491
	# $tot
	# [1] 1736

	# g$tstv(vcfresults$REF[vcfresults$fst$fst >.95],vcfresults$ALT[fst$fst > 0.99])
	
	vcfresults <- vcfresults

	# cat("\nInterim save\n")
	# save(vcfresults, file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = ""))
	# load(file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = ""))
	# vcfresults <- vcfresults

	} # end if(includeFst}
	
# --------------

if(includePsd){
	
	vcfresults <- vcfresults

	cat("\nCalculating 'percent sequence divergence' between each pair of individuals at each variant
		(indels are treated simply as alleles)\n")
	
	psd <- g$psd(vcfresults$GT, groups)
	vcfresults$psd <- psd
	rm(psd)

	vcfresults <- vcfresults
	} # end if(includePsd)
	
# -------------


cat("\nSaving results\n")
save(vcfresults, file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = ""))
# load(file = paste(vcfdir, project, ".", chrname, ".vcfresults.rdd", sep = "")) # vcfresults
# vcfresults <- vcfresults
	
