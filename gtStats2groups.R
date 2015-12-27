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

# setwd("~/Desktop")
# git("genome.r")

# groupnames must uniquely be substrings of the fishnames (ignoring case)
# They are used in a "grep" to divide the fish uniquely into groups

args <- commandArgs(TRUE) # project chrname groupnames[vector]
# args <- c("BenlimPax22pacMar7", "chrUn", "paxl", "paxb")
# args <- c("BenlimAllMarine", "chrXXI", "marine-pac", "qryb")
# args <- c("BenlimAllMarine", "chrXXI", "paxl", "paxb")

Glazerize 		<- TRUE # using new assembly coordinates; use the results in vcfresultsNew files

includePfisher 	<- FALSE
includeFst     	<- TRUE
includePsd		<- TRUE
trueSnpOnly		<- FALSE # Note: this must be evaluated separately for each pair of populations

# note: An indel relative to REF genome, if fixed and shared between limnetic and benthic, 
# is not an indel between limnetic and benthic, but a fixed difference. 
# However, REF indel that is polymorphic within or between limnetics and benthics is an indel.
# So handle indels carefully.
# Q: What if two alt alleles that are indels relative to REF are the same width. 
#	 Might they be true snps between benthic and limnetic? I have not investigated.

project <- args[1]
chrname <- args[2]
groupnames <- args[3:length(args)] # Here, refers only to the two under consideration
if(Glazerize){
	vcfresultsfile <- paste(project, ".", chrname, ".vcfresultsNew.rdd", sep = "")
	} else {
	vcfresultsfile <- paste(project, ".", chrname, ".vcfresults.rdd", sep = "")}

gtstatsfile 	<- paste(project, chrname, paste(groupnames, collapse = "."), "gtstats.rdd", sep = ".")
gtstats <- list()

if(length(groupnames) > 2 ) stop("Provide names of only two groups")

cat("\nProject, chromosome, 2groups:", project, chrname, groupnames, "\n")

cat("\nLoading vcfresults file\n")
load(file = vcfresultsfile)   # object is "vcfresults"
# lapply(vcfresults, object.size)

gc()

# names(vcfresults)
 # [1] "groupnames"        "groupcodes"        "control"          
 # [4] "nInd"              "nMin"              "vcf"              
 # [7] "newChr"            "newPos"            "altUsedList"      
# [10] "snpTypeList"       "alleleFreqByGroup"

control <- vcfresults$control
GTminFrac <- control$GTminFrac
control$gtstatsOptions <- c(includePfisher = includePfisher, includeFst = includeFst, 
	includePsd = includePsd, trueSnpOnly = trueSnpOnly)
control$Glazerize <- Glazerize

library(VariantAnnotation)
# library(GenomicFeatures)

# Group codes for the two groups of interest here
# 1 and 2 will indicate the two groups of interest; all other fish are assigned a code of 0
if(Glazerize){  # grab fish names
	fishnames <- samples(header(vcfresults$vcf[[1]]))
	} else {
	fishnames <- samples(header(vcfresults$vcf))
	}

groupcodes <- rep(0, length(fishnames))		 # initialize
for(i in 1:length(groupcodes)){
	x <- grep(groupnames[i], fishnames, ignore.case = TRUE)
	groupcodes[x] <- i
	}
cat("\ngroupcodes:\n")
print(groupcodes)
 # [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2
# [39] 2 2 2 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 2
# [77] 2 2 2 2 1 1 1 1 1 1
groups <- groupcodes[groupcodes > 0]
# groups
 # [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

# nInd <- as.vector(table(groups)) 	# n individuals genotyped in each group eg 11 11
# nMin <- floor(GTminFrac*nInd) 		# minimum number required in each group eg 7 7
# nMin <- mapply(rep(5, length(groupnames)), nMin, FUN = min) # criterion is 5 or nMin, whichever is smaller, eg 5 5
nInd <- vcfresults$nInd[groupnames]
nMin <- vcfresults$nMin[groupnames]

cat("\nMinimum sample size criterion based on GTminFrac = ", GTminFrac, "*nInd or 5, whichever is smaller\n", sep = "")
print( data.frame(groupnames=groupnames, nInd, nMin) )

# ----------
# Pull out the genotypes for just the two groups being analyzed. 
# 	geno(vcfresults$vcf)$GT has already been edited in "saveSnpsByGroup.R"
# Bring over coding annotations too, if present
if(Glazerize){
	genotypes <- lapply(vcfresults$vcf, function(x){geno(x)$GT[ , groupcodes > 0]})
	genotypes <- do.call("rbind", genotypes)
	} else genotypes <- geno(vcfresults$vcf)$GT[ , groupcodes > 0]

# Wipe out the genotypes in vcf to save memory?
if(Glazerize){
	for(i in 1:length(vcfresults$vcf)) geno(vcfresults$vcf[[i]])$GT <- NULL
	} else geno(vcfresults$vcf)$GT <- NULL

# gcinfo(TRUE)
gc()

cat("\nExtracting allele frequencies for each locus for the two populations being analyzed\n")

alleleFreqByGroup <- lapply(vcfresults$alleleFreqByGroup, function(x){x[groupnames,]})
vcfresults$alleleFreqByGroup <- NULL

cat("\nDone extracting allele frequencies for each locus\n")
gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells 14022358 748.9   25310187 1351.8  25310187 1351.8
# Vcells 57018471 435.1  197794295 1509.1 245339434 1871.8

# ----------
# Drop the rows with insufficient numbers of alleles - these will not be seen again, whether variant or invariant

cat("\nIdentifying loci with sufficient numbers of individuals genotyped\n")

# Identify loci with enough alleles in both groups (the minimum number of alleles is nMin*2)
sufficient.alleles.per.group <- sapply(alleleFreqByGroup, function(x){
	z <- rowSums(x, na.rm = TRUE)
	z[1] >= 2*nMin[1] & z[2] >= 2*nMin[2]
	})

cat("\nDropping loci with insufficient numbers of individuals genotyped\n")

# Drop cases with insufficient numbers of individuals
# Keep snpTypeList too, need it if want to remove non-polymorphic indels later

genotypes <- genotypes[sufficient.alleles.per.group, ]
alleleFreqByGroup <- alleleFreqByGroup[sufficient.alleles.per.group]

# if Glazerize is TRUE, then vcfresults$vcf is a list (containing blocks from 1 or more former chromosomes)
# So need to drop cases with insufficient numbers of individuals separately for each element of vcfresults$vcf
if(Glazerize){
	nPerChr <- sapply(vcfresults$vcf, nrow) # nrow for each vcf block in vcfresults$vcf
	oldChr <- rep(names(nPerChr), nPerChr)
	gtstats$vcf <- vector("list", length(nPerChr))
	names(gtstats$vcf) <- names(vcfresults$vcf)
	for(i in 1:length(vcfresults$vcf)){
		gtstats$vcf[[i]] <- vcfresults$vcf[[i]][ sufficient.alleles.per.group[oldChr == names(gtstats$vcf)[i]] ]
		}
	} else gtstats$vcf <- vcfresults$vcf[sufficient.alleles.per.group]

snpTypeList <- vcfresults$snpTypeList[sufficient.alleles.per.group]
if(Glazerize) newPos <- vcfresults$newPos[sufficient.alleles.per.group]

rm(sufficient.alleles.per.group)
rm(vcfresults) # assuming we have extracted all the useful bits

cat("\nDone dropping loci with insufficient numbers of individuals genotyped\n")

gc()
           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
# Ncells 12626268 674.4   25310187 1351.8  25310187 1351.8
# Vcells 53339500 407.0  158235436 1207.3 245339434 1871.8

# nrow(genotypes)
# [1] 871120

# table(genotypes[rownames(genotypes) == "chrXXI:54638_T/TG",])
# 0/0 1/1 
 # 20   1 

# ---------------------------------
# Determine which loci are variable, ie have more than one allele across both groups.
# Those that are not polymorphic we want to keep as a record of an invariant between these two groups

# "i" = invariant
# "v" = variant
# "n" = missing (later step, the snp is dropped)

status <- rep("v", length(alleleFreqByGroup)) # initialize
is.invariant <- sapply(alleleFreqByGroup, function(x){
	z <- colSums(x, na.rm = TRUE)
	sum( z > 0 ) < 2
	})
status[is.invariant] <- "i"
rm(is.invariant)

# table(status)
     # i      v 
# 627449 243671

# status[rownames(genotypes) == "chrXXI:54638_T/TG"]
# [1] "v"

# Don't drop the invariant sites - these need to be counted with the invariants when doing block stats

# -----------------
# If trueSnpOnly = TRUE, *drop everything except true snp and invariants* 
# Sites have multiple alt alleles, and only some are true snps

# Remember: if an indel w.r.t. REF is fixed for the same allele in both limnetic and benthic, it is not an indel.
# 	so we drop only the polymorphic loci that are classified as indels
if(trueSnpOnly){

	# tGT is transposed genotypes but does NOT include invariant sites between the two groups
	tGT <- as.data.frame(t(genotypes[status == "v", ]), stringsAsFactors = FALSE)
	# ncol(tGT)
	# [1] 153997

	# Identify the ALT alleles that are indels (NAs are ignored)
	which.alt.not.snp <- lapply(snpTypeList[status == "v"], function(x){
		which(x != "snp")
		})
	# test <- c("chrXXI:75572_CACTG/C", "chrXXI:75793_CTTG/C", "chrXXI:76075_T/C")
	# genotypes[test, ]
	# snpTypeList[test]
	# which.alt.not.snp[test]
	
	# which.alt.not.snp["chrXXI:54638_T/TG"]
	# $`chrXXI:54638_T/TG`
	# [1] 1
	
	# Sets tGT individual genotypes to NA if not a true SNP.
	z <- mapply(tGT, which.alt.not.snp, FUN = function(x, i){
		if(sum(i) > 0) x[grep(paste("[",paste(i, collapse = ""),"]", sep = ""), x)] <- NA
		return(x)
		})
	z <- t(z) # is a matrix again, has the right rownames but colnames are absent
	colnames(z) <- rownames(tGT)
	# z[test,]
	# table(z["chrXXI:54638_T/TG",])
	# 0/0 
	 # 20 
	
	genotypes[status == "v", ] <- z
	rm(which.alt.not.snp)
	
	# groupcodes
	 # [1] 0 0 0 0 0 0 0 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
	# groups
	 # [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

	# redo allele frequencies by group - use groups as argument not groupcodes
	alleleFreqByGroup[status == "v"] <- g$tableAlleleFreqByGroup(genotypes[status == "v", ], 
		groupnames, groups) 

	# alleleFreqByGroup["chrXXI:54638_T/TG"]
	
	gc()


	# Filter once more to change the status of rows having insufficient numbers of alleles after indels deleted
	sufficient.alleles.per.group <- sapply(alleleFreqByGroup, function(x){
		z <- rowSums(x, na.rm = TRUE)
		z[1] >= 2*nMin[1] & z[2] >= 2*nMin[2]
		})
	# table(sufficient.alleles.per.group)
	 # FALSE   TRUE 
	 # 12629 256720
	
	status[!sufficient.alleles.per.group] <- "n"
	rm(sufficient.alleles.per.group)
	# table(status)
	     # i      n      v 
	# 115352  12629 141368
	
	# Who is still polymorphic (among those formerly with status "v")
	is.invariant <- sapply(alleleFreqByGroup[status == "v"], function(x){
	z <- colSums(x, na.rm = TRUE)
	sum( z > 0 ) < 2
	})
	status[status == "v"][is.invariant] <- "i"
	rm(is.invariant)

	# table(status)
	# status
	     # i      n      v 
	# 126841  12629 129879
		
	# From now on, analyze only the snp with status = "v"

	} # end if(trueSnpOnly)

gc() 

# ------------------
# Genotype frequencies
# if(getGenotypeFreq){
	# # Generate all possible genotypes for multi-allelic SNP 
	# # Allele 0 always refers to the REF allele, 1-3 to ALT alleles
	# z1 <- outer(c(0,1),c(0,1),paste, sep="/")
	# z2 <- outer(c(0,1,2),c(0,1,2),paste, sep="/")
	# z3 <- outer(c(0,1,2,3),c(0,1,2,3),paste, sep="/")
	# gtypes1 <- z1[upper.tri(z1,diag=TRUE)] # "0/0" "0/1" "1/1"
	# gtypes2 <- z2[upper.tri(z2,diag=TRUE)] # "0/0" "0/1" "1/1" "0/2" "1/2" "2/2"
	# gtypes3 <- z3[upper.tri(z3,diag=TRUE)] # "0/0" "0/1" "1/1" "0/2" "1/2" "2/2" "0/3" "1/3" "2/3" "3/3" 
	
	# # Genotype contingency table for every SNP: genotype by population
	# # Inside, change gtypes to a factor so that the same columns are present in every table
	
	# cat("\nCalculating genotype frequencies at every snp (maximum of 3 ALT alleles assumed)\n")
	# genotypeFreqTable <- lapply(geno,function(x){
					# z <- table(groups, factor(x, levels = gtypes3))
					# }) 
					
	# # head(genotypeFreqTable) # 0 is REF allele
	# # $`1`   
	# # pop 0/0 0/1 1/1 0/2 1/2 2/2 0/3 1/3 2/3 3/3
	  # # 1   0   1   1   0   0   0   0   0   0   0
	  # # 2   0   0   0   0   0   0   0   0   0   0
	# # $`2`
	# gtstats$genotypeFreqTable <- genotypeFreqTable
	# rm(genotypeFreqTable)
	# }
	

# Start saving key results
gtstats$groupnames <- groupnames 	# "paxl" "paxb"
gtstats$groups <- groups 			# 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
gtstats$nInd <- nInd
gtstats$nMin <- nMin
gtstats$control <- control
gtstats$genotypes <- genotypes		# rows are loci
rm(genotypes)
gtstats$alleleFreqByGroup <- alleleFreqByGroup
rm(alleleFreqByGroup)
gtstats$status <- status
rm(status)
if(Glazerize){
	gtstats$newPos <- newPos
	rm(newPos)
	}

# -------------------------------------------------------------
# Group differences in allele frequencies using Fisher exact test
# P = 1 if there are no differences or locus is monomorphic

if(includePfisher){ 
	# Fisher exact tests of allele and genotype frequencies using ALLELE (not genotype) frequencies

	# Set to NA non-polymorphic sites with sufficiently many genotypes
	gtstats$pfisher <- rep(NA, length(gtstats$alleleFreqByGroup))
	names(gtstats$pfisher) <- names(gtstats$alleleFreqByGroup)
	
	cat("\nCalculating Fisher exact test p-values for group allele freq differences (maximum of 3 ALT alleles assumed)\n")
	# Loci that are not divergent will have a pfisher = 1
	gtstats$pfisher[gtstats$status == "v"] <- sapply(gtstats$alleleFreqByGroup[gtstats$status == "v"], function(x){
		pfisher <- fisher.test(x)$p.value
		})
	
	gc()
	
	# Write to a file in bedGraph format
	# g$write.bedGraph(paste(chrname,".pfisher.bedGraph", sep=""), chrinfo = chrname, vcfresults$POS, vcfresults$POS+1, round(-log10(pfisher),3),
		# graphname = "-log10Pfisher", col = "41,155,241", altcol = "100,100,100", ylinemark = 4,
		# description = "-log10Pfisher") 

	# Plot where the significant SNPs reside
	
	# pdf(file = paste(project, ".", chrname, ".pfisherPlot", ".pdf", sep=""))	
	# x <- start(rowData(gtstats$vcf))/(10^6)
	# y <- (-log10(gtstats$pfisher) + rnorm(length(x), mean=0, sd=.005))
	# colour <- (round(-log10(gtstats$pfisher)/1.2)+1)
	# plot( y ~ x, pch = ".", ylab="-log10 Pfisher", col = colour, xlab = chrname,
		# main = "Fisher exact test group allele freq differences (max 3 ALT alleles)") 
	# dev.off()
	
	# cat("\nTransition-transversion ratio for pfisher < 0.01\n")
	}

gc()

# --------------

if(includeFst){
	
	# Prepare genotypes for Fst

	# fst <- data.frame(row.names = rownames(gtstats$genotypes) ) # initialize
	# geno <- as.data.frame(t(gtstats$genotypes[gtstats$status == "v", ]), stringsAsFactors = FALSE) 

	library(hierfstat)
	
	cat("\nCalculating Fst (slowly) between groups\n")
	
	# Per-locus Fst using Weir and Cockrham estimates
	# Fst is calculated as tsiga/sum(c(tsiga, tsigb, tsigw)), where
	# tsiga <- sum(lsiga)
	# tsigb <- sum(lsigb)
	# tsigw <- sum(lsigw)
	# lsiga, lsigb, and lsigw are the variance components: for each locus among populations,
	#   among individuals within populations, and within individuals, respectively
	# wc command does not like duplicate row names

	# groups
	# [1] 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1

	pop <- as.character(groups)

	# Convert genotypes to hierfstat format:
	# Columns are genotypes indicated as 11 (0/0), 12 (0/1), 22 (1/1), or NA (NA)
	# temp <- unlist(geno)
	# temp <- g$recode(temp, c("0/0","0/1","1/1","0/2","1/2","2/2","0/3","1/3","2/3","3/3"), 
							# c("11","12", "22", "13", "23", "33", "14", "24", "34", "44"))
	# temp <- as.data.frame(matrix(as.integer(temp), nrow = length(pop)))
	# colnames(temp) <- colnames(geno) # otherwise names are lost

	# head(names(temp)) 

	#table(pop)
	# pop
	 # 1  2 
	# 11 11 

	# rm(geno)
	
	# system.time( z <- g$wc.revised(cbind(groups, temp[,1:1000])) ) # 5 sec

	# Must use groups instead of pop in wc command because must be a number

	# Calculate Fst - # took about 40 min
	# z <- wc(cbind(groups, temp))

	# Use this slightly faster version - took about 20 minutes
	# z <- g$wc.revised(cbind(groups,temp))

	# ---
	# Alternate strategy:
	# Break up into smaller blocks to reduce memory needed by hierfstat

	# Initialize
	fst <- matrix(nrow=nrow(gtstats$genotypes), ncol = 5, 
			dimnames = list(rownames(gtstats$genotypes), c("lsiga","lsigb","lsigw","fst","fis")) ) # initialize

	geno <- gtstats$genotypes[gtstats$status == "v", ]
	# geno <- geno[1:1000,]

	blocks <- cut(1:nrow(geno), 100) # break up the cases into many smaller blocks and run on each
	# z <- by(geno, blocks, FUN=function(x){
	# system.time({
	z <- list()
	for(i in 1:length(levels(blocks))){
			# i <- 1
			x <- t( geno[blocks == levels(blocks)[i], ] )
			x1 <- apply(x, 2, function(x){
				as.integer( g$recode(x, c("0/0","0/1","1/1","0/2","1/2","2/2","0/3","1/3","2/3","3/3"), 
								c("11","12", "22", "13", "23", "33", "14", "24", "34", "44")) )
								})
			x1 <- as.data.frame(x1, stringsAsFactors = FALSE)
			# Must use groups in wc command because must be a number
			z1 <- g$wc.revised( cbind(groups, x1) )
			z[[i]] <- cbind( lsiga  = z1$sigma.loc[,"lsiga"], 
						lsigb = z1$sigma.loc[,"lsigb"], 
						lsigw  = z1$sigma.loc[,"lsigw"], 
						fst = z1$per.loc$FST, 
						fis = z1$per.loc$FIS )
			}
			# }) # 27 minutes
			# })
	# }
	z <- do.call("rbind", z)
	fst[gtstats$status == "v", ] <- z
	rm(geno)
	rm(z)
	
	gc() # when trueSnpOnly = FALSE
	          # used  (Mb) gc trigger   (Mb)  max used   (Mb)
	# Ncells 13134522 701.5   25310187 1351.8  25310187 1351.8
	# Vcells 61577380 469.8  158235436 1207.3 245339434 1871.8

	cat("\nWhole-chromosome Fst:\n")
	# print(z$FST)
	#[1] 0.4119456
	
	z <- colSums(fst, na.rm = TRUE)
	FST <- z["lsiga"]/sum(c(z["lsiga"], z["lsigb"], z["lsigw"]))
	print(unname(FST))
	#[1] 0.4119456
	
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
	
	# head(fst, 10)
	                          # lsiga        lsigb      lsigw         fst         fis
	# chrUn:5108405_C/A             NA           NA         NA          NA          NA
	# chrUn:5108407_C/G    0.002492877 -0.001282051 0.06666667  0.03672613 -0.01960784
	# chrUn:5108413_C/T             NA           NA         NA          NA          NA
	# chrUn:5108418_C/A             NA           NA         NA          NA          NA
	# chrUn:5108425_C/T             NA           NA         NA          NA          NA
	# chrUn:5108430_G/A             NA           NA         NA          NA          NA
	# chrUn:5109443_ATC/A           NA           NA         NA          NA          NA
	# chrUn:5109445_C/A             NA           NA         NA          NA          NA
	# chrUn:5109559_A/G   -0.015202530  0.036868687 0.30000000 -0.04726183  0.10944528
	# chrUn:5109564_A/G             NA           NA         NA          NA          NA
			
	gtstats$fst <- fst
	rm(fst)
	
	gc()
	           # used  (Mb) gc trigger   (Mb)  max used   (Mb)
	# Ncells 13136759 701.6   25310187 1351.8  25310187 1351.8
	# Vcells 65349107 498.6  158235436 1207.3 245339434 1871.8

	# pdf(file = paste(project, ".", chrname, ".FstPlot", ".pdf", sep=""))
	# Plot Fst
	# add a tiny random number to reduce discrete bands
	# y <- vcfresults$fst$fst + rnorm(length(start(vcfresults$rowdata)), mean = 0, sd = .005)
	# x <- start(vcfresults$rowdata)/(10^6)
	# colour <- abs(round(abs(vcfresults$fst$fst) * 4.8)) + 1
	# plot( y ~ x, pch = ".", ylab="Fst", col = colour, xlab="MB",
		# main = "Fst (max 3 ALT alleles)" ) 

	# dev.off()

	# Interim save
	# save(gtstats, file = gtstatsfile)
	
	} # end if(includeFst}

# load(file = gtstatsfile) # object is "gtstats"

gc()

# --------------

if(includePsd){
	
	cat("\nCalculating 'percent sequence divergence' between each pair of individuals at each variant
		(indels are treated simply as alleles)\n")
	
	z <- g$psd(gtstats$genotypes[gtstats$status == "v", ], groups)
	
	# nrow(z)
	# head(z)
	
	# initialize
	psd <- as.data.frame( matrix(nrow = nrow(gtstats$genotypes), ncol = ncol(z)), stringsAsFactors = FALSE )
	rownames(psd) <- rownames(gtstats$genotypes)
	
	psd[ rownames(z) , ] <- z
	colnames(psd) <- colnames(z) # duplicate names were ok
	rm(z)
	
	cat("\nDone calculating percent sequence divergence\n")
	gc()

	gtstats$psd <- psd
	rm(psd)

	} # end if(includePsd)
	
# -------------

gc()

cat("\nSaving results\n")
save(gtstats, file = gtstatsfile)
# load(file = gtstatsfile) # object is "gtstats"
	
