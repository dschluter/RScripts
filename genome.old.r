#To save memory overall, this now moved to blockStats2groups
# # g$blockstats <- function(gtstats, stepsize, goodInvariants, chrvec, control){	
	
	# # This function calculates preliminary results for sliding window analyses on a (repeat-masked) chromosome
	# # It breaks the chromosome into blocks or bins of size stepsize and calculates a value of interest in each bin

	# # gtstats is the gtstats list containing the variables of interest for a single chromosome
	# # Function asks whether gtstats$fst and gtstats$psd exist and computed blockstats accordingly.

	# # GATK 3.4 no longer gives multiple rows of snps at the same value of POS

	# Glazerize <- control$Glazerize
	# # psdMissingAction <- control$psdMissingAction
	# if( !is.element(psdMissingAction, c("meanAll", "meanBW", "meanGroup")) )
		# stop("psdMissingAction must be one of meanAll, meanBW, or meanGroup")

	# # psdMissingAction is for g$blockstats, how to average when psd values (pairwise seq divergence) are missing.
	# # Must be one of the following if task include "css":
	# # 	"meanAll", then psd = NA replaced by the mean pairwise distance for all non-missing psd values
	# # 	"meanBW", then psd = NA replaced by mean psd, calculated separately for Between-group and Within-group pairs.
	# #   "meanGroup" then psd = NA replaced by mean psd, calculated separately for every unique type of pair
	# #		For example, pairs that are "1,1", "1,2" and "2,2" are treated separately, whereas "meanBW" treats
	# #		"1,1" and "2,2" as the same (i.e., both are within-group).

	# # chrvec is a vector of the whole MASKED chromosome, all nucleotides as distinct elements, obtained by:
	# #       Masked bases are indicates with an "M"

	# # invariants contains numbers indicating chromosome positions meeting the minimum-number-of genotypes criterion
	# #		Probably in file "processedInvariants"

	# # trueSnpOnly = TRUE, only true snps (not indels) were used.
	
# if(Glazerize){
	# goodInvariantsFile <- paste(project, ".", chrname, ".goodInvNew.rdd", sep="")
	# chrvecfile = paste("chrvec.", chrno, ".glazer.rdd", sep = "")
	# } else {
	# goodInvariantsFile <- paste(project, ".", chrname, ".goodInv.rdd", sep="") # object is goodInvariants
	# chrvecfile = paste("chrvec.", chrno, ".masked.rdd", sep = "")
	# }

	# pop <- as.character(gtstats$groups)
	# k <- as.integer(stepsize) # this is the size of the block
	# nbases <- length(chrvec)
		
	# # Establish the break points of the bins into which nucleotides will be grouped (e.g., k = 500 bases per bin)
	# # The last bin goes from the final bin of fully k nucleotides to the last nucleotide. 
	# # This last bin might be small.
	
	# ibase <- c( seq(1, nbases, k), nbases + 1) # bases marking breaks between steps of size k
	# midbase <- ibase + (stepsize)/2 # If want something to indicate midpoint of blocks instead

	# # head(midbase)
	# # [1]  251  751 1251 1751 2251 2751
	# # head(ibase)
	# # [1]    1  501 1001 1501 2001 2501
	# # tail(ibase)
	# # [1] 17355501 17356001 17356501 17357001 17357501 17357773 # * last one is less than k

  # # 1) Count number of non-M bases in each nucleotide bin of stepsize "k". 
	
	# # Break nucleotide index of chr into bins of stepsize k
	# chrbins <- findInterval(1:nbases, ibase) # indicates in which "bin" of size k each base belongs
											 # # if none missing, there should be k 1's, k 2's, etc.
	# # chrbins[499:505]
	# # [1] 1 1 2 2 2 2 2
	
	# # count up the number of unmasked nucleotides in each bin - These are the ones coded as ACGT
	# nUnmasked <- tapply(chrvec, chrbins, function(x){length(x[x %in% c("A","C","G","T")])})
	# # head(nUnmasked)
	# # nUnmasked[1:20]
	 # # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
	 # # 0  0  0  9  0  0  0  0  0  0  0  0  0  0  4  0  0 11  0  0
	# rm(chrvec)

  # # 2) Count up the number of snp and number of invariants in each bin
	# # Need to split SNP counts into same bins of stepsize k and then sum up

	# # Check that this works if not glazerizing (ie gtstats$vcf is a vcf not a list)
	# if(Glazerize) pos <- gtstats$newPos else {
			# library(VariantAnnotation)
			# pos <- start(rowData(gtstats$vcf))
			# }

	# # head(pos) # 
	# # [1] 1794949 1794951 1794957 1794962 1794969 1794974
	
	# # All start positions are unique 
	# #	- snps and indels at the same start POS are just different alleles
	# # length(pos)
	# # [1] 871120
	# # length(unique(pos))
	# # [1] 871120
		
	# # Count up the number of good snp in each bin (cut works because ibase intervals are made as factors)
	# snpBins <- cut(pos, breaks = ibase, right = FALSE)
	# nSnp <- table( snpBins )

	# # head(nSnp)
    # # [1,501)  [501,1001) [1001,1501) [1501,2001) [2001,2501) [2501,3001) 
    # #       0           0           0           0           0           0 
    
    # # Data still includes some invariants, sites not polymorphic in these two populations
    # # The following is now redundant -- I think this information is already in gtstats$status
    # # is.polymorphic <- sapply(gtstats$alleleFreqByGroup, function(x){
		# # z <- colSums(x, na.rm = TRUE)
		# # sum( z > 0 ) >= 2
		# # })
	# # table(is.polymorphic, gtstats$status)
	# # is.polymorphic      i      v
	         # # FALSE 627449      0
	         # # TRUE       0 243671

	# is.monomorphic <- split(gtstats$status == "i", snpBins)
	# # is.monomorphic <- split(!is.polymorphic, snpBins)
	# nMonomorphic <- sapply(is.monomorphic, sum)


	# # Count up the number of good invariants in each bin
	# if(Glazerize){ nInvariants <- table( cut(goodInvariants$newPos, breaks = ibase, right = FALSE) )
	# } else { nInvariants <- table( cut(goodInvariants$POS, breaks = ibase, right = FALSE) ) }
	
	# # head(nInvariants)
    # # [1,501)  [501,1001) [1001,1501) [1501,2001) [2001,2501) [2501,3001)  
    # #       0           0           0           0           0           0 
    
	# rm(goodInvariants)
	# gtstats$genotypes <- NULL
	# gtstats$alleleFreqByGroup <- NULL
	# gtstats$vcf <- NULL

	# # all the following items have the same length

	# results <- data.frame(ibase = ibase[-length(ibase)], midbase = midbase[-length(midbase)], nUnmasked = nUnmasked, 
					# nSnp = as.vector(nSnp) - nMonomorphic, nInvariants = as.vector(nInvariants) + nMonomorphic)

	# # results[110:130,]
	    # # ibase midbase nUnmasked nSnp nInvariants
	# # 110 54501   54751         0    0           0
	# # 111 55001   55251         0    0           0
	# # 112 55501   55751        28    0           0
	# # 113 56001   56251         0    0           0
	# # 114 56501   56751         9    0           0
	# # 115 57001   57251         9    0           0
	# # 116 57501   57751         5    0           0
	# # 117 58001   58251        99    4          60
	# # 118 58501   58751       193   33         160
	# # 119 59001   59251        89    6          83
	# # 120 59501   59751       411   19         389
	# # 121 60001   60251       398    3         194
	# # 122 60501   60751       343    5         222
	# # 123 61001   61251       314    3         266
	# # 124 61501   61751       456    5         180
	# # 125 62001   62251       336    8          88
	# # 126 62501   62751         4    0           0
	# # 127 63001   63251         8    0           0
	# # 128 63501   63751       106    0           0
	# # 129 64001   64251       283    0           0
	# # 130 64501   64751       363    0           6
		
	
  # # 3) Calculate Fst summary stats
	# if( "fst" %in% names(gtstats) ){
		
		# cat("\nBeginning Fst calculations\n")
		# # break the data frame into the bins
		# fstBinned <- split(gtstats$fst, snpBins) 
		
		# # Sum Fst components for each bin
		# z <- lapply(fstBinned, function(z){
				# VARa <- sum(z$lsiga, na.rm = TRUE) # among populations
				# VARb <- sum(z$lsigb, na.rm = TRUE) # among individuals within populations
				# VARw <- sum(z$lsigw, na.rm = TRUE) # within individuals
				# FST <- VARa/(VARa + VARb + VARw)
				# return(c(VARa = VARa, VARb = VARb, VARw = VARw, fst = FST))
				# })
		# # z[1]
		# # $`[1,501)`
		# # VARa VARb VARw  fst 
		   # # 0    0    0  NaN 
		   
		# # z[20000]
		# # $`[9.9995e+06,1e+07)`
		     # # VARa      VARb      VARw       fst 
		# # 2.6595041 0.3818182 2.6363636 0.4684134

		# z <- as.data.frame(do.call("rbind", z))
		
		# # Drop the NaN's
		# z$fst[ is.nan(z$fst) ] <- NA
		
		# # Put results into the results data frame. 
		# results <- cbind.data.frame(results, z)
		# rm(z)
		# rm(fstBinned)

		# cat("\nDone Fst calculations\n")
		
		# } # end fst

  			
  # # 4) Calculate CSS quantities of interest within bins
  
  	# if( "psd" %in% names(gtstats) ){
  		
  		# cat("\nBeginning CSS calculations\n")

		# gtpairs <- names(gtstats$psd)
		# # gtpairs <- combn(pop, 2, FUN=function(x){paste(x, collapse=",")})
		
		# # Eliminate missing pairwise differences by replacing with averages
			
		# # 1. Figure out categories for "psdMissingAction" behavior
		# wbPairs <- rep("w", length(gtpairs))
		# wbPairs[sapply(strsplit(gtpairs, split=","), function(x){length(unique(x))}) == 2] <- "b"
		# # wbPairs
		  # # [1] "w" "w" "w" "w" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w" "b" "b" "b"
		 # # [19] "b" "b" "b" "w" "w" "w" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w" "b"
		 # # [37] "b" "b" "b" "b" "b" "w" "w" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w"
		 # # [55] "b" "b" "b" "b" "b" "b" "w" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w"
		 # # [73] "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w" "b"
		 # # [91] "b" "b" "b" "b" "b" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w"
		# # [109] "w" "w" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w"
		# # [127] "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w" "w" "w" "w" "b" "b" "b"
		# # [145] "b" "b" "b" "w" "w" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w"
		# # [163] "w" "w" "w" "w" "w" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "w" "w"
		# # [181] "b" "b" "b" "b" "b" "b" "w" "w" "w" "b" "b" "b" "b" "b" "b" "w" "w" "b"
		# # [199] "b" "b" "b" "b" "b" "w" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b"
		# # [217] "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w" "w"
	
		# if(psdMissingAction == "meanBW") psdGroups <- wbPairs 
		# else if(psdMissingAction == "meanAll") psdGroups <- rep(1, length(wbPairs)) 
		# else if(psdMissingAction == "meanGroup") psdGroups <- gtpairs 
		# else stop("Unrecognizeable value for 'psdMissingAction'")
		
		# # 2. Assign average pairwise distance to missing pairwise distance values at every marker separately
		# # Function "ave" modified to allow missing values
		# ave.rm <- function(x, ..., FUN = mean){
		    # if(missing(...)) 
		        # x[] <- FUN(x)
		    # else{
		        # g <- interaction(...)
		        # split(x, g) <- lapply(split(x, g), FUN, na.rm = TRUE)
		    	# }
		    # x
			# }
  		# cat("\nAssigning average pairwise distances\n")

		# # This is slightly wasteful because it still includes monomorphic sites
		# # It uses a lot of memory
		# z1 <- apply(gtstats$psd, 1, function(x){
			# z1 <- ave.rm(x, psdGroups) 
			# x[is.na(x)] <- z1[is.na(x)]     # assign averages to missing pairs
			# x
			# })
		# gtstats$psd <- NULL
		# psd <- as.data.frame(t(z1))
		# names(psd) <- gtpairs
		# rm(z1)

		# psdBins <- split(psd, snpBins)
		# rm(psd)
  		# cat("\nEnded assigning average pairwise distances\n")
		
  		# cat("\nDropping monomorphic sites\n")
		# # 3. Drop the monomorphic sites to minimize confusion
		# # Can't get this mapply to work, no reason
		# # z <- mapply(is.monomorphic, psdBins, FUN = function(i, x){
			# # # x <- psdBins[[10000]]; i <- is.monomorphic[[10000]]
			# # z <- x#[!i,]
			# # return(z)
			# # }, SIMPLIFY = FALSE)
		# # This worked instead
		# z <- list()
		# for(i in 1:length(psdBins)){
			# z[[i]] <- psdBins[[i]][!is.monomorphic[[i]],]
			# }
		# psdBins <- z
		# rm(z)
		# #}		

  		# cat("\nSumming the psd's within bins\n")
		# # Sum the psd's within bins or blocks
		# # colSum yields a value of 0 if there's no data! Fix? Or just pay attention to nSnp in results
		# psdSum <- lapply(psdBins, colSums)
		# psdSum <- as.data.frame( do.call("rbind", psdSum), stringsAsFactors = FALSE )

		# # psdSum[120:125,]
		          # # 2,2       2,2       2,2       2,2      2,1       2,1       2,1
		# # 120 4.1061111 3.8134125 5.6071429 4.6071429 6.000000 5.6071429 4.0000000
		# # 121 0.7600000 0.6052842 0.6052842 0.6052842 0.750000 0.5238095 0.5863095
		# # 122 0.3838589 0.2116101 0.2116101 1.3838589 0.225000 0.3678571 0.2250000
		# # 123 0.4913219 2.0363636 0.4913219 0.3286713 1.100000 1.1000000 1.1000000
		# # 124 0.3929746 0.9286889 0.9286889 0.3284585 1.963333 2.5633333 2.5633333
		# # 125 3.5000000 2.8092589 2.8092589 2.5000000 2.000000 3.7000000 2.5000000
		         # # 2,1      2,1       2,2       2,2       2,2        2,2       2,2
		# # 120 6.500000 4.868056 6.5149063 4.2113717 5.6220492 6.00000000 4.6000000
		# # 121 0.750000 0.750000 0.7600000 0.6052842 0.7600000 0.76000000 0.6052842
		# # 122 1.225000 0.725000 0.4747680 0.1207010 0.4747680 0.12070099 1.1207010
		# # 123 1.100000 1.100000 1.0363636 0.4913219 1.3286713 0.03636364 0.3286713
		# # 124 1.563333 1.463333 0.3284585 0.3284585 0.3284585 0.32845850 0.9286889
		# # 125 3.583692 0.500000 3.1202266 2.9264464 3.4200000 1.50000000 3.5490323
		# # ....

		# names(psdSum) <- paste("psd", names(psdSum), sep = ".")

		# results <- cbind.data.frame(results, psdSum, stringsAsFactors = FALSE)
  		# cat("\nEnding CSS calculations\n")
		# }
		# results
# }

g$sam2gvcf <- function(inputfish = "", chrname = "", mem = 4, walltime = 72, recalibrate = TRUE,
		recalPlot = FALSE, Rversion = "3.1.2", GATKversion = "3.4.0", 
		genome = "gasAcu1pitx1new.fa", knownsites = "knownSnpsAllchrPitx1new.vcf", run = TRUE){
	# Takes .sam file all the way to gvcf file for a single fish.
	# Based on "sam2haplotypeCallerWithRecalibration.pbs" -- including arguments, for simplicity
	# RGPU (see below) is fixed as "Benlim-q15"
	# Recalibration is optional.
	# if recalPlot = TRUE, certain R libraries must be installed, see GATK help for AnalyzeCovariates.

	if(inputfish == "") stop("You need to provide inputfish")
	if(chrname != "") stop("Code not working for a single chromosome yet")
	if(genome != "gasAcu1pitx1new.fa") stop("Code currently works only with gasAcu1pitx1new.fa")
	if(knownsites != "knownSnpsAllchrPitx1new.vcf") stop("Code works only with knownSnpsAllchrPitx1new.vcf")

	# Attach date and time to name of pbs file to make unique
	hour <- gsub("[ :]", "-", Sys.time())
	pbsfile <- paste("sam2gvcf-", inputfish, chrname, "-", hour, ".pbs", sep = "")
	outfile <- file(pbsfile, "w")
	
	writeLines(			"#!/bin/bash", outfile)
	writeLines(			"#PBS -S /bin/bash", outfile)
	writeLines(paste(	"#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste(	"#PBS -l mem=", mem, "gb", sep = ""), outfile)

	writeLines(	"\n# pbs file to run sam2gvcf, based on 'sam2haplotypeCallerWithRecalibration.pbs'", outfile)
	writeLines(	  "# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)

	writeLines('
	samfile="${inputfish}.sam"
	outfile="${inputfish}.out"
	RGID="BENLIM.${inputfish}"
	RGLB="${inputfish}.SB"
	RGSM="${inputfish}" 
	RGPL="ILLUMINA" 
	RGPU="Benlim-q15"
	fastafile="gasAcu1pitx1new.fa"
	faifile="gasAcu1pitx1new.fa.fai"
	dictfile="gasAcu1pitx1new.dict"
	intervals="${inputfish}.intervals"
	sortedbam="${inputfish}.sorted.bam"
	realignedbam="${inputfish}.realigned.bam"
	mkdupmetrics="${inputfish}.mkdup.metrics"
	mkdupbam="${inputfish}.mkdup.bam"
	mkdupbai="${inputfish}.mkdup.bai"
	reducedbam="${inputfish}.reduced.bam"
	vcffile="${inputfish}.vcf"
	recalbam="${inputfish}.recal.bam"
	recaltable="${inputfish}.recal.table"
	afterrecaltable="${inputfish}.after.recal.table"
	recalplots="${inputfish}.recal.pdf"
	
	module load samtools	
	', outfile)

	writeLines(paste('module load gatk', GATKversion, sep = "/"), outfile)
	writeLines(paste('module load R', Rversion, sep = "/"), outfile)

	writeLines('
	if [ ! -f "$dictfile" ]; then
	    java -Xmx2g -jar /global/software/picard-tools-1.89/CreateSequenceDictionary.jar R=$fastafile O=$dictfile
	fi
	
	if [ ! -f "$faifile" ]; then
	    samtools faidx $fastafile
	fi
	
	# 1. SortSam
	java -Xmx2g -jar /global/software/picard-tools-1.89/SortSam.jar I=$samfile O=$sortedbam \\
		SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
	
	# 2. MarkDuplicates
	java -jar /global/software/picard-tools-1.89/MarkDuplicates.jar \\
		I=$sortedbam O=$mkdupbam M=$mkdupmetrics \\
		VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=FALSE ASSUME_SORTED=TRUE
	
	# 3. AddOrReplaceReadGroups
	java -Xmx2g -jar /global/software/picard-tools-1.89/AddOrReplaceReadGroups.jar \\
		RGID=$RGID RGLB=$RGLB RGSM=$RGSM RGPL=$RGPL RGPU=$RGPU I=$mkdupbam O=$sortedbam \\
		SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
	
	# 4. RealignerTargetCreator
	gatk.sh -Xmx4g -T RealignerTargetCreator -R $fastafile -I $sortedbam -o $intervals \\
		--allow_potentially_misencoded_quality_scores
	
	# 5. IndelRealigner:
	gatk.sh -Xmx4g -T IndelRealigner -R $fastafile -I $sortedbam -targetIntervals $intervals \\
		-o $realignedbam -LOD 0.4 --allow_potentially_misencoded_quality_scores
	', outfile)
	
	if(recalibrate){
		writeLines('
		# 6. BaseRecalibrator
		gatk.sh -Xmx4g -T BaseRecalibrator -R $fastafile -I $realignedbam \\
			-knownSites knownSnpsAllchrPitx1new.vcf -o $recaltable
		gatk.sh -Xmx4g -T PrintReads -R $fastafile -I $realignedbam -BQSR $recaltable -o $recalbam 
		', outfile)
		
		if(recalPlot){	
			writeLines('
			gatk.sh -Xmx4g -T BaseRecalibrator -R $fastafile -I $realignedbam \\
				-knownSites knownSnpsAllchrPitx1new.vcf -BQSR $recaltable -o $afterrecaltable
			gatk.sh -Xmx4g -T AnalyzeCovariates -R $fastafile -before $recaltable \\
				-after $afterrecaltable -plots $recalplots
			', outfile)
			}
		
		writeLines('
		# 7. HaplotypeCaller
		gatk.sh -Xmx4g -T HaplotypeCaller -R $fastafile -I $recalbam \\
		     --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \\
		      -o $vcffile --allow_potentially_misencoded_quality_scores
		', outfile)
	
		} else{
		writeLines('
		# 7. HaplotypeCaller
		gatk.sh -Xmx4g -T HaplotypeCaller -R $fastafile -I $realignedbam \\
		     --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \\
		      -o $vcffile --allow_potentially_misencoded_quality_scores
		', outfile)
		}

	writeLines('\nexit 0', outfile)

	writeLines('\necho \"Job finished with exit code $? at: \`date\`\"', outfile)
	
	close(outfile)
	
	# Use as "qsub -v inputfish=PRIL17 sam2haplotypeCallerWithRecalibration.pbs
	if(run){
		qsub <- paste("qsub -v inputfish=", inputfish, sep = "")
		system(paste(qsub, pbsfile))
		}
	}

g$makePbsRunRscript <- function(Rscript = "", Rversion = "3.1.2", mem = 2, walltime = 24, run = TRUE){
	# Creates a *.pbs file to run the full Rscript command "Rscript"
	# The Rscript command executes a particular *.R file and provides any needed arguments
	# This version also downloads the .R file from github so that the latest version is available
	# This version does not include pbsfile as an argument, it will make from the .R root
	# This version automatically submits the pbs file to the queue
	# Example: (note that all the arguments to .R script are included in a single quotation
	# 	g$makePbsRunRscript(Rscript = "countInvariantsByGroup.R BenlimPax22pacMar7 chrXXI 1 paxl paxb marine-pac")
	
	if(Rscript == "") stop("You need to provide an Rscript command")
		
	# Remove and repaste initial "Rscript" command if present in case it is missing from the Rscript text submitted
	# Also grab the name of the .R file argument, store in dotRfile
	Rscript <- sub("^[ ]*Rscript[ ]*|[ ]*", "", Rscript, ignore.case = TRUE)
	dotRfile <- strsplit(Rscript, split = " ")[[1]][1]
	dotRfile <- gsub("[.]R", "", dotRfile, ignore.case = TRUE)
	Rscript <- paste("Rscript", Rscript)
	
	# Attach date and time to name of file to make unique
	# pbsfile <- gsub(".pbs$", "", pbsfile)
	hour <- gsub("[ :]","-",Sys.time())
	pbsfile <- paste(dotRfile, "-", hour, ".pbs", sep = "")
	
	outfile <- file(pbsfile, "w")
	writeLines("#!/bin/bash", outfile)
	writeLines("#PBS -S /bin/bash", outfile)
	writeLines(paste("#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste("#PBS -l mem=", mem, "gb", sep = ""), outfile)
	#writeLines("#PBS -l epilogue=epilogue.script", outfile)

	#writeLines(paste("#PBS -N", jobname), outfile)
	writeLines("\n# Generic torque file to run Rscript", outfile)
	writeLines("# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	writeLines("# Run as \"qsub pbsfile.pbs\"\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)

	# writeLines(paste("curl -o " , dotRfile, ".R", 
		# " https://raw.githubusercontent.com/dschluter/RScripts/master/", dotRfile, ".R", sep = ""), outfile)

	writeLines(paste("module load R/", Rversion, sep = ""), outfile)
	writeLines(Rscript, outfile)

	writeLines("\necho \"Job finished with exit code $? at: \`date\`\"", outfile)
	
	close(outfile)
	if(run) system(paste("qsub", pbsfile))
	}
