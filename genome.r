g<-list()

# qsub -I -l walltime=04:00:00 -l mem=4gb 
# module load R/3.1.2
# R

# g$summaryFunction <- function(seqname, bamFile, ...) {
	# # S function that takes a single sequence name (chromosome) as input, 
	# # 	reads in specific information from the BAM file, and calculates coverage over that sequence
	# # Taken from Martin Morgan's "Exercises: An Introduction to Rsamtools, 29-30 July, 2010"
	# # "BAM files can be large, containing more information on more genomic regions than
	# # are of immediate interest or than can fit in memory. The first strategy for dealing
	# # with this is to select, using the what and which arguments to ScanBamParam, just
	# # those portions of the BAM file that are essential to the current analysis, e.g.,
	# # specifying what=c(‘rname’, ’qname’, ’pos’) when wishing to calculate coverage
	# # of ungapped reads.
	# # When selective input of BAM files is still too memory-intensive, the file can
	# # be processed in chunks, with each chunk distilled to the derived information of
	# # interest. Chromosomes will often be the natural chunk to process. For instance,
	# # here we write a summary function that takes a single sequence name (chromosome)
	# # as input, reads in specific information from the BAM file, and calculates
	# # coverage over that sequence."
	# library(Rsamtools)
	# param <- ScanBamParam(what = c("pos", "qwidth"),
		# which = GRanges(seqname, IRanges(1, 1e+07)),
		# flag = scanBamFlag(isUnmappedQuery = FALSE))
	# x <- scanBam(bamFile, ..., param = param)[[1]]
	# coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))
	# }


g$geno2snpgds <- function(geno, pos = NULL, chr = NULL, gdsOutfile = "geno.gds"){
	# Converts a genotype data frame "geno" of the following format into  a snpgds object in SNPRelate
	# 	creating a gds file in the local directory at the same time.
	# "geno" must have the following format, with the rownames as id's and the colnames as marker names
	         # chrI:11963492 chrI:1245655 chrI:14261764 chrI:1549902 chrI:18548172
	# GEF1fem1            AC           AG            AG           CC            GG
	# GEF1fem3            AC           AG            AG           CG            AG
	# GEF1fem4            AC           GG            AG           CC            GG
	# GEF1fem6            AC           AG            AG           CC            AG
	# GEF1fem7            AC           AA            GG           CC            AG
	
	# Load the R packages: gdsfmt and SNPRelate
	library(SNPRelate, quietly = TRUE)
	library(stringr, quietly = TRUE)

	# Ensure it is a data frame
	if(class(geno)!="data.frame") geno <- as.data.frame(geno, stringsAsFactors = FALSE)
	
	# Identify the unique alleles at every marker (at most 2 are allowed)
	alleles <- lapply(geno, function(x){
		genotypes <- unique(x)
		alleles <- unique( unlist( strsplit(genotypes, split = "") ) )
		alleles <- alleles[!is.na(alleles)] # drop the NA
		})
	# head(alleles)
	# $`chrI:11963492`
	# [1] "A" "C"
	# $`chrI:1245655`
	# [1] "A" "G"
	# $`chrI:14261764`
	# [1] "A" "G"
	# $`chrI:1549902`
	# [1] "C" "G"
	# $`chrI:18548172`
	# [1] "G" "A"
	# $`chrI:20584613`
	# [1] "A" "C"
	
	# Check that there are exactly 2 alleles for each genotype
	nalleles <- sapply(alleles, length)
	bad <- which(nalleles != 2)
	if(length(bad) > 0){
		geno <- geno[,-bad]
		alleles <- alleles[-bad]
		if(!is.null(pos)) pos <- pos[-bad] 
		if(!is.null(chr)) chr <- chr[-bad] 
		chr = NULL
		cat(length(bad), "markers lacking the requisite 2 alleles were dropped")
		}
	
	# extract the id and marker names (do this after the check for exactly 2 alleles, in case somem markers are dropped)
	id <- rownames(geno)
	snpname <- colnames(geno)

	# Need to convert genotypes to 0, 1, 2 indicating number of alleles
	genoscore <- mapply(geno, alleles, FUN = function(x,y){
		z2 <- str_count(x, y[2])
		})

	# Need to transpose so that markers are rows and individuals are columns
	# result will be a matrix (required for next step)
	genoscore <- t(genoscore) 
	
	# create a variable retaining the original genotypes
	snpAllele <- sapply(alleles, function(x){paste(x, collapse = "/")})

	# Put it all into a list
	genoList <- list(sample.id = id, snp.id = snpname, snp.position = pos, snp.chromosome = chr, 
				snp.allele = snpAllele, genotype = genoscore)
	
	# Create a GDS file of genotypes from the geno matrix
	with(genoList, snpgdsCreateGeno(gdsOutfile, genmat=genotype, sample.id=sample.id, 
			snp.id=snp.id, snp.chromosome=snp.chromosome, snp.position=snp.position, 
			snp.allele=snp.allele, snpfirstdim=TRUE))

	# Open the GDS file
	# genofile <- snpgdsOpen(gdsOutfile)
	# 
	# genofile
	
	invisible()
	
	}
	
g$genotypeGVCFs <- function(gvcffiles, outvcfname, GATKversion = "3.4.0", 
	mem = 4, walltime = 24, genome = "gasAcu1pitx1new.fa", maxAltAlleles = 3, run = TRUE){
	# Generates "genotypeGVCFs.pbs" to carry out genotypeGVCFs to call snps
	#	from multiple gvcf files inputted, for a SINGLE chromosome
	# chrname is extracted from gvcf file names.
	# The fasta file is now a whole genome, not "chromosome.fa", because
	#	we are splitting only after haplotypeCaller, which needs the whole genome fasta
	# Includes the '--intervals chr??" to prevent the whole genome being represented in the gvcf file.
	
	if( !all(grepl("[.]vcf$", gvcffiles)) ) stop("Provide only .vcf files as arguments")

	cat("\n", length(gvcffiles), " gvcf files included:\n", sep = "")
	for(i in 1:length(gvcffiles)){
		cat(gvcffiles[i],"\n")
		}

	# Identify the chromosome number from gvcffiles names
	chrname <- gsub(".*[.](chr[A-z1]+)[.]vcf", "\\1", gvcffiles)
	chrname <- unique(chrname)
	if(length(chrname) != 1) stop("chr missing from some filenames, or more than one chromosome represented")
	

	# Attach date and time to name of pbs file to make unique
	hour <- gsub("[ :]", "-", Sys.time())
	pbsfile <- paste("genotypeGVCFs-", chrname, "-", hour, ".pbs", sep = "")
	outfile <- file(pbsfile, "w")

	writeLines(			"#!/bin/bash", outfile)
	writeLines(			"#PBS -S /bin/bash", outfile)
	writeLines(paste(	"#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste(	"#PBS -l mem=", mem, "gb", sep = ""), outfile)

	writeLines(	"\n# pbs file to run GATK GenotypeGVCFs on multiple gvcf files for one chromosome", outfile)
	writeLines(	paste("\n# pbsfile =", pbsfile), outfile)
	writeLines(	  "\n# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)
		
	writeLines(paste('module load gatk', GATKversion, sep = "/"), outfile)
	
	writeLines(paste("gatk.sh -Xmx", mem, "g -R ", genome, " -T GenotypeGVCFs \\", sep = ""), outfile)
	vcfarguments <- gvcffiles
	for(i in 1:length(gvcffiles)){
		gvcffiles[i] <- paste("     --variant", gvcffiles[i], "\\")
		writeLines(gvcffiles[i], outfile)
		}
	writeLines(paste("     --intervals", chrname, "--includeNonVariantSites --max_alternate_alleles",
		maxAltAlleles, "-o", outvcfname), outfile)
	
	close(outfile)
	if(run) system(paste("qsub", pbsfile))
	}


g$fixBaseQualityScores <- function(samfile = "", mem = 4, walltime = 24, GATKversion = "3.4.0",
		samtoolsVersion = "0.1.19", genome = "gasAcu1pitx1new.fa", outputSam = FALSE, run = TRUE){

	# Standalone function to take a sorted sam or bam file and fix the base quality scores using PrintReads
	# Testing as follows led to an error; need to add read group step for this standalone to work.
	# g$fixBaseQualityScores(samfile = "DenmarkBS27File7.sam.sorted.sam", outputSam = TRUE, run = TRUE)
	# "DenmarkBS27File7.sam.sorted.bam is malformed: SAM file doesn't have any read groups defined in the header. 
	#  The GATK no longer supports SAM files without read groups"

	# if outSam = TRUE, result is converted from bam file to a sam file

	if(samfile=="") stop("Provide samfile on input")
	root <- gsub(".[bs]+am$", "", samfile)
	filetype <- casefold( gsub(".*([bs]+am$)", "\\1", samfile) ) # must be bam or sam
	
	# Attach date and time to name of pbs file to make unique
	hour <- gsub("[ :]", "-", Sys.time())
	pbsfile <- paste("fixBaseQualityScores-", samfile, "-", hour, ".pbs", sep = "")
	outfile <- file(pbsfile, "w")
	
	writeLines(			"#!/bin/bash", outfile)
	writeLines(			"#PBS -S /bin/bash", outfile)
	writeLines(paste(	"#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste(	"#PBS -l mem=", mem, "gb", sep = ""), outfile)

	writeLines(	"\n# pbs file to run GATK PrintReads to fix wrong base quality scores", outfile)
	writeLines(	paste("\n# Samfile =", samfile), outfile)
	writeLines(	  "\n# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)

	parameters <- '
		samfile="${root}.sam"
		bamfile="${root}.bam"
		fixedmisencoded="${root}.fixed-misencoded.bam"
		outsam="${root}.fixed-misencoded.sam"
		fastafile="gasAcu1pitx1new.fa"
		'
	if(genome != "gasAcu1pitx1new.fa"){
		 parameters <- gsub("gasAcu1pitx1new.fa", genome, parameters)
		 parameters <- gsub("gasAcu1pitx1new.dict", gsub("[.][fasta]+",".dict",genome), parameters)
		 }

	convertsam2bam <- '
		samtools view -bS -o $bamfile $samfile
		'
	fixqualityscores <- '
		gatk.sh -Xmx4g -T PrintReads -R $fastafile -I $bamfile -o $fixedmisencoded --fix_misencoded_quality_scores
			'
	convertbam2sam <- '
		samtools view -h -o $outsam $fixedmisencoded
		'
	
	writeLines(parameters, outfile)

	writeLines(paste('module load gatk', GATKversion, sep = "/"), outfile)
	writeLines(paste('module load samtools', samtoolsVersion, sep = "/"), outfile)

	if(filetype == "sam") writeLines(convertsam2bam, outfile)
	writeLines(fixqualityscores, outfile)
	if(outputSam) writeLines(convertbam2sam, outfile)
		
	writeLines('\nexit 0', outfile)
	writeLines('\necho \"Job finished with exit code $? at: \`date\`\"', outfile)
	
	close(outfile)
	
	# Run as "qsub -v root=samfilename pbsfile.pbs
	if(run){
		qsub <- paste("qsub -v root=", root, sep = "")
		system(paste(qsub, pbsfile))
		}
	}


g$qualityScoreDistribution <- function(samfile = "", mem = 2, walltime = 24, Rversion = "3.1.2", 
		sorted = FALSE, run = TRUE){
	# Shows the frequency distribution of bas quality scores in a bam or sam file
	# Useful to check which scoring system is in use ( see http://drive5.com/usearch/manual/quality_score.html )
	# Use "sorted = TRUE" if already sorted; otherwise picard's SortSam is run first
	# Output of SortSam is determined by input: sam -> sam or bam -> bam
	# QualityScoreDistribution also accepts sam or bam as input
	
	if(samfile=="") stop("Provide samfile on input")
	root <- gsub(".[bs]+am$", "", samfile)
	filetype <- casefold( gsub(".*([bs]+am$)", "\\1", samfile) ) # must be bam or sam

	hour <- gsub("[ :]", "-", Sys.time())
	pbsfile <- paste("qualScoreDist-", samfile, "-", hour, ".pbs", sep = "")
	outfile <- file(pbsfile, "w")
	
	writeLines(			"#!/bin/bash", outfile)
	writeLines(			"#PBS -S /bin/bash", outfile)
	writeLines(paste(	"#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste(	"#PBS -l mem=", mem, "gb", sep = ""), outfile)

	writeLines(	"\n# pbs file to run Picard's QualityScoreDistribution'", outfile)
	writeLines(	paste("\n# Samfile =", samfile), outfile)
	writeLines(	  "\n# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)
	
	parameters <- '
		samfile="${root}.sam"
		qualscoredist="${root}.basequalscores.txt"
		qualscorechart="${root}.basequalscores.pdf"
		sortedsam="${root}.sorted.sam"
		'
	sortsam <- '
		java -Xmx2g -jar /global/software/picard-tools-1.89/SortSam.jar I=$samfile O=$sortedsam \\
			SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
			'
	qualityscoredistribution <- '
		java -Xmx2g -jar /global/software/picard-tools-1.89/QualityScoreDistribution.jar \\
			I=$sortedsam O=$qualscoredist CHART=$qualscorechart
			'
	if(sorted) parameters <- gsub('sortedsam="${root}.sorted.sam"', 
		'sortedsam="${root}.sam"', parameters, fixed = TRUE)
	if(filetype == "bam") parameters <- gsub(".sam", ".bam", parameters, fixed = TRUE)

	writeLines(parameters, outfile)

	writeLines(paste('module load R', Rversion, sep = "/"), outfile)

	if(!sorted) writeLines(sortsam, outfile)
	writeLines(qualityscoredistribution, outfile)

	writeLines('\nexit 0', outfile)
	writeLines('\necho \"Job finished with exit code $? at: \`date\`\"', outfile)
	
	close(outfile)
	
	# Run as "qsub -v root=Marine-Atl-Denmark-BS27-Fuelner.sorted pbsfile"
	if(run){
		qsub <- paste("qsub -v root=", root, sep = "")
		system(paste(qsub, pbsfile))
		}
	}

g$haplotypeCaller <- function(gatkBamfile = "", mem = 4, walltime = 72, GATKversion = "3.4.0", 
		samtoolsVersion = "0.1.19", run = TRUE){
	# Takes bamfile file -- recal.bam or realigned.bam output of g$gatk --
	#	and makes a gvcf file using gatk haplotype caller.
	# Input file gatkBamfile must be "root.bam" or "root.sam"
	# If samfile is supplied, it is first converted to bam

	if(gatkBamfile == "") stop("You need to provide bamfile (or samfile) name")
	root <- gsub(".[bs]+am$", "", gatkBamfile)
	filetype <- casefold( gsub(".*([bs]+am$)", "\\1", gatkBamfile) ) # must be bam or sam
	if( !(filetype %in% c("sam","bam")) ) stop("You need to provide .bam (or .sam) file")

	# Attach date and time to name of pbs file to make unique
	hour <- gsub("[ :]", "-", Sys.time())
	pbsfile <- paste("haplotypeCaller-", gatkBamfile, "-", hour, ".pbs", sep = "")
	outfile <- file(pbsfile, "w")
	
	writeLines(			"#!/bin/bash", outfile)
	writeLines(			"#PBS -S /bin/bash", outfile)
	writeLines(paste(	"#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste(	"#PBS -l mem=", mem, "gb", sep = ""), outfile)

	writeLines(	"\n# pbs file to run haplotypeCaller on gatk output, bam or sam file", outfile)
	writeLines(	paste("\n# gatkBamfile =", gatkBamfile), outfile)
	writeLines(	  "\n# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)

	parameters <- '
		bamfile="${root}.bam"
		samfile="${root}.sam"
		baifile="${root}.bai"
		fastafile="gasAcu1pitx1new.fa"
		vcffile="${root}.vcf"
		'
				
	writeLines(parameters, outfile)

	writeLines(paste('module load gatk', GATKversion, sep = "/"), outfile)
	writeLines(paste('module load samtools', samtoolsVersion, sep = "/"), outfile)

	convertsam2bam <- '
		samtools view -bS -o $bamfile $samfile
		samtools index -b $bamfile $baifile
		'
	haplotypecaller <- '
		gatk.sh -Xmx4g -T HaplotypeCaller -R $fastafile -I $bamfile \\
		     --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \\
		      -o $vcffile --allow_potentially_misencoded_quality_scores
		      '	
	if(filetype == "sam") 
		writeLines(convertsam2bam, outfile)
	
	writeLines(haplotypecaller, outfile)
		
	writeLines('\nexit 0', outfile)
	writeLines('\necho \"Job finished with exit code $? at: \`date\`\"', outfile)
	
	close(outfile)
	
	# Run as "qsub -v root=Marine-Pac-Salmon-01-Sara.recal.chrM pbsfile
	if(run){
		qsub <- paste("qsub -v root=", root, sep = "")
		system(paste(qsub, pbsfile))
		}
	}


g$bwaMem <- function(inputfish = "", mem = 2, walltime = 24, BWAversion = "0.7.7", 
				fastq1 = "", fastq2 = "", genome = "gasAcu1pitx1new.fa", run = TRUE){
	# Creates pbs file that runs BWA-mem on paired read fastq files
	# inputfish determines the name of the sam output file
	# Include chrname in inputfish name if running on a single chromosome
	# Based on "bwa077mem.pbs"
	
	if(inputfish == "") stop("You need to provide inputfish")
	if(fastq1 == "") stop("You need to provide fastq1 file name")
	if(fastq2 == "") stop("You need to provide fastq2 file name")
	if(genome != "gasAcu1pitx1new.fa") stop("Code currently works only with gasAcu1pitx1new.fa")
	
	samfile <- paste(inputfish, ".sam", sep = "")

	# Attach date and time to name of pbs file to make unique
	hour <- gsub("[ :]", "-", Sys.time())
	pbsfile <- paste("bwaMem-", inputfish, "-", hour, ".pbs", sep = "")
	outfile <- file(pbsfile, "w")

	writeLines(			"#!/bin/bash", outfile)
	writeLines(			"#PBS -S /bin/bash", outfile)
	writeLines(paste(	"#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste(	"#PBS -l mem=", mem, "gb", sep = ""), outfile)

	writeLines(	"\n# pbs file to run bwa mem, based on 'bwa077mem.pbs'", outfile)
	writeLines(	paste("\n# Inputfish =", inputfish), outfile)
	writeLines(	  "\n# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	
	writeLines("cd $PBS_O_WORKDIR", outfile)
	writeLines('\necho \"Current working directory is `pwd`\"',outfile)
	writeLines('\necho \"Starting run at: \`date\`\"', outfile)

	writeLines(paste("module load bwa", BWAversion, sep = "/"), outfile)
	
	bwaCommand <- 'bwa mem -M -t 16'			

	# bwa mem -M -t 16 $genome $fastq1 $fastq2 > $samfile
	writeLines(
		paste(bwaCommand, genome, fastq1, fastq2, ">", samfile),
		outfile)

	# writeLines('\nexit 0', outfile)

	writeLines('\necho \"Job finished with exit code $? at: \`date\`\"', outfile)
	
	close(outfile)
	
	if(run){
		qsub <- paste("qsub", pbsfile)
		system(qsub)
		}
	}

g$gatk <- function(inputfish = "", mem = 4, walltime = 72, recalibrate = TRUE, recalPlot = FALSE, 
		makegvcf = FALSE, bam2sam = FALSE, Rversion = "3.1.2", GATKversion = "3.4.0", samtoolsVersion = "0.1.19",
		genome = "gasAcu1pitx1new.fa", knownsitesvcf = "knownSnpsAllchrPitx1new.vcf", 
		fixBaseQualityScores = FALSE, qualityScoreDistribution = FALSE, run = TRUE){
	# Takes .sam file (output of bwa) through the gatk pipeline to realigned and/or recalibrated bam.
	# If makegvcf = TRUE, program runs HaplotypeCaller on the resulting bam file.
	# Instead, optionally converts final bam to sam to allow (later) splitting the sam file by chromosome.
	# Pipeline is based on "sam2haplotypeCallerWithRecalibration.pbs"
	# RGPU (see below) is fixed as "Benlim-q15"
	# if recalPlot = TRUE, certain R libraries must be installed, see GATK help for AnalyzeCovariates.
	# if fixQualityScores = TRUE, runs PrintReads with --fix_misencoded_quality_scores before BaseRecalibrator
	# 	This sometimes works to stop errors in BaseRecalibrator about misencoded quality scores 
	#	that happen despite "--allow_potentially_misencoded_quality_scores" argument.
	if(inputfish == "") stop("You need to provide inputfish")

	# Attach date and time to name of pbs file to make unique
	hour <- gsub("[ :]", "-", Sys.time())
	pbsfile <- paste("gatk-", inputfish, "-", hour, ".pbs", sep = "")
	outfile <- file(pbsfile, "w")
	
	writeLines(			"#!/bin/bash", outfile)
	writeLines(			"#PBS -S /bin/bash", outfile)
	writeLines(paste(	"#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste(	"#PBS -l mem=", mem, "gb", sep = ""), outfile)

	writeLines(	"\n# pbs file to run gatkpipeline based on 'sam2haplotypeCallerWithRecalibration.pbs'", outfile)
	writeLines(	paste("\n#Inputfish =", inputfish), outfile)
	writeLines(	  "\n# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)

	parameters <- '
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
		knownsitesvcf="knownSnpsAllchrPitx1new.vcf"
		recalbam="${inputfish}.recal.bam"
		recaltable="${inputfish}.recal.table"
		afterrecaltable="${inputfish}.after.recal.table"
		recalplots="${inputfish}.recal.pdf"
		realignedsam="${inputfish}.realigned.sam"
		recalsam="${inputfish}.recal.sam"
		fixedmisencoded="${inputfish}.fixed-misencoded.bam"
		qualscoredist="${inputfish}.qual-score-dist.txt"
		qualscorechart="${inputfish}.qual-score-dist.pdf"
		'
	if(genome != "gasAcu1pitx1new.fa"){
		 parameters <- gsub("gasAcu1pitx1new.fa", genome, parameters)
		 parameters <- gsub("gasAcu1pitx1new.dict", gsub("[.][fasta]+",".dict",genome), parameters)
		 }
	if(knownsitesvcf != "knownSnpsAllchrPitx1new.vcf"){
		parameters <- gsub("knownSnpsAllchrPitx1new.vcf", knownsitesvcf, parameters)
		}
	writeLines(parameters, outfile)

	writeLines(paste('module load samtools', samtoolsVersion, sep = "/"), outfile)
	writeLines(paste('module load gatk', GATKversion, sep = "/"), outfile)
	writeLines(paste('module load R', Rversion, sep = "/"), outfile)

	writeLines('
		if [ ! -f "$dictfile" ]; then
		    java -Xmx2g -jar /global/software/picard-tools-1.89/CreateSequenceDictionary.jar R=$fastafile O=$dictfile
		fi
		
		if [ ! -f "$faifile" ]; then
		    samtools faidx $fastafile
		fi
		', outfile)

	sortsam <- '
		java -Xmx2g -jar /global/software/picard-tools-1.89/SortSam.jar I=$samfile O=$sortedbam \\
			SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
			'
	qualityscoredistribution <- '
		java -Xmx2g -jar /global/software/picard-tools-1.89/QualityScoreDistribution.jar \\
			I=$sortedbam O=$qualscoredist CHART=$qualscorechart
			'
	
	markduplicates <- '
		java -jar /global/software/picard-tools-1.89/MarkDuplicates.jar \\
			I=$sortedbam O=$mkdupbam M=$mkdupmetrics \\
			VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=FALSE ASSUME_SORTED=TRUE
			'
		
	addorreplacereadgroups <- '
		java -Xmx2g -jar /global/software/picard-tools-1.89/AddOrReplaceReadGroups.jar \\
			RGID=$RGID RGLB=$RGLB RGSM=$RGSM RGPL=$RGPL RGPU=$RGPU I=$mkdupbam O=$sortedbam \\
			SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT
			'
	
	realignertargetcreater <- '
		gatk.sh -Xmx4g -T RealignerTargetCreator -R $fastafile -I $sortedbam -o $intervals \\
			--allow_potentially_misencoded_quality_scores
			'
	indelrealigner <- '
		gatk.sh -Xmx4g -T IndelRealigner -R $fastafile -I $sortedbam -targetIntervals $intervals \\
			-o $realignedbam -LOD 0.4 --allow_potentially_misencoded_quality_scores
			'
			
	fixqualityscores <- '
		gatk.sh -Xmx4g -T PrintReads -R $fastafile -I $realignedbam -o $fixedmisencoded --fix_misencoded_quality_scores
		mv $fixedmisencoded $realignedbam
			'
	
	baserecalibrator <- '
		gatk.sh -Xmx4g -T BaseRecalibrator -R $fastafile -I $realignedbam \\
			-knownSites knownSnpsAllchrPitx1new.vcf -o $recaltable \\
			--allow_potentially_misencoded_quality_scores
		gatk.sh -Xmx4g -T PrintReads -R $fastafile -I $realignedbam -BQSR $recaltable -o $recalbam 
			'
			
	analyzecovariates <- '
		gatk.sh -Xmx4g -T BaseRecalibrator -R $fastafile -I $realignedbam \\
			-knownSites $knownsitesvcf -BQSR $recaltable -o $afterrecaltable \\
			--allow_potentially_misencoded_quality_scores
		gatk.sh -Xmx4g -T AnalyzeCovariates -R $fastafile -before $recaltable \\
			-after $afterrecaltable -plots $recalplots		
			'
	# baserecalibratorFix <- '
		# gatk.sh -Xmx4g -T BaseRecalibrator -R $fastafile -I $realignedbam \\
			# -knownSites knownSnpsAllchrPitx1new.vcf -o $recaltable \\
			# --allow_potentially_misencoded_quality_scores --fix_misencoded_quality_scores
		# gatk.sh -Xmx4g -T PrintReads -R $fastafile -I $realignedbam -BQSR $recaltable -o $recalbam 
			# '
	# analyzecovariatesFix <- '
		# gatk.sh -Xmx4g -T BaseRecalibrator -R $fastafile -I $realignedbam \\
			# -knownSites $knownsitesvcf -BQSR $recaltable -o $afterrecaltable \\
			# --allow_potentially_misencoded_quality_scores  --fix_misencoded_quality_scores
		# gatk.sh -Xmx4g -T AnalyzeCovariates -R $fastafile -before $recaltable \\
			# -after $afterrecaltable -plots $recalplots		
			# '
	haplotypecaller <- '
		gatk.sh -Xmx4g -T HaplotypeCaller -R $fastafile -I $recalbam \\
		     --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \\
		      -o $vcffile --allow_potentially_misencoded_quality_scores
		      '
	convertbam2sam <- '
			samtools view -h -o $recalsam $recalbam
			'
		
	writeLines(sortsam, outfile)
	
	if(qualityScoreDistribution) writeLines(qualityscoredistribution, outfile)
	
	writeLines(markduplicates, outfile)
	writeLines(addorreplacereadgroups, outfile)
	writeLines(realignertargetcreater, outfile)
	writeLines(indelrealigner, outfile)
	
	if(fixBaseQualityScores) writeLines(fixqualityscores, outfile)
	
	if(recalibrate){
		writeLines(baserecalibrator, outfile)
		if(recalPlot){	
			writeLines(analyzecovariates, outfile)
			}
		}
	
	if(makegvcf){
		if(recalibrate) writeLines(haplotypecaller, outfile) else{
			haplotypecaller <- gsub("recalbam", "realignedbam", haplotypecaller)
			writeLines(haplotypecaller, outfile)
			}
		}
		
	if(bam2sam){
		if(recalibrate) writeLines(convertbam2sam, outfile) else{
			convertbam2sam <- gsub("recal", "realigned", convertbam2sam)
			writeLines(convertbam2sam, outfile)
			}
		}
		
	writeLines('\nexit 0', outfile)
	writeLines('\necho \"Job finished with exit code $? at: \`date\`\"', outfile)
	
	close(outfile)
	
	# Run as "qsub -v inputfish=PRIL17 sam2haplotypeCallerWithRecalibration.pbs
	if(run){
		qsub <- paste("qsub -v inputfish=", inputfish, sep = "")
		system(paste(qsub, pbsfile))
		}
	}


g$submitRscript <- function(Rscript = "", Rversion = "3.1.2", mem = 2, walltime = 24, run = TRUE){
	# Replaces makePbsRunRscript
	# Creates a *.pbs file to run the full Rscript command "Rscript"
	# Download .R script file from github before using
	# The Rscript command executes a particular *.R file and provides any needed arguments
	# if run=TRUE, automatically submits the pbs file to the queue
	# Example: (note that all the arguments to .R script are included in a single quotation
	# 	g$submitRscript(Rscript = "countInvariantsByGroup.R BenlimPax22pacMar7 chrXXI 1 paxl paxb marine-pac")
	
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
	
	writeLines(			"#!/bin/bash", outfile)
	writeLines(			"#PBS -S /bin/bash", outfile)
	writeLines(paste(	"#PBS -l walltime=", walltime, ":00:00", sep=""), outfile)
	writeLines(paste(	"#PBS -l mem=", mem, "gb", sep = ""), outfile)
	#writeLines(		"#PBS -l epilogue=epilogue.script", outfile)

	#writeLines(paste("#PBS -N", jobname), outfile)
	writeLines(	"\n# Generic torque file to run Rscript", outfile)
	writeLines(paste("\n# ", Rscript), outfile)
	writeLines(	  "\n# THIS FILE IS GENERATED BY R, DO NOT EDIT\n", outfile)
	
	writeLines("\necho \"Starting run at: \`date\`\"", outfile)

	writeLines(paste("module load R/", Rversion, sep = ""), outfile)
	writeLines(Rscript, outfile)

	writeLines("\necho \"Job finished with exit code $? at: \`date\`\"", outfile)
	
	close(outfile)
	if(run) system(paste("qsub", pbsfile))
	}

g$doesThisGenotypeMatchParent <- function(offspring, parent){
	# offspring <- "AG"; parent <- "AA"
	# Tests whether at least one allele in the offspring could have come from this parent
	# "offspring" and parent are vectors of genotypes.
	# doesThisGenotypeMatchParent("AG", "AA") # TRUE
	# doesThisGenotypeMatchParent("GG", "AA") # TRUE
	# doesThisGenotypeMatchParent("AG", as.character(NA)) # FALSE
	# doesThisGenotypeMatchParent(as.character(NA), as.character(NA)) # TRUE - both NA!

	# This is how "outer" lines up the comparisons:
	# outer(c(off1=1,off2=2), c(par1=3,par2=4), FUN=paste)
	     # par1  par2 
	# off1 "1 3" "1 4"
	# off2 "2 3" "2 4"
	#
	offspring <- strsplit(offspring, split = "")
	parent <- strsplit(parent, split = "")
	z <- mapply(offspring, parent, FUN = function(offspring, parent){
		z <- outer( offspring, parent, FUN = is.element)
		# print(z)
		      # [,1]  [,2]
		# [1,]  TRUE  TRUE
		# [2,] FALSE FALSE
		# Just 1 TRUE in the array is needed to say that offspring genotype matches parent
		sum(z) > 0
		})
	return(z)
	}

g$convertScaf2Chr <- function(scafNumber, pos, scafFile = NULL, scafTable = NULL){
	# Converts coordinates for a single scaffold number to 'old' chr coordinates (Jones et al 2012 genome assembly)
	# Function is based on "g$glazerConvertOld2New" below
	
	# scafNumber is a single element string or number, e.g. "186" as given in the first column of "glazerFileS4NewScaffoldOrder.csv"
	# pos is a vector indicating nucleotide positions to be converted.
	# scafFile gives the path and file name to the file 'glazerFileS4NewScaffoldOrder.csv'
	# scafTable is scafFile read into a data frame
	# If scafTable is provided, then this is used, otherwise scafFile is read from "glazerFileS4NewScaffoldOrder.csv"
	# If scafFile NULL, then the file is grabbed from my github.

	# Returns a data frame of [oldChr, oldPos], ie the old coordinates

	if(is.null(scafTable)) {
		if(is.null(scafFile)){
			library(RCurl)
			scafFile <- getURL( "https://raw.githubusercontent.com/dschluter/genomeScripts/master/glazerFileS4NewScaffoldOrder.csv")
			scafTable <- read.csv(text = scafFile, stringsAsFactors = FALSE)
			} else
		scafTable <- read.csv(scafFile, stringsAsFactors = FALSE)
		}
  
	chrNumeric <- as.integer(scafNumber)
		
	# Find those rows in the table that correspond to the scaffold
	chrRows <- which( scafTable$Scaffold == chrNumeric )
	
	# scafTable[chrRows, ]
	
	# To deal with skips, all of which are 1000 bases (probably all NNNN) set all newPos that are skipped to NA
	# So, to initiate
    oldPos <- rep(NA, length(pos)) 
    oldChr <- rep(NA, length(pos))

	# Repeat the coversion using each row of scafTable corresponding to the old chromosome, one at a time
	for(i in chrRows){
		# i <- chrRows[1]
		x <- scafTable[i, ]
		k <- which( pos <= x$Length)
		oldChr[k] <- x$OldChr
		oldPos[k] <- x$OldStart + (pos - 1)
		}
	
	# Put results in a new data frame
	oldCoords <- data.frame(oldChr = oldChr, oldPos = oldPos, stringsAsFactors = FALSE)
	nrowsBeforeNAdrop <- nrow(oldCoords)
	
	# Clean the data frame of NA's
	oldCoords <- na.omit(oldCoords)
	nrowsAfterNAdrop <- nrow(oldCoords)
	
	if(nrowsBeforeNAdrop != nrowsAfterNAdrop) cat("\nSome POS were dropped -- not included in scafTable -- probably all NNNNNN\n")

	return(oldCoords)
	}

g$fasta2vec <- function(fastafile){
	# Reads a single chromosome fast file and converts to a vector
	# z <- g$fasta2vec("chrVIIpitx1.fa")
	chr <- scan(fastafile, what=character())  
	chr <- chr[-1]  # drops the line with ">chr??"
	chr <- paste(chr, collapse="") # combined all the lines into a single word
	chrvec <- strsplit(chr, split="")[[1]] # break apart genome into individual bases
	}
	
g$fasta2list <- function(fastafile){
	# Reads a fast file and converts to a list of vectors, one element per chromosome
	# Alternatives
	# ape:
	#	chrlist <- read.dna(fastafile, format="fasta", as.character=TRUE, as.matrix=FALSE) # returns a list
	# 	write.dna(chrlist, fastafile, format="fasta", colsep = "") # beware of blanks between columns
	# 	chr <- read.FASTA(fastafile) # returns a list of class "DNAbin"
	# sequinr:
	#	chrlist <- read.fasta(fastafile) # returns a list
	#	write.fasta(chrlist, fastafile)
	# Biostrings: (see "chrVIIpitx1.r" for examples of how to use)
	# 	chr <- read.DNAStringSet(fastafile, "fasta") # returns an object of class "DNAStringSet"
	# 	write.XStringSet(chr, file = fastafile)
	library(ape)
	chrlist <- read.dna(fastafile, format="fasta", as.character=TRUE, as.matrix=FALSE)
	}

g$mergePDF <- function(..., file, gsversion = NULL, in.file = NULL) {
	# Found on github
	# Uses ghostscript
    if (is.null(in.file)) {
      in.file <- substitute(...())
    } 
    infiles <- paste(unlist(lapply(in.file, function(y) as.character(y))), 
        collapse = " ")
    if (is.null(gsversion)) {
      gsversion <- names(which(Sys.which(c("gs", "gswin32c", "gswin64c")) != ""))
      if (length(gsversion) == 0) 
        stop("Please install Ghostscript and ensure it is in your PATH")
      if (length(gsversion) > 1)
        stop("More than one Ghostscript executable was found:", 
             paste(gsversion, collapse = " "), 
             ". Please specify which version should be used with the gsversion argument")
    }   
    pre = " -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="
    system(paste(paste(gsversion, pre, file, sep = ""), infiles, collapse = " "))
}

g$plotSlidewinInterestingPairsByChr <- function(project, chrname, method, interestingPairs, ymax = 0,
			stepsize = 500, nsteps.per.window = 5, windowNmin = 100, orderChr = TRUE,
			Glazerize = TRUE, scafFile = "glazerFileS4 NewScaffoldOrder.csv"){
	# Carries out sliding window analysis and plots the results to a pdf file, one chromosome per page.
	# All interestingPairs are plotted on the same page, whose length is adjusted according to their number

	# interestingPairs must be a list, 
	#	eg interestingPairs <- list(c("paxl", "paxb"), c("pril", "prib"))

	# methods: 
	# 	"vara" is totVarA, variance among, per base
	# 	"fst" is weir-cockerham Fst
	# 	"css" is my slightly modified version of felicity's CSS score, per base
	
	# stepsize is the size of the block in the corresponding blockstats file
	# nsteps.per.window 	# window size is (nsteps.per.window)*(stepsize), e.g., 5*500 = 2500
	# windowNmin	 		# minimum number of good bases in window to include, otherwise not included
	
	# Glazerize = TRUE results in line segments being included in plot to locate the old assembly of that chr
	# 		"glazerFileS4 NewScaffoldOrder.csv" must be in local directory

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
	pdf( paste(project, "slidewin", method, "pdf", sep = "."), height = max(11, round(npairs * 1.5)) )
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
				if(ymax > 0) VarPerBase[VarPerBase > ymax] <- ymax
				ibaseMillions <- FSTwin$ibase/10^6
				plot(VarPerBase ~ ibaseMillions, type="l", lwd = 0.5, main = header)
				if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(VarPerBase, na.rm=TRUE), col = "blue", lwd = 2)
				} else
	
			if(tolower(method) == "fst"){
				FSTwin <- g$slidewin(blockstats, method = c("FST"), nsteps.per.window = nsteps.per.window, 
					windowNmin = windowNmin)
				Fst <- FSTwin$fst
				ibaseMillions <- FSTwin$ibase/10^6
				if(ymax > 0) Fst[Fst > ymax] <- ymax
				plot(Fst ~ ibaseMillions, type="l", lwd = 0.5, main = header)
				if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(FSTwin$fst, na.rm=TRUE), col = "blue", lwd = 2)
				} else
	
			if(tolower(method) == "css"){
				CSSwin <- g$slidewin(blockstats, method = c("CSS"), nsteps.per.window = nsteps.per.window, 
					windowNmin = windowNmin)
				cssPerBase <- CSSwin$CSS/CSSwin$nbases
				ibaseMillions <- CSSwin$ibase/10^6
				if(ymax > 0) cssPerBase[cssPerBase > ymax] <- ymax
				plot(cssPerBase ~ ibaseMillions, type="l", lwd = 0.5, main = header)
				if(drawOldAssembly)	segments(x0 = zstart, x1 = zend, y0 = min(cssPerBase, na.rm=TRUE), col = "blue", lwd = 2)
				} else stop("Method must be vara, fst, or css")
					
			})
		}
	dev.off()
	}

g$chrname2numeric <- function(chrname){
	# Convert chromosome name to Glazer code format (a number if it is a roman numeral, character otherwise)
	# Detects whether chrno refers to an actual number by ensuring there are no lower case letters, no U and no M
	# chrno will work too (ie without the "chr" prefix)
	# Modified to allow chrname to be a vector
	
	# g$chrname2numeric("chrXXI")
	# [1] 21
	# g$chrname2numeric("chrUn")
	# [1] "Un"
	# g$chrname2numeric("chrM")
	# [1] "M"
	# g$chrname2numeric("chrVIIpitx1")
	# [1] "VIIpitx1"
	# g$chrname2numeric("I")
	# [1] 1
	# g$chrname2numeric("XX")
	# [1] 20

	chrno <- gsub("^chr|^group", "", chrname)
	chrNumeric <- chrno
	
	# Convert to an actual number if name of chromosome is a roman numeral
	#	(detects by ensuring there are no lower case letters, no U and no M)
	isNumber <- !grepl("[a-zMU]+", chrno)
	chrNumeric[isNumber] <- as.numeric( as.roman( chrNumeric[isNumber] ) )
	# if( !grepl("[a-zMU]+", chrno) ) chrNumeric <- as.numeric( as.roman( chrNumeric ) )
	chrNumeric
	}

g$numeric2chrname <- function(chrNumeric){
	# Convert chromosome number (a number if it is a roman numeral, character otherwise) to roman
	# Detects whether chrNumeric refers to an actual number by testing for characters
	# Adds the "chr" prefix

	# Modified to allow chrNumeric to be a vector
	
	# g$numeric2chrname(4)
	# [1] "chrIV"
	# g$numeric2chrname(21)
	# [1] "chrXXI"
	# g$numeric2chrname("Un")
	# [1] "chrUn"

	chrname <- chrNumeric
	isNumber <- !grepl("[a-zMU]+", chrNumeric)
	chrname[isNumber] <- as.character( as.roman( as.integer( chrname[isNumber] ) ) )

	# if( !grepl("[a-zMU]+", chrNumeric) ) chrname <- as.roman( as.integer(chrNumeric) ) else chrname <- chrNumeric
	chrname <- paste("chr", chrname, sep = "")
	
	chrname
	}

g$glazerConvertOld2New <- function(chrname, pos, scafFile = NULL, scafTable = NULL){
	# Converts coordinates for a single chromosome from 'old' to 'new'
	# 'Old' is the Jones et al 2012 'old' genome assembly 
	# 'New' is the Glazer et al 2015 'new' assembly coordinates
	# Function is extensively modified from file "glazerconvertCoordinate.R" to allow pos to be vector
	
	# Requires access to the "glazerFileS4NewScaffoldOrder.csv"
	# Returns a list of [newChr, newPos]

	# Inputs:
	# chrname is an element string e.g. "chrXXI" or "chrUn"
	# pos is a vector indicating nucleotide positions to be converted.
	# scafFile gives the path and file name to the file 'glazerFileS4NewScaffoldOrder.csv'
	# 	if NULL, then the file is grabbed from my github.

	if(is.null(scafTable)) {
		if(is.null(scafFile)){
			library(RCurl)
			scafFile <- getURL( "https://raw.githubusercontent.com/dschluter/genomeScripts/master/glazerFileS4NewScaffoldOrder.csv")
			scafTable <- read.csv(text = scafFile, stringsAsFactors = FALSE)
			} else
		scafTable <- read.csv(scafFile, stringsAsFactors = FALSE)
		}
  
	# Convert chromosome name to Glazer code format (a number if it is a roman numeral, character otherwise)
	chrno <- gsub("^chr", "", chrname)
	chrNumeric <- chrno
	# Convert to an actual number if name of chromosome is a roman numeral
	#	(detects by ensuring there are no lower case letters, no U and no M)
	if( !grepl("[a-zMU]+", chrno) ) chrNumeric <- as.numeric( as.roman( chrNumeric ) )
        
	translate <- function(pos, startA, startB, endB, orientation){
	    # Translates from one coordinate system to a second
		if(!orientation == 'reverse'){
	    	pos2 = startB + (pos - startA)
	    	} else{
	      	pos2 = endB - (pos - startA)
	    	}
	    	return(pos2)
	  		}
	
	# Find those rows in the table that correspond to the current chromosome
	oldChrRows <- which( as.character(scafTable$OldChr) == as.character(chrNumeric) )
	# [1] 178 180 181 183 305 # chrXXI
	
	# scafTable[oldChrRows, ]
	    # Scaffold  Length NewChr NewStart   NewEnd NewOrientation OldChr OldStart   OldEnd
	# 178      144  297070     21  1494736  1791805        unknown     21        1   297070
	# 180       16 9196662     21  4442219 13638880        forward     21   298071  9494732 * skips 1000 here to next Oldstart
	# 181       43 1902836     21 13639881 15542716        forward     21  9495733 11398568 * here too
	# 183      127  266767     21 16438087 16704853        forward     21 11399569 11666335 * here too
	# 305      257   50152     Un 10615548 10665699        unknown     21 11667336 11717487
	
	# To deal with skips such as the one above, all of which are 1000 bases (probably all NNNN) set all newPos that are skipped to NA
	# So, to initiate
    newPos <- rep(NA, length(pos)) 
    newChr <- rep(NA, length(pos))

	# Repeat the coversion using each row of scafTable corresponding to the old chromosome, one at a time
	for(i in oldChrRows){
		# i <- oldChrRows[1]
		x <- scafTable[i, ]
		k <- which( pos >= x$OldStart & pos <= x$OldEnd)
		newChr[k] <- x$NewChr
		newPos[k] <- translate(pos[k], x$OldStart, x$NewStart, x$NewEnd, x$NewOrientation)
		}
	
	# Put results in a new data frame
	newCoords <- data.frame(newChr = newChr, newPos = newPos, stringsAsFactors = FALSE)
	nrowsBeforeNAdrop <- nrow(newCoords)
	
	# Clean the data frame of NA's
	newCoords <- na.omit(newCoords)
	nrowsAfterNAdrop <- nrow(newCoords)
	
	if(nrowsBeforeNAdrop != nrowsAfterNAdrop) cat("\nSome POS were dropped -- not included in scafTable -- probably all NNNNNN\n")

	return(newCoords)
	}

# This works but is too slow
g$glazerConvertCoordinate <- function(chr, pos, direction = 'old2new', scafFile = "glazerFileS4 NewScaffoldOrder.csv"){
	# Converts Jones et al 2012 'old' genome assembly coordinates to Glazer et al 2015 'new' assembly coordinates
	# 	or the reverse direction ( direction = 'new2old' )
	# Requires access to the "glazerFileS4 NewScaffoldOrder.csv"
	# Command is modified from file "glazerconvertCoordinate.R" to allow chr and pos to be vectors
	# Returns a list of [chromosome, position]

	# Inputs:
	# chr is a number or string, e.g. c(1,'1','Un') or c("IV",'XXI','Un') of the starting chromosome.
	# pos is a number of the starting position.
	# direction is either 'old2new' or 'new2old'.
	# scafFile gives the path and file name to the file 'FileS4 NewScaffoldOrder.csv'
  
	# Examples:
	# g$glazerConvertCoordinate(3, 1538202) 			# same position
	# g$glazerConvertCoordinate('Un', 37499024) 			# now on chr 1
	# g$glazerConvertCoordinate('Un', 23343225) 			# now on chr 2
	# g$glazerConvertCoordinate(c('Un','Un',"1"), c(23343225, 37499024,541084)) # "1" instead of "I"
	# g$glazerConvertCoordinate(c('Un','Un',"I"), c(23343225, 37499024,541084)) # "I" instead of "1"
	# g$glazerConvertCoordinate("1", 541084) 			# different location on chr 1
	# g$glazerConvertCoordinate('1', 680442, 'new2old') # reverse of previous line
	# g$glazerConvertCoordinate(12, 594205, 'new2old') 	# used to be on Un
	# g$glazerConvertCoordinate(1, 540083) 				# in between contigs in original assembly, so NA
    
    if(length(chr) != length(pos)) 
    	warning("In glazerConvertCoordinate, length of chr and pos must be the same")

    # Convert roman to numbers (won't hurt if already numbers)
    # Won't work if both numeric and roman chromosomes are in the same vector
    chr[chr != "Un"] <- as.numeric( as.roman(chr[chr != "Un"]) ) 
        
	translate <- function(pos, startA, startB, endB, orientation){
	    # Translates from one coordinate system to a second
		if(!orientation == 'reverse'){
	    	pos2 = startB + (pos - startA)
	    	} else{
	      	pos2 = endB - (pos - startA)
	    	}
	    	return(pos2)
	  		}
	  
	scafTable <- read.csv(scafFile, header=TRUE, stringsAsFactors = FALSE)

	if(direction == 'old2new'){
		
		# repeat for every chr and pos: 
		z <- mapply(chr, pos, FUN = function(chr, pos){
		
		    # Pull out right scaffold
		    x <- scafTable[scafTable$OldChr == chr & scafTable$OldStart <= pos & scafTable$OldEnd >= pos, ]
	
		    # Make sure there's exactly 1 scaffold meeting criteria
			if( !(nrow(x) == 1) ){
		    	return (list(NA,NA))
		    	} else {
		    	# Calculate new coordinate
				newChr = x[1, 'NewChr']
				newPos = translate(pos, x[1,'OldStart'], x[1,'NewStart'], x[1,'NewEnd'], x[1, 'NewOrientation'])
				return(list(newChr = newChr, newPos = newPos))
		    	}
		    })
		return(t(z))
	    	
	  	} else if(direction == 'new2old'){

		# repeat for every chr and pos:
		z <- mapply(chr, pos, FUN = function(chr, pos){

		    # Pull out right scaffold
		    x <- scafTable[scafTable$NewChr == chr & scafTable$NewStart <= pos & scafTable$NewEnd >= pos, ]
	
		    # Make sure there's exactly 1 scaffold meeting criteria
		    if( !(nrow(x) == 1) ){
		    	return (list(oldChr = NA, oldPos = NA))
		    	} else {
				# Calculate new coordinate
				oldChr <- x[1,'OldChr']
				oldPos <- translate(pos, x[1,'NewStart'], x[1,'OldStart'], x[1,'OldEnd'], x[1,'NewOrientation'])
				return(list(oldChr = oldChr, oldPos  = oldPos))
				}
			})
		return(t(z))
			
		} else{
	    	print('Direction must be old2new or new2old')
	    	return(list(Chr = NA, Pos = NA))
	  		}
	}

g$genDistList <- function(alleleFreqByGroup, nMin, method = "nei"){
	# Genetic distance calculator between all pairs of groups or populations
	# "alleleFreqByGroup" is a list of tabulated allele frequenies, locus by locus
	# "nMin" is a vector with the minimum number of genotypes per group to make a calculation
	# 	nMinAlleles will be twice nMin
	# For a given locus, a genetic distance is NA if one of the members of the pair has too few alleles ( < nMinAlleles)

	# methods:
	#   "nei" (Nei 1972, see http://evolution.genetics.washington.edu/phylip/doc/gendist.html) and
	#	"pd" (proportional dissimilarity, ie 1 - proportional similarity")
	
	nMinAlleles <- 2 * nMin
	groupnames <- rownames(alleleFreqByGroup[[1]])
	# groupnames
	# [1] "paxl"       "paxb"       "pril"       "prib"       "qryl"       "qryb"       "ensl"      
	# [8] "ensb"       "marine-pac" "marine-atl" "marine-jap" "solitary"  

	namesAllpairs <- combn(groupnames, 2, FUN=function(x){paste(x, collapse=".")})
	# head(namesAllpairs)
	# [1] "paxl.paxb" "paxl.pril" "paxl.prib" "paxl.qryl" "paxl.qryb" "paxl.ensl"
	
	# Columns refer to rows of alleleFreqByGroup elements, NOT groupcodes
	# So the first column indicates 1 2, meaning row 1 vs row 2 of alleleFreqByGroup (paxl vs paxb)
	rowpairs <- combn(1:length(groupnames), 2)
	# rowpairs
	     # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17]
	# [1,]    1    1    1    1    1    1    1    1    1     1     1     2     2     2     2     2     2
	# [2,]    2    3    4    5    6    7    8    9   10    11    12     3     4     5     6     7     8
	# ...
	
	# FUNCTION: Proportional similarity
	ps <- function(x, i, j){
		# x is the full table of allele proportions, with populations in the rows and alleles in the columns,
		# e.g.,:
		              # 0  1  2  3
		  # paxl       12  0  0  0
		  # paxb       15  1  0  0
		  # pril       12  0  0  0
		  # prib        4  0  0  0
		# i and j indicate the two selected rows of the table of proportions to compute distance for
		P <- x[c(i, j), ]
		ps <- sum( apply(P, 2, min) )
		# print(ps)
		}
		
	# FUNCTION: Nei's (1972) genetic distance
	neidist1972 <- function(x, i, j){
		# x is the full table of allele proportions, with populations in the rows and alleles in the columns,
		# e.g.,:
		              # 0  1  2  3
		  # paxl       12  0  0  0
		  # paxb       15  1  0  0
		  # pril       12  0  0  0
		  # prib        4  0  0  0
		# i and j indicate the two selected rows of the table of proportions to compute distance for	.
		# Using Wixipedia notation:
		# Nei (1972) distance is -log( Jxy / (sqrt(Jx * Jy) ) where each J is summer over all loci
		Jxy <- sum(x[i, ] * x[j, ])
		Jx  <- sum(x[i, ] * x[i, ])
		Jy  <- sum(x[j, ] * x[j, ])
		return(c(Jxy = Jxy, Jx = Jx, Jy = Jy))
		}
	
	genDistList <- lapply(alleleFreqByGroup, function(x){
		# x <- alleleFreqByGroup[[50]]
		# x
		              # 0  1  2  3
		  # paxl       12  0  0  0
		  # paxb       15  1  0  0
		  # pril       12  0  0  0
		  # prib        4  0  0  0
		  # qryl        0  0  0  0
		  # qryb        0  0  0  0
		  # ensl        0  0  0  0
		  # ensb        4  0  0  0
		  # marine-pac  8  2  0  0
		  # marine-atl  2  0  0  0
		  # marine-jap  8  0  0  0
		  # solitary    8  0  0  0
		
		nAlleles <- rowSums(x)
		      # paxl       paxb       pril       prib       qryl       qryb       ensl       ensb marine-pac 
		        # 12         16         12          4          0          0          0          4         10 
			  # marine-atl marine-jap   solitary 
		         # 2          8          8 
		         
		# Which populations fail to meet the nMinAlleles criterion at the locus? 
		z <- nAlleles < nMinAlleles

		# Set failing populations to missing
		nAlleles[z] <- NA
		
		# Convert allele frequencues to proportions for those meeting nMinAlleles criterion
		alleleProp <- sweep(x, 1, nAlleles, FUN = "/")
		# alleleProp  
		                  # 0      1      2      3
		  # paxl       1.0000 0.0000 0.0000 0.0000
		  # paxb       0.9375 0.0625 0.0000 0.0000
		  # pril       1.0000 0.0000 0.0000 0.0000
		  # prib                                  
		  # qryl                                  
		  # qryb                                  
		  # ensl                                  
		  # ensb       1.0000 0.0000 0.0000 0.0000
		  # marine-pac 0.8000 0.2000 0.0000 0.0000
		  # marine-atl 1.0000 0.0000 0.0000 0.0000
		  # marine-jap 1.0000 0.0000 0.0000 0.0000
		  # solitary   1.0000 0.0000 0.0000 0.0000
		 
		if(method == "pd"){
			# Calculate proportional dissimilarity for all pairs meeting nMinALleles criterion
			pd <- apply(rowpairs, 2, function(z){
								pd <- 1 - ps(alleleProp, z[1], z[2])
								})
			pd[is.nan(pd)] <- NA
			names(pd) <- namesAllpairs
			
			# Confirm that the missing pairs are associated with those populations failing to meet nMinAlleles criterion
			# pd
			            # paxl.paxb             paxl.pril             paxl.prib             paxl.qryl 
			               # 0.0625                0.0000                    NA                    NA 
			# ...
			}
		
		if(method == "nei"){
			nei <- apply(rowpairs, 2, function(z){
				nei <- neidist1972(alleleProp, z[1], z[2])
				})
			colnames(nei) <- namesAllpairs
			# return(nei)
			}
			
		# How to put the results into a matrix (use in a blockstats version of this)
		# Use "njs" in ape package, which can handle some missing values
		# library(gdata) # needed for lowerTriangle()
		# d <- matrix(0, nrow = nrow(alleleProp), ncol = nrow(alleleProp))
		# rownames(d) <- rownames(alleleProp)
		# colnames(d) <- rownames(alleleProp)
		# lowerTriangle(d) <- pd
		# d <- t(d)
		# lowerTriangle(d) <- pd
		
		return(nei)
		})
	}

g$vcfTsTv <- function(REF, ALTlist, snpTypeList){
	# Calculates the raw transition-transversion ratio from the REF and list of ALT alleles, making
	# use of the snp types already provided in snpTypeList
	#
	n1 <- sapply(ALTlist, length) # length including NAs
	n2 <- sapply(snpTypeList, length)
	if( !all(n1 == n2) ) stop("ALTlist and snpTypeList have unequal element numbers")
	snp <- unlist(snpTypeList)
	ref <- rep(as.vector(REF), n1)
	alt <- unlist(ALTlist)

	cat("Table of variant types used (<*:DEL> and unused are NA)\n")
	print(table(snp, useNA="always"))
	
	# Transition-transversion ratios - true snp only (not indels)
	cat("\nTransition-transversion ratio - all true snp\n")
	snp1 <- snp[!is.na(snp) & snp == "snp"]
	ref1 <- ref[!is.na(snp) & snp == "snp"]
	alt1 <- alt[!is.na(snp) & snp == "snp"]
	
	ref1 <- substr(ref1, 1, 1) # keep the first letter of REF
	alt1 <- substr(alt1, 1, 1) # keep the first letter of ALT
	
	# c( length(ref1[ref1 != alt1]), length(alt1[ref1 != alt1]) )
	# [1] 216418 216418
	
	tstv <- g$tstv(ref1[ref1 != alt1], alt1[ref1 != alt1])
	return(tstv)
	}

g$tableAlleleFreqByGroup <- function(GT = geno(vcf)$GT, groupnames, groupcodes, nMaxAlt = 3, split = "/"){
	# Generates a table of allele frequencies at each locus by group
	# Group is specified by the integer groupcodes, eg 3 3 3 3 3 3 3 2 2 2 2 2 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 1
	# 
	groupnames <- factor(groupnames, levels = groupnames)
	tGT <- as.data.frame(t(GT), stringsAsFactors = FALSE) # annoying but apply won't return a list
	alleleFreqByGroup <- lapply(tGT, function(locus){
		# columns of tGT are loci, so apply function locus by locus
		# locus <- geno(vcf)$GT[6, ]
		z1 <- split(locus, groupcodes) # the genotypes for each group at the locus
		z2 <- lapply(z1, function(x){ 
			unlist(strsplit(x, split = split))
			}) 		# the alleles for each group at the locus
		z3 <- z2	# initialize
		for(i in 1:length(z3)) z3[[i]] <- rep( groupnames[i], length(z3[[i]]) ) # replace z3 with corresponding group name
		table( unlist(z3), factor(unlist(z2), levels=0:nMaxAlt) ) # factor so that all alleles are counted in each table
		})
	}

g$makeSnpTypeList <- function(REF = ref(vcf), ALTlist = altUsedList, AltDelSymbol = "<*:DEL>"){
	# Creates a list indicating what snp type is each ALT alleles: snp, ins, or del
	i <- nchar(REF) # length of the reference sequence
	i <- split(i, 1:length(i)) # split REF allele size into a list
	j <- lapply(ALTlist, function(x){
		j <- nchar(x)            	# length of ALT alleles in altList; 
									# is integer(0) if altUsedList[[i]] is character(0)
									# but is 2 if altUsedList[[i]] is NA
		j[x == AltDelSymbol] <- NA  # length of "<*:DEL>" set to NA
		j[is.na(x)] <- NA
		return(j)
		})
	snpTypeList <- mapply(i, j, FUN = function(i, j){ # i and j refer to nchar of ref and alt
		snptype <- rep(NA, length(j)) # initialize with NAs
		# snptype[i == 1 & j == 1] <- "snp" # first run
		# snptype[i < j  & i == 1] <- "ins" # first run
		# snptype[i > j  & j == 1] <- "del" # first run
		#
		# 2nd run below uses broader criteria to designate variant type
		# Note: j is based on ALTlist, so snpTypeList is also
		snptype[i == j] <- "snp" # 2nd run
		snptype[i <  j] <- "ins" # 2nd run
		snptype[i >  j] <- "del" # 2nd run
		snptype[length(j) == 0] <- NA # snp type is NA if length of altUsedList is 0, ie no snp used
		return(snptype)
		})
	names(snpTypeList) <- names(ALTlist)
	return(snpTypeList)
	}

g$makeAltUsedList <- function(GT, ALT = alt(vcf)){
	# Returns a list of the ALT alleles actually used in the genotypes in the sample
	# The list is based on alt(vcf) but unused ALT alleles are set to NA 
	# 	- NB they are NOT DELETED to preserve indices
	# GT is geno(vcf)$GT or a subset with rows as loci and columns as individuals
	#
	if(nrow(GT) != length(ALT)) stop("Number of rows of GT must equal length of ALT")
	z <- as.list(ALT) # mapply takes forever if use alt(vcf) instead of as.list(alt(vcf))
	whichAltAllelesUsed <- g$whichAltAllelesUsed(GT)
	altUsedList <- mapply(whichAltAllelesUsed, z, FUN = function(x, y){
		# x <- whichAltAllelesUsed[[7]]; y <- alt(vcf)[[7]]
		ialt <- 1:length(y)
		y[!(ialt %in% x)] <- NA
		return(y)
		})
	return(altUsedList)
	}

g$whichAltAllelesUsed <- function(GT, split = "/"){
	# returns a list indicating which alleles are actually used in the genotypes 
	# 	in the sample by their index (1, 2, ...
	# GT is geno(vcf)$GT or a subset with rows as loci and columns as individuals
	whichAltAllelesUsed <- apply(GT, 1, function(x){ 		# is a list
			# x <- GT[1, ]
			x1 <- strsplit(x[!is.na(x)], split = split)
			x1 <- sort( as.integer(unique(unlist(x1))) )
			whichAltAllelesUsed <- x1[x1 > 0] # drops the REF allele
			})
	}

g$makeChrvecMasked <-  function(chrname, windowsmaskername){
	chr <- scan(paste(chrname, ".fa", sep = ""), what=character())  # get genome
	chr <- chr[-1] # drops the line with ">chr??"
	chr <- paste(chr, collapse="") # combined all the lines into a single word
	chrvec <- strsplit(chr, split="")[[1]] # break apart genome into individual bases
	
	#table(chrvec) # lower case refers to masked bases from the initial assembly (? maybe not in gasAcu1.fa?)
	#      A       C       G       N       T 
	#3133419 2580881 2578205  275208 3149774
	
	# Using repeatmasker output Felicity gave me on Oct 31, 2012
	chrmaskinfo <- read.table(windowsmaskername, stringsAsFactors = FALSE) 
	names(chrmaskinfo) <- c("chrno","from","to")
	
	#head(chrmaskinfo)
	   # chrno from   to
	# 1 chrXXI    0  423
	# 2 chrXXI  435 3277
	# 3 chrXXI 3285 6311
	# 4 chrXXI 6319 7888
	# 5 chrXXI 7976 8026
	# 6 chrXXI 8082 8874
	
	# ---
	# Get the positions of all the masked bases
	# NOTE: windowmasker start position is 0-based, so add 1 (end position is 1-based!)
	# Note: chrmaskinfo[ ,-1] used because otherwise each vector x is converted to all character and "seq" generates error
	
	z <- apply(chrmaskinfo[ ,-1], 1, function(x){seq(from=x[1] + 1, to=x[2])}) 
	m <- unlist(z) 
	
	chrvec[m] <- "M" # convert all these bases to "M" to indicate masked regions 
	
	cat("Summary of masked chromosome:")
	print(table(chrvec)) # Notice that "N" is not present anymore, so these were also masked in windowsmaskerSdust.chrXXI.txt
	# chrvec
	  # A       C       G       M       T 
	# 2464777 2130503 2128953 2521828 2471426 
	
	# Save to .rdd file in vcfdir
	chrno <- gsub("^chr", "", chrname)
	chrmaskfile <- paste("chrvec.", chrno, ".masked.rdd", sep = "") # chrvec.XXI.masked.rdd
	cat("\nSaving masked chrvec object to", chrmaskfile, "\n")
	save(chrvec, file = chrmaskfile)
	#return(chrvec)
	}        
  
g$downloadSra <- function(SraSample, fishName){ # eg g$downloadSra("SRS639967", "KODK")
	# R commands for Westgrid to download SRA files from NCBI
	# based on "Using the SRAdb Package to Query the Sequence Read Archive", Zhu & Davis

	# SraSample <- "SRS639967"; fishName <- "KODK"
	# source("http://bioconductor.org/biocLite.R"); biocLite("SRAdb")
	library(SRAdb)
	sqlfile <- 'SRAmetadb.sqlite'
	if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()
	
	z <- list.files(pattern = "^SRR.")
	if(length(z) > 0) stop("Method assumes there are no files beginning with 'SRR' in working directory")
	
	# Open a connection
	sra_con <- dbConnect(SQLite(), sqlfile) # opens a connection
	
	# Obtain all the experiment accessions, run accessions, etc, associated with the sample identifier
	sra_info <- sraConvert( SraSample, sra_con = sra_con )
	#sra_info
	      # sample submission     study experiment        run
	# 1  SRS639967  SRA164437 SRP042012  SRX621222 SRR1449383
	# 2  SRS639967  SRA164437 SRP042012  SRX608534 SRR1437855
	# 3  SRS639967  SRA164437 SRP042012  SRX608534 SRR1462291
	# 4  SRS639967  SRA164437 SRP042012  SRX608534 SRR1462120
	# 5  SRS639967  SRA164437 SRP042012  SRX621222 SRR1449714
	# 6  SRS639967  SRA164437 SRP042012  SRX608534 SRR1449659
	# 7  SRS639967  SRA164437 SRP042012  SRX608534 SRR1462754
	# 8  SRS639967  SRA164437 SRP042012  SRX621222 SRR1449492
	# 9  SRS639967  SRA164437 SRP042012  SRX608534 SRR1462185
	# 10 SRS639967  SRA164437 SRP042012  SRX608534 SRR1422274
	# 11 SRS639967  SRA164437 SRP042012  SRX608534 SRR1462825
	
	sra_list <- apply(sra_info, 2, unique)
	# $sample
	# [1] "SRS639967"
	# $submission
	# [1] "SRA164437"
	# $study
	# [1] "SRP042012"
	# $experiment
	# [1] "SRX621222" "SRX608534"
	# $run
	 # [1] "SRR1449383" "SRR1437855" "SRR1462291" "SRR1462120" "SRR1449714"
	 # [6] "SRR1449659" "SRR1462754" "SRR1449492" "SRR1462185" "SRR1422274"
	# [11] "SRR1462825"
	
	# Lists all the ftp addresses of the sra files
	# listSRAfile( SraSample, sra_con, fileType = 'sra' )
	
	# Gets the ftp addresses of the fastq files corresponding to the run
	# getFASTQinfo( sra_list$run, srcType = 'ftp' )
	# Downloads them!
	getSRAfile( sra_list$run, sra_con, fileType = 'fastq' )
	
	# Rename the files to match population for easy recognition using file.rename(from, to)
	z <- list.files(pattern = "^SRR.")
	for(i in 1:length(z)){
		file.rename(z[i], paste(fishName, z[i], sep="_"))
		}
	cat("Files downloaded and renamed:\n")
	list.files(pattern=paste("^", fishName,".", sep = ""))
	}

g$recode<-function (..., ret = c("numeric", "factor"), none = if (ret == "numeric") 0 else "none", na){
	# stolen from the Hmisc library
    ret <- match.arg(ret)
    w <- list(...)
    if (!is.logical(w[[1]]) && length(w) == 3) {
        z <- w[[3]][match(w[[1]], w[[2]])]
        if (!missing(none)) 
            z[if (is.numeric(none)) 
                is.na(z)
            else z == ""] <- none
        return(z)
    }
    nam <- names(w)
    if (missing(ret)) 
        ret <- if (all.is.numeric(nam)) 
            "numeric"
        else "factor"
    result <- rep(none, length(w[[1]]))
    for (i in 1:length(w)) result[w[[i]]] <- if (ret == "numeric") 
        numnam[i]
    else nam[i]
    if (ret == "factor") 
        result <- as.factor(result)
    if (!missing(na)) 
        result[is.na(na)] <- NA
    result
}

g$geno.diff <- function(GT1, GT2){
	# calculates the divergence between two genotype vectors of equal length (e.g., columns of GT matrix of vcf file)
	# Genotypes of individuals are in VCF "GT" format:
		# NA   "0/0"   "1/1"      NA   "0/1"   "0/1"   "0/0"      NA   "1/1"   "1/1"   "0/0"   "0/0"
	#
	# Tests that this is done correctly:
	# g$geno.diff("0/0", "0/1") # [1] 0.5
	# g$geno.diff("0/0", "1/1") # [1] 1
	# g$geno.diff("0/0", "0/0") # [1] 0
	# g$geno.diff("0/0", "0/1") # [1] 0.5
	# g$geno.diff("0/1", "0/0") # [1] 0.5
	# g$geno.diff("0/0", "1/1") # [1] 1
	# g$geno.diff("0/1", "2/3") # [1] 1
	# g$geno.diff("0/1", "1/2") # [1] 0.5
	# g$geno.diff(c("0/1","0/0"), c(NA,"0/0")) 	# [1] NA  0
	#
	# GT1 <- c(NA, "0/0", "0/0", "0/0", "0/1", "0/1", "0/1", "1/1", "1/1", "0/1", "0/1", "0/1")
	# GT2 <- c(NA, "0/0", "0/1", "1/1", "0/0", "0/1",    NA, "0/0", "0/1", "1/2", "2/3", "1/4")
	x <- strsplit(GT1, split = "/") 
	y <- strsplit(GT2, split = "/") 
	diff <- mapply(x, y, FUN = function(xi, yi){
		# xi <- x[[2]]; yi <- y[[2]]
		diff <- ceiling((4-sum(is.element(xi,yi), is.element(yi,xi)))/2)/2
		if(any(is.na(c(xi,yi)))) diff <- NA
		diff
		})
	diff
}

g$psd <- function(genotypes, groups){
	# Calculate the number of genotypic differences between all pairs of individuals, one pair at a time.
	# Pairwise sequence difference "psd" is 0, 0.5, and 1 according to alleles shared between individuals at the locus.
	# "genotype" is a matrix of genotypes of individuals in VCF "GT" format:
		# NA   "0/0"   "1/1"      NA   "0/1"   "0/1"   "0/0"      NA   "1/1"   "1/1"   "0/0"   "0/0"
	# "groups" indicates corresponding group membership of individuals.	
	# genotypes <- rbind( c(NA, "0/0", "0/0", "0/0", "0/1", "0/1", "0/1", "1/1", "1/1", "0/1", "0/1", "0/1"),
	#                     c(NA, "0/0", "0/1", "1/1", "0/0", "0/1",    NA, "0/0", "0/1", "1/2", "2/3", "1/4") )
	# groups <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)
	# genotypes <- x1$GT

	pop <- as.character(groups)  # population names

	gtpairs <- combn(pop, 2, FUN=function(x){paste(x, collapse=",")})
	# gtpairs
	 # [1] "1,1" "1,1" "1,1" "1,1" "1,1" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,1" "1,1" "1,1" "1,1" "1,2" "1,2" "1,2" "1,2"
	# [20] "1,2" "1,2" "1,1" "1,1" "1,1" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,1" "1,1" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2"
	# [39] "1,1" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "1,2" "2,2" "2,2" "2,2" "2,2" "2,2" "2,2"
	# [58] "2,2" "2,2" "2,2" "2,2" "2,2" "2,2" "2,2" "2,2" "2,2"
	
	colpairs <- combn(1:ncol(genotypes), 2)
	colpairs
	     # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17]
	# [1,]    1    1    1    1    1    1    1    1    1     1     1     2     2     2     2     2     2
	# [2,]    2    3    4    5    6    7    8    9   10    11    12     3     4     5     6     7     8
	     # [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
	# [1,]     2     2     2     2     3     3     3     3     3     3     3     3     3     4     4     4
	# [2,]     9    10    11    12     4     5     6     7     8     9    10    11    12     5     6     7
	     # [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49]
	# ...
	
	psd <- apply(colpairs, 2, function(z){
			psd <- g$geno.diff(genotypes[ ,z[1] ], genotypes[ ,z[2] ])
			})
	
	psd <- as.data.frame(psd, stringsAsCharacters = FALSE)
	names(psd) <- gtpairs
	psd
	}


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


g$slidewin <- function(blockstats, method = "FST", nsteps.per.window, windowNmin){

	# Combines results of g$blockstats across nsteps.per.window of bins at a time in sliding window.

	# "ibase" in blockstats is the start of the sliding window of 5 blocks (e.g.).
	# rollapply adds padding to the left of the data range near the start of the data set and
	# 	to the right of the data range near the end; padding is determined by "fill" option (NA is the default).

	# In all cases, windowNmin is the smallest number of good bases needed in each full sliding window to yield result, 
	# 	otherwise window's value is set to NA. 

	# FST: returns the Weir-Cockeram VARiances obtained by summing variance component across all bases in the window
	#		and then calculating Fst as the usual ratio of the variance components. Monomorphic
	#		bases in a window do not influence this Fst, because they have no variance. 
	#		However, the sliding windows for VARiance components include the 0's where there are data.
	#
	# CSS: returns Felicity's CSS score, following the formula in Jones et al 2012 Supplement Eq. 1.
	# 		I assume that she used metric multidimensional scaling on the pairwise sequences divergences between
	#			pairs of fish within each sliding window. Calculation of CSS uses Euclidean distances between 
	#			pairs of fish on just the first 2 axes from this MDS analysis.
	#		The formula for CSS given in the Jones et al Suppl is
	#				CSS <- ( sij/(m*n) ) - (m+n)*( sii/((m-1)*(m^2)) + sjj/((n-1)*(n^2)) )
	# 		but the correct formula is
	#				CSS <- ( sij/(m*n) ) - (1/(m+n))*( sii/((m-1)/2) + sjj/((n-1)/2) ) 
	# 		where sij/(m*n) is the mean of all pairwise distances between individuals of differenc species, lim and ben
	# 			(1/(m+n))*( sii/((m-1)/2) ) here works out to be the mean pairwise distance between individuals within ben
	# 			(1/(m+n))*( sjj/((n-1)/2) ) here works out to be the mean pairwise distance between individuals within lim
	#					
	#		To handle missing values, Felicity et al calculated 21x21 pairwise distances
	# 		between the 21 populations of fish in each sliding window. Any pairwise distance based on fewer
	#		that 100 good bases (out of 2500 bases in the window, ie 4 %) was tossed and replaced with the
	#		matrix average of the remaining pairwise comparisons. If more than 50% of the pairwise distances
	#		had to be tossed, the whole result for the window was set to NA.		
	# 		My value for the minimum number of validated positions is windowNmin
	# 		Note: I have already replaced NA's with the mean psd (separately for each locus.
	# 		So there should be no missing psd's within an informative locus by this point

	#		I handled missing values differently from Felicity et al. 
	#		1. For a given pair of fish I calculated psd (perc seq div) at every nucleotide (snp and non-snp).
	#		2. If a psd value was NA at a nucleotide, I used an average psd from the remaining pairs of fish
	#			(according to the rule specified by psdMissingAction -- see blockstats)
	#		3. For the given pair of fish, I took the average psd across all nucleotides in the block of 500 (eg),
	#			(averaging over both the real and substituted "missing" psd values. 
	# 		4. Good bases already use a criterion for a minimum number of genotypes per group, for variant and invariant
	#			sites. This ensures that most pairwise divergence values are present for every good base.
	#		5. I also employ a minimum number of good bases per window to generate a result, otherwise the result
	#			from a window is set to NA.

  	# Calculate sliding window values, weighting by n and N.
	library(zoo)
	
	stepsize <- blockstats$ibase[2] - blockstats$ibase[1]
	k <- stepsize # this is the size of the block
	wink <- nsteps.per.window # number of blocks per sliding window
	
	if( casefold(method) == "fst" ){
		FST <- rollapply(blockstats, width = wink, by.column = FALSE, fill = NA, 
				FUN = function(x){    # remember: argument "x" is a window, so has multiple rows
					#print(x[,1:9])
					VARa <- x[, "VARa"] 
					VARb <- x[, "VARb"]
					VARw <- x[, "VARw"]
					
					nbases <- sum(x[, "nSnp"] + x[, "nInvariants"], na.rm = TRUE)
					if(nbases < windowNmin){
						totVARa <- NA
						totVARb <- NA
						totVARw <- NA
						FST <- NA
						}
					else{
						totVARa <- sum(VARa, na.rm=TRUE)
						totVARb <- sum(VARb, na.rm=TRUE)
						totVARw <- sum(VARw, na.rm=TRUE)
						FST <- totVARa/(totVARa + totVARb + totVARw)
						FST[is.nan(FST)] <- NA
						}
					#print(c(nbases, totVARa, totVARb, totVARw, FST))
					c(nbases = nbases, totVARa = totVARa, totVARb = totVARb, totVARw = totVARw, fst = FST)
					})
		#print(FST)
		results <- cbind.data.frame(ibase = blockstats$ibase, FST)		
		}


	if( casefold(method) == "css" ){
				
		# indicate identity of pairs, ie "1,1", "1,2", "2,2"
		# gtpairs <- combn(pop, 2, FUN = function(x){paste(x, collapse=",")}) # 1,1 1,2 2,2 indicates groups

		is.psd <- grep("psd.", names(blockstats))
		gtpairs <- names(unlist(blockstats[1,]))[is.psd] # only way I found to stop R adding suffix to duplicate names
		gtpairs <- gsub( "psd.", "", gtpairs )
		npop <- ceiling(sqrt(2*length(gtpairs)))
		
		# This next bit assumes just two populations types
		#poptype <- unique(pop)
		poptype <- sort( unique(unlist(strsplit(gtpairs, split = ","))) )
		if(length(poptype) != 2) stop("CSS only works for two population types")
		
		ii <- gtpairs == paste(poptype[1], poptype[1], sep = ",") # cases of "1,1"
		jj <- gtpairs == paste(poptype[2], poptype[2], sep = ",") # cases of "2,2"
		ij <- gtpairs == paste(poptype[1], poptype[2], sep = ",") # cases of "1,2"
		ji <- gtpairs == paste(poptype[2], poptype[1], sep = ",") # cases of "2,1"
		m <- ceiling(sqrt(2 * length(gtpairs[ii])))
		n <- ceiling(sqrt(2 * length(gtpairs[jj])))
		
		# check
		if( length(c(gtpairs[ij], gtpairs[ji])) != m*n ) stop("Number of pairs of individuals doesn't match up")
			
		# A <- grep(poptype[1],colnames(euclid)) # e.g., benthics
		# B <- grep(poptype[2],colnames(euclid)) # e.g., limnetics
		# m <- length(A)
		# n <- length(B)

		library(gdata) # needed for lowerTriangle()
		CSS <- rollapply(blockstats, width = wink, by.column = FALSE, fill = NA, 
				FUN = function(x){
				# FUN = function(x,  # remember: argument "x" is a window, so has multiple rows
						# is.psd = is.psd,
						# gtpairs = gtpairs,
						# npop = npop,
						# poptype = poptype,
						# ii = ii,
						# jj = jj,
						# m = m,
						# n = n
						# ){   
				# x <- blockstats[106:110, ]
					
					# xpairs <- x[, gtpairs]
					xpairs <- x[ ,is.psd ]
					
					sumdist <- apply(xpairs, 2, sum, na.rm = TRUE) # sum of summed pairwise dist over blocks in window
					names(sumdist) <- gtpairs
					
					# total number of bases in window, assume that CSStrueSnps = TRUE
					nbases <- sum(x[, "nSnp"] + x[, "nInvariants"], na.rm = TRUE) 
										
					if(nbases < windowNmin) CSS <- NA
					else{
						# Here I'm just calculating the sum of psd for all pairs, for now ignoring
						# snpCountByFish, ie whether some individuals may have too few genotypes

						# Put sumdist into matrix
						# first get [i,j] of the pairs for the matrix d
												
						d <- matrix(0, nrow = npop, ncol = npop)
						# Test that this will work. 
						# groups <- c(2,2,2,2,2,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1) # from "paxl" "paxb" analysis
						# lowerTriangle(d) <- gtpairs; d <- t(d); lowerTriangle(d) <- gtpairs
						# rownames(d) <- groups
						# colnames(d) <- groups
						# d 
						# This shows that lowerTriangle() puts the gtpairs into the matrix down the columns, ie
						#	it fills up the first column, then continues on the second, and so on. So the gtpair "2,1"
						#   refers to "column,row". This is why we use lowerTriangle() to fill the matrix, so that
						#	the results match gtpairs.

						lowerTriangle(d) <- sumdist
						d <- t(d)
						lowerTriangle(d) <- sumdist
						# rownames(d) <- groups
						# colnames(d) <- groups
						
						# Following the formula in Jones et al 2012 Supplement Eq. 1
						# Assume that she used metric multidimensional scaling
						# The following command generates a lot of warnings of the sort:
						#  "1: In cmdscale(d, k = 2, add = TRUE) : only 0 of the first 2 eigenvalues are > 0"
						# 	so I have suppresed them in case they halt execution.
						suppressWarnings(
							z <- cmdscale(d, k = 2, add = TRUE)$points  # it crashed if add=FALSE used
						)
												
						# calculate pairwise distance between observations on the cmdscale axes
						euclid <- as.matrix(dist(z, method = "euclidean")) 	# matrix of distances;
						euclid.lower <- euclid[lower.tri(euclid)] 			# pairwise distances, order matches gtpairs
						
						# sij <- sum(euclid[A,B])    # is the same as when use euclid[B,A] instead
						# sii <- sum(euclid[A,A][lower.tri(euclid[A,A])])
						# sjj <- sum(euclid[B,B][lower.tri(euclid[B,B])])
						
						sii <- sum(euclid.lower[ii])
						sjj <- sum(euclid.lower[jj])
						sij <- sum( sum(euclid.lower[ij]), sum(euclid.lower[ji]), na.rm=TRUE )
						
						CSS <- ( sij/(m*n) ) - (1/(m+n))*( sii/((m-1)/2) + sjj/((n-1)/2) ) 
						}
						
					return(c(nbases = nbases, CSS = CSS))
					})
		results <- cbind.data.frame(ibase = blockstats$ibase, CSS)			
		detach(2) # to avoid conflicts among libraries
		}
		
	results
}

g$plot.reads <- function(x1, from = 1, to = nrow(x1), scalefactor = 1, minqual = 30, keepunmappedmate=TRUE){
	# x1 is a data frame made from scanBam object
	# If they don't already exist, the following LOGICAL variables are created:
	# "unmapped" "unmappedmate", "bothforward","bothreverse"
	# Positive insert size: 
	# 	draws a line connecting the start of a read (top) to the end of a read (bottom)
	# 	draws only reads whose position (pos) at the START of read is between "from" and "to"
	# Negative insert size: 
	#	draws a line connecting the end of a read (bottom) backwards to the start of a read (top)
	# 	draws only reads whose position (pos) at the END of read is between "from" and "to"
	# scalefactor is the amount that the x-axis in plot is extended beyond "for" and "to", in multiples of to-from

	if(!is.data.frame(x1)) stop("Error: object must be a data frame\n")
	interval <- to - from
	if(interval < 0) stop("Error: to must be greater than from\n")
	library(bitops)
	if(!exists("unmapped")) unmapped <- bitAnd(x1$flag, 0x0004) == 0x0004
	if(!exists("unmappedmate")) unmappedmate <- bitAnd(x1$flag, 0x0008) == 0x0008
	if(!exists("bothforward")) bothforward <- 
			bitAnd(x1$flag, 0x0010) == 0x0000 & bitAnd(x1$flag, 0x0020) == 0x0000 & !unmapped & !unmappedmate
	if(!exists("bothreverse")) bothreverse <- 
			bitAnd(x1$flag, 0x0010) == 0x0010 & bitAnd(x1$flag, 0x0020) == 0x0020 & !unmapped & !unmappedmate

	# regular isizes, positive insert sizes
	if(keepunmappedmate){
		i <- !bothforward & !bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize > 0} 
	else{
		i <- !unmappedmate & !bothforward & !bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize > 0}
	dat <- x1[i & x1$pos>from & x1$pos<to,]
	plot(c(from-interval*scalefactor,to+interval*scalefactor),c(0,1), type="n", 
		xlab="positive inserts",ylab="")
	segments(dat$pos, rep(1,nrow(dat)), dat$pos+71+dat$isize, rep(0,nrow(dat)), lwd=0.2)

	# regular isizes, negative insert sizes
	if(keepunmappedmate){
		j <- !bothforward & !bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize < 0} 
	else{
		j <- !unmappedmate & !bothforward & !bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize < 0}
	datj <- x1[j & x1$pos>from & x1$pos<to,]
	plot(c(from-interval*scalefactor,to+interval*scalefactor),c(0,1), type="n", 
		xlab="negative inserts",ylab="")
	segments(datj$pos, rep(0,nrow(datj)), datj$pos+datj$isize, rep(1,nrow(datj)), lwd=0.2)

	# both reads point same way, positive insert sizes
	i <- bothforward & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize > 0
	j <- bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize > 0
	dati <- x1[i & x1$pos>from & x1$pos<to,]
	datj <- x1[j & x1$pos>from & x1$pos<to,]
	plot(c(from-interval*scalefactor,to+interval*scalefactor),c(0,1), type="n", 
		xlab="positive inserts, both reads point same direction",ylab="")
	segments(dati$pos, rep(1,nrow(datj)), dati$pos+71+dati$isize, rep(0,nrow(dati)), lwd=0.5, col="red")
	segments(datj$pos, rep(1,nrow(datj)), datj$pos+71+datj$isize, rep(0,nrow(datj)), lwd=0.5, col="blue")

	# both reads point same way, negative insert sizes
	i <- bothforward & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize < 0
	j <- bothreverse & !is.na(x1$mapq) & x1$mapq >= minqual & x1$isize < 0
	dati <- x1[i & x1$pos>from & x1$pos<to,]
	datj <- x1[j & x1$pos>from & x1$pos<to,]
	plot(c(from-interval*scalefactor,to+interval*scalefactor),c(0,1), type="n", 
		xlab="negative inserts, both reads point same direction",ylab="")
	segments(dati$pos, rep(0,nrow(datj)), dati$pos+dati$isize, rep(1,nrow(dati)), lwd=0.5, col="red")
	segments(datj$pos, rep(0,nrow(datj)), datj$pos+datj$isize, rep(1,nrow(datj)), lwd=0.5, col="blue")

	invisible()
	}


g$peel <- function(x,i){
	# Peels off the i'th elements from a comma-separated character-based set of numbers
	# 	e.g., "peel(PV4,1)" peels the first element "0.34" of PV4 "0.34,1.5e-20,1,1"
	# Repeats for every entry if x is a vector
	# Returns a data frame if i is a vector
	# 	e.g., "peel(PV4,c(1,3))" returns the first and third elements of PV4
	z1 <- strsplit(x, split=",")
	z2 <- sapply(z1,function(x,i){x[i]},i[1])
	if(length(i) > 1){
		z3 <- data.frame(z2)
		for(j in 2:length(i)){
				z3[j] <- sapply(z1,function(x,i){x[i]},i[j])
			}
		names(z3) <- paste("peel",1:length(i), sep="")
		return(z3)
		}
	else return(z2)
	}


g$flag.readable <- function(i){
	# i must be a vector of integers, which are samtools flags
	# i is converted to a human-readable version of the flag
	# For example, i-69 is converted to:
	# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    #   01   00   01   00   00   00   01   00
	#    p    P    u    U    r    R    1    2  samtools flag codes
	# interpretation of the 8 binary numbers (1=yes, 0=no)
	# 1	read paired?
	# 2	read mapped in a proper pair?
	# 3   query sequence unmapped?
	# 4	mate is unmapped?
	# 5 	query is reverse strand?
	# 6   mate is reverse strand?
	# 7	read is first in pair?
	# 8	read is second in pair?
	# Note also that it is possible to query integer i using the bitops library
	# E.g., to query if mate is unmapped, use "bitAnd(x1$flag, 0x0004) == 0x0004"
	# (not implemented)
	# library(bitops)
	# 0x0001	p	the read is paired in sequencing
	# 0x0002	P	the read is mapped in a proper pair
	# 0x0004	u	the query sequence itself is unmapped
	# 0x0008	U	the mate is unmapped
	# 0x0010	r	strand of the query (1 for reverse)
	# 0x0020	R	strand of the mate
	# 0x0040	1	the read is the first read in a pair
	# 0x0080	2	the read is the second read in a pair
	z <- sapply(i, function(i){
		z <- as.logical(intToBits(i)[1:8])
		z1 <- c("p","P","u","U","r","R","1","2")[z]
		z2 <- paste(z1,collapse="")
		})
	return(z)
	}


g$gtf2thirdpositions <- function(gtffilename, refseqname, directory = "", check=TRUE, dropNonTriplets=TRUE,
	collapse=TRUE){
	# Analyzes a single chromosome
	# gtffilename is name of "group.gtf" file (in quotations)
	# refseqname is name of "chr.fa" file (in quotations)
	# if collapse=FALSE, returns a list by transcript, otherwise the list is collapsed into a single
	#	vector with unique bases
	# If collapse=FALSE, result names are a paste of transcript name : gene name
	# use as "z <- g$gtf2thirdpositions("groupXXI.gtf", "chrXXI.fa", check=TRUE, dropNonTriplets=TRUE, collapse=TRUE)"
		
	library(GenomicFeatures)
	library(Biostrings)
	
	# setwd("/Users/schluter/Documents/Research/genomics/reference genomes")
	# gtffilename <- "groupXXI.gtf"
	# refseqname <- "chrXXI.fa"
	
	# make the transcript data base
	gtf <- makeTranscriptDbFromGFF(file= gtffilename, format="gtf") # has class "TranscriptDb"
	
	# Pull out gene names corresponding to each transcript
	z <- select(gtf, keys = keys(gtf), keytype="GENEID",
		cols = c("CDSSTART", "CDSEND", "EXONID", "EXONSTRAND", "EXONSTART","EXONEND", "TXNAME"))
	# Warning message:
	# In .generateExtraRows(tab, keys, jointype) :
  	# 'select' resulted in 1:many mapping between keys and return rows
	transcriptGeneName <- tapply(z$GENEID, z$TXNAME, unique)


	# Extract GRanges from transcript data base
	grange <- transcripts(gtf)
	
	#head(grange)
      # seqnames           ranges strand |     tx_id            tx_name
   	     # <Rle>        <IRanges>  <Rle> | <integer>        <character>
  	# [1] groupXXI [240732, 252486]      + |         1 ENSGACT00000002162
  	# [2] groupXXI [280690, 284889]      + |         2 ENSGACT00000002166
  	# [3] groupXXI [308111, 318180]      + |         3 ENSGACT00000002174
  	# [4] groupXXI [308262, 317945]      + |         4 ENSGACT00000002176
  	# [5] groupXXI [614128, 678875]      + |         5 ENSGACT00000002201
  	# [6] groupXXI [698931, 733565]      + |         6 ENSGACT00000002210
	#length(grange)
	# [1] 599
	if(check) cat("Number of transcripts = ", length(grange), "\n")

	# Pull out all the coding regions by transcript
	# NOTE: expect different transcripts of the same gene to overlap?
	# The method uses an internal ordering to ensure that transcripts with the same name are unique, just in case,
	# but this makes it hard to identify the gene and transcript in the end. 
	# Can do "use.names=TRUE" if no duplicate names - this worked for groupXXI
	cdslist <- cdsBy(gtf, by = "tx", use.names=TRUE)

	#head(cdslist)
	# $ENSGACT00000002162   # these are transcript names
    	  # seqnames           ranges strand |    cds_id    cds_name exon_rank
    	     # <Rle>        <IRanges>  <Rle> | <integer> <character> <integer>
  	# [1] groupXXI [240732, 240840]      + |         1        <NA>         1
  	# [2] groupXXI [241074, 241090]      + |         2        <NA>         2
  	# [3] groupXXI [242286, 242437]      + |         3        <NA>         3
  	# [4] groupXXI [242593, 242612]      + |         4        <NA>         4
  	# [5] groupXXI [243731, 243789]      + |         5        <NA>         5
  	# [6] groupXXI [243941, 244036]      + |         6        <NA>         6
  	# [7] groupXXI [244147, 244257]      + |         7        <NA>         7
  	# [8] groupXXI [244463, 244599]      + |         8        <NA>         8
  	# [9] groupXXI [251813, 251912]      + |         9        <NA>         9
	#
	# $ENSGACT00000002166 
    	  # seqnames           ranges strand | cds_id cds_name exon_rank
  	# [1] groupXXI [281806, 281887]      + |     10     <NA>         2
  	# [2] groupXXI [281980, 282131]      + |     11     <NA>         3
  	# [3] groupXXI [282801, 283077]      + |     12     <NA>         4
  	# [4] groupXXI [284468, 284502]      + |     13     <NA>         5
	# ...

	# Number of transcripts represented
	# This is not 599 beause some transcripts are missing data?
	#length(cdslist) 
	# [1] 586
	#setdiff( grange$tx_name,  names(cdslist)) # these are the missing transcripts
 	# [1] "ENSGACT00000029571" "ENSGACT00000029368" "ENSGACT00000002915" "ENSGACT00000029474" "ENSGACT00000029188"
 	# [6] "ENSGACT00000029179" "ENSGACT00000029436" "ENSGACT00000028331" "ENSGACT00000029353" "ENSGACT00000029344"
	# [11] "ENSGACT00000028203" "ENSGACT00000029513" "ENSGACT00000028077"
	if(check) cat("Number of transcripts with coding region groups = ", length(cdslist), "\n")

	# Check that all cda's belonging to each transcript have the same strand ("+" or "-")
	z1 <- sapply(cdslist, function(x){
		length(unique(strand(x)))
		})
	if(check){
		cat("Frequency table indicating the number of unique strand directions per CDS group","\n",
			"(should be all 1)","\n")
		print(table(z1))
		}
		# z1
  		# 1
		# 586 

	# Check that the CDS are divisible by 3
	z1 <- sapply(cdslist, function(x){
		sum(width(x))/3
		})
		#head(z1)
		# ENSGACT00000002162 ENSGACT00000002166 ENSGACT00000002174 ENSGACT00000002176 ENSGACT00000002201 ENSGACT00000002210 
	               # 267                182                552                510                573               1039 
	z2 <- z1[z1!=round(z1,0)]
	if(check){
		cat("These transcripts do not have an integer number of triplets","\n")
		if(dropNonTriplets) cat("they have been removed","\n")
		print(z2)
		}	
		# ENSGACT00000003716 ENSGACT00000006614 
        		  # 441.3333           174.6667 
	
	#cdslist['ENSGACT00000003716'] # nothing unusual about this transcript
	
	if(dropNonTriplets & length(z2) > 0) cdslist <- cdslist[-which(is.element(names(cdslist),names(z2)))]

	# Check whether exons are always arranged along chr in increasing order (or decreasing in the case of strand = "-")
	if(check){
		 z1 <- sapply(cdslist, function(z){
			if(as.logical(strand(z)[1]=="+")) z1 <- all(start(z) == sort(start(z))) else 
			z1 <- all(start(z) == sort(start(z), decreasing=TRUE)) 
			})
	cat("Test of whether exons occur in sequential order along chromosome","\n")
	print(table(z1))
		# z1
		# TRUE 
 		# 584
		}

	# Identify bases in first, second, or third base positions of transcripts
	chromseq <- readDNAStringSet(refseqname, "fasta")[[1]]

	# Find the base pairs corresponding to third positions
	thirdpositions <- lapply(cdslist, function(z){
	# z <- cdslist[['ENSGACT00000002162']] # strand is "+"
	# z <- cdslist[['ENSGACT00000003737']] # strand is "-"
		z1 <- start(z) # start position for each range
		z2 <- end(z)   # end position of each range
		if(as.logical(strand(z)[1]=="+")){ 
			z3 <- mapply(FUN = seq, z1, z2) # seq(from=z1, to=z2) to get all nucleotides for each range
			z3 <- unlist(z3) 				# collapse nucleotide positions from all ranges to a single vector
			# all(z3==sort(z3))
			z4 <- z3[seq(3, length(z3), by = 3)] # start at the 3rd base and increment to end by 3 bases at a time
			} else{
			z3 <- mapply(FUN = seq, z2, z1) # seq(from=z2, to=z1) as above but reversing the order
			z3 <- unlist(z3) 				# collapse nucleotide positions from all ranges to a single vector
			# all(z3==sort(z3, decreasing=TRUE))
			# remember: z3 starts high and counts down, 
			# 	i.e., if head(z3) is (7102135 7102134 7102133 7102132 7102131 7102130), 
			# 	then the first instance of a third base position is 7102133
			# so this is the correct sequence of third-position bases:
			z4 <- z3[seq(3, length(z3), by = 3)] # starts at the 3rd from the end and increments down by 3
			}
		return(z4)
		})
	# thirdpositions['ENSGACT00000002162']
	# thirdpositions['ENSGACT00000003737']
	if(collapse) thirdpositions <- unique(unlist(thirdpositions))
	else{
		# Change names of the list to include the gene name, so that this information is retained
		z1 <- transcriptGeneName[match(names(thirdpositions), names(transcriptGeneName))]
		names(thirdpositions) <- paste(names(thirdpositions), z1, sep=":")
		}
	return(thirdpositions)
	}
	


g$tstv <- function(REF,ALT){
	# calculates the raw transition to transversion ratio R
	if(!(length(REF)==length(ALT))) stop("Vectors not of same length")
	ref <- rep(NA,length(REF))
	alt <- rep(NA,length(ALT))
	ref[is.element(REF,c("A","G"))] <- "pur"
	ref[is.element(REF,c("C","T"))] <- "pyr"
	alt[is.element(ALT,c("A","G"))] <- "pur"
	alt[is.element(ALT,c("C","T"))] <- "pyr"
	z <- table(ref,alt)
	R <- list(R = (z[1,1]+z[2,2]) / (z[1,2]+z[2,1]), tstv = z, tot = sum(z))
	return(R)
}
	
g$write.bedGraph <- function(filename, chrinfo = "chrXXI:1-11,717,487", from, to, stat,
	graphname = "DolphsGraph", col = "41,155,241", altcol = "100,100,100", ylinemark = 5,
	description = "chrXXI sliding window")
	{
	# bedGraph allows scientific notation in the stat variable, but the base lmits must be integers.
	z <- cbind.data.frame(chr = rep(gsub("[:][0-9,-]*", "",chrinfo), length(stat)), from = from, to = to, stat = stat)
	z <- z[!is.na(z$stat),]
	lim <- round(range(stat, na.rm=TRUE),5)
	lim <- paste(lim[1],lim[2],sep=":")
	write(file=filename, paste("browser position",chrinfo, sep=" ")) 
	write(file=filename, append=TRUE, 
		paste("track type=bedGraph name=", graphname, ' description="', description, '"',
	" color=",col, " altColor=",altcol, " viewLimits=",lim, " yLineOnOff=on yLineMark=", ylinemark, 
	" windowingFunction=mean autoScale=off visibility=full priority=20", sep=""))
	# "format(z)" suppresses scientific notation in the integer bases
	write.table(format(z), file=filename, append=TRUE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
	invisible()	
}

g$wc.revised <- function(ndat){
	# modified "wc" from hierfstat version 0.04-14 to speed it up
	library(hierfstat)
	# library(parallel)
	# ndat <- cbind(xpop,geno[,1:10000]) # test

    diploid <- TRUE # no longer an argument
    pol <- 0        # no longer an argument
    pop <- ndat[, 1]
    ni <- length(pop)
    dat <- ndat
    loc.names <- names(dat)[-1]

	# ---
	# The following code to count number of genotyped individuals per population is too slow, 
	# replace with subsequent block of code

    # system.time({
    # n <- t(ind.count(dat)) 
	# }) # system time 4-7 sec with 10000 markers analyzed on 2 groups size 6 each
    
    # system.time({
    popf <- factor(pop)
    z <- split(pop, popf)
    z1 <- lapply(z, function(z){
		z1 <- apply(dat[ which(pop %in% z), -1], 2, function(x){sum(!is.na(x))})
		})
	n <- do.call("cbind", z1)
	# }) # system time is about 1 second
	# (n[n[,1]==0 | n[,2]==0,])[1:100,] # many lines with one of the columns having a 0
	        # 1 2
	# V1      2 0
	# V22     1 0
	# V25     1 0
	# V42     0 4
	# V43     0 6
# ---
    
    nt <- apply(n, 1, sum, na.rm = TRUE)
    z <- table(nt)
    if(!is.na(z["0"])) warning("some markers have 0 individuals, which will affect length of output")
    untyped.loc <- which(nt == 0) # integer(0)
    typed.loc <- which(nt != 0)
    if (length(untyped.loc) > 0) {
        dat <- dat[, -(untyped.loc + 1)]
        n <- t(ind.count(dat))
        nt <- apply(n, 1, sum, na.rm = TRUE)
    	}
    
	# -------
	# The following code to count total number of alleles is too slow, 
	# replace with subsequent block of code

    # system.time({
    # alploc <- nb.alleles(cbind(rep(1, ni), dat[, -1]))
    # }) # takes about 10 seconds
    	
    # system.time({
    # Break up genotypes into alleles and keep the result -- slows things here but speeds up a later step
    # Missing alleles are set to "9" and then back to NA so I can split them
    allAlleles <- lapply(dat[ ,-1], function(x){
    		x[is.na(x)] <- 99
    		z <- unlist(strsplit(as.character(x), split = ""))
    	z[z == "9"] <- NA
    	return(z)
    		})
    	nalleles <- lapply(allAlleles, function(x){length(unique(na.omit(x)))})
    	alploc <- matrix(nalleles, ncol = 1, dimnames = list(loc.names, NULL))
    # }) # takes about 2 seconds
	# -------
    
    np <- dim(n)[2]
    npl <- apply(n, 1, tempfun <- function(x) sum(!is.na(x)))
    nl <- dim(n)[1]
    
    # -----
    # system.time({
    # p <- pop.freq(dat, diploid=diploid)
    # pb <- pop.freq(cbind(rep(1, length(pop)), dat[, -1]), diploid)
    # }) # took 22 seconds
    
    # system.time({
    popf <- factor(pop) 							# [1] 1 1 1 1 1 1 2 2 2 2 2 2
    popalleles <- rep(popf, rep(2, length(popf)))	# [1] 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2
    alleleFreqTable <- lapply(allAlleles, function(x){ 
#    alleleFreqTable <- mclapply(allAlleles, FUN=function(x){ # parallel lapply - does it work on hermes?
    		# x <- allAlleles[[175]]
		z1 <- table(x, popalleles)
		})	
    p <- lapply(alleleFreqTable, function(x){
		z1 <- sweep(x, 2, colSums(x), FUN = "/")
		})
	pb <- lapply(alleleFreqTable, function(x){
		z <- rowSums(x)
		z2 <- matrix(z/sum(z), ncol = 1, dimnames = list(names(z), "1"))
		})
    # }) # took 8 seconds
    
    # -----
    
    n <- matrix(unlist(n), ncol = np)
    nal <- n[rep(1:nl, alploc), ]
    nc <- (nt - apply(n^2, 1, sum, na.rm = TRUE)/nt)/(npl - 1)
    ntal <- rep(nt, alploc)
    ncal <- rep(nc, alploc)
    p <- matrix(unlist(lapply(p, t)), ncol = np, byrow = TRUE)
    pb <- matrix(unlist(pb), ncol = 1)
	dum <- getal.b(dat[, -1]) # slow
	all.loc <- apply(dum, 2, tempfun1 <- function(y) as.numeric(dimnames(table(y))[[1]])) # slow
	hetpl <- apply(dum, 2, fun <- function(z) { # slow
		lapply(as.numeric(dimnames(table(z))[[1]]), who.is.het <- function(y) apply(z == 
			y, 1, ind.is.het <- function(x) xor(x[1], x[2])))
			})
	mho <- lapply(hetpl, tempfun2 <- function(x) matrix(unlist(lapply(x, 
			tempfun3 <- function(y) tapply(y, pop, sum, na.rm = TRUE))), 
			ncol = np))
	mho <- matrix(unlist(mho), ncol = np, byrow = TRUE)
	mhom <- (2 * nal * p - mho)/2
    SSG <- apply(nal * p - mhom, 1, sum, na.rm = TRUE)
    dum <- nal * (p - 2 * p^2) + mhom
    SSi <- apply(dum, 1, sum, na.rm = TRUE)
    dum1 <- nal * (sweep(p, 1, pb))^2
    SSP <- 2 * apply(dum1, 1, sum, na.rm = TRUE)
    ntalb <- rep(npl, alploc)
    MSG <- SSG/ntal
    MSP <- SSP/(ntalb - 1)
    MSI <- SSi/(ntal - ntalb)
    sigw <- MSG
    sigb <- 0.5 * (MSI - MSG)
    siga <- 1/2/ncal * (MSP - MSI)
    Fxy <- function(x) x[1]/sum(x, na.rm = TRUE)
    FST.pal <- apply(cbind(siga, sigb, sigw), 1, Fxy)
    FIS.pal <- apply(cbind(sigb, sigw), 1, Fxy)
    loc <- rep(1:nl, alploc)
    lsiga <- tapply(siga, loc, sum, na.rm = TRUE)
    lsigb <- tapply(sigb, loc, sum, na.rm = TRUE)
    lsigw <- tapply(sigw, loc, sum, na.rm = TRUE)
    sigloc <- cbind(lsiga, lsigb, lsigw)
    if (length(untyped.loc) > 0) {
        x <- order(c(typed.loc, untyped.loc))
        dum1 <- matrix(numeric(3 * length(untyped.loc)), ncol = 3)
        dum <- rbind(sigloc, dum1, deparse.level = 0)[x, ]
        sigloc <- dum
    	}
    sigloc <- data.frame(sigloc)
    names(sigloc) = c("lsiga", "lsigb", "lsigw")
    rownames(sigloc) <- loc.names
    lFST <- apply(cbind(lsiga, lsigb, lsigw), 1, Fxy)
    lFIS <- apply(cbind(lsigb, lsigw), 1, Fxy)
    tsiga <- sum(siga, na.rm = TRUE)
    tsigb <- sum(sigb, na.rm = TRUE)
    tsigw <- sum(sigw, na.rm = TRUE)
    tFST <- Fxy(c(tsiga, tsigb, tsigw))
    tFIS <- Fxy(c(tsigb, tsigw))
    # res <- list(call = cl, sigma = cbind(loc, siga, sigb, sigw), 
        # sigma.loc = sigloc, per.al = list(FST = FST.pal, FIS = FIS.pal), 
        # per.loc = list(FST = lFST, FIS = lFIS), FST = tFST, FIS = tFIS)
    res <- list(sigma = cbind(loc, siga, sigb, sigw), 
        sigma.loc = sigloc, per.al = list(FST = FST.pal, FIS = FIS.pal), 
        per.loc = list(FST = lFST, FIS = lFIS), FST = tFST, FIS = tFIS)
    #}
    class(res) <- "wc"
    res
}

