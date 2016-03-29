# Combined the snp data from Jones et al Current Biol paper: 
#	jones et al 2012 curr biol - global and species-pair patterns of snp divergence - supplement.xls
# with the SNP data at NCBI from the Kinglsey lab:
#	http://www.ncbi.nlm.nih.gov/SNP/snp_batchSearch.cgi?org=69293&type=SNP&from=&to=&condh=starts+with&handle=&condb=starts+with&lbid=&Search+Batch=Search+Batch

# Page 1 (Table S2) of jones et al 2012 curr biol - global and species-pair patterns of snp divergence - supplement.xls is in .csv file
# It contains info on the 1159 SNP used in the Jones et al study.

# The file "dbsnp_snpbatch_concatenated.txt" is the concatenated, downloaded data from NCBI of the Kingley SNPs deposited
# 	( files are at http://www.ncbi.nlm.nih.gov/SNP/snp_batchSearch.cgi?org=69293&type=SNP&from=&to=&condh=starts+with&handle=&condb=starts+with&lbid=&Search+Batch=Search+Batch )

# rsnps package didn't give much information:
# install.packages("rsnps", dependencies = TRUE)
# library(rsnps)
# NCBI_snp_query('rs119103371')
        # Query Chromosome      Marker Class Gene Alleles Major Minor MAF BP
# 1 rs119103371       <NA> rs119103371   snp <NA>     A,G     A     G  NA NA
# Warning message:
# In NCBI_snp_query("rs119103371") :
  # No chromosomal information for rs119103371; may be unmapped
# -----

setwd("/Volumes/schluter/Dolph_data/genomics/dpSNP")
setwd("/Users/schluter/Documents/Research/genomics/dpSNP")

# Kingsley all snps. Contains the allele data and rs codes
# chr position is contained in the submitter.snp.id 
x <- read.table("dbsnp_snpbatch_concatenated.txt", header = TRUE, stringsAsFactors = FALSE)
head(x)
           # ss          submitter.snp.id allele samplesize          rs ss2rsOrien  chr contig.acc
# 1 ss120258411 GasAcu1.0_LG01_00,913,033    A/G         20 rs119103371          + N.D.       N.D.
# 2 ss120258412 GasAcu1.0_LG01_03,494,580    A/G         20 rs119103372          + N.D.       N.D.
# 3 ss120258413 GasAcu1.0_LG01_07,955,458    A/C         20 rs119103373          + N.D.       N.D.
# 4 ss120258414 GasAcu1.0_LG01_09,345,491    A/T         20 rs119103374          + N.D.       N.D.
# 5 ss120258415 GasAcu1.0_LG01_11,963,492    A/C         20 rs119103375          + N.D.       N.D.
# 6 ss120258416 GasAcu1.0_LG01_12,038,660    A/C         20 rs119103376          + N.D.       N.D.

# drop some of the unneeded variables
x$ss2rsOrien <- NULL
x$chr <- NULL
x$contig.acc <- NULL

# Jones et al snps. Contains chr positions and ss codes (rs codes incomplete)
y <- read.csv("jones et al 2012 curr biol - global and species-pair patterns of snp divergence - supplement p1.csv", stringsAsFactors = FALSE)
head(y)
        # snpname              local.snp.id     ss.id     rs.id snp.group
# 1 chrI:11963492 GasAcu1.0_LG01_11,963,492 120258415 119103375         2
# 2 chrI:12038660 GasAcu1.0_LG01_12,038,660 120258416 119103376         2
# 3 chrI:12199899 GasAcu1.0_LG01_12,199,899 418641996        na         2
# 4 chrI:12231274 GasAcu1.0_LG01_12,231,274 418641997        na         2
# 5  chrI:1245655 GasAcu1.0_LG01_01,245,655 418641981        na         1
# 6  chrI:1320011 GasAcu1.0_LG01_01,320,011 418641982        na         1

# Transfer over the SNP names in the Jones et al spreadsheet (a subset of the Kingsley SNPs)
x$ss <- gsub("ss","", x$ss)
x$snpname.jones <- matchcol(x$ss, y$ss.id, y$snpname)
head(x)
         # ss          submitter.snp.id allele samplesize          rs snpname.jones
# 1 120258411 GasAcu1.0_LG01_00,913,033    A/G         20 rs119103371   chrI:913033
# 2 120258412 GasAcu1.0_LG01_03,494,580    A/G         20 rs119103372  chrI:3494580
# 3 120258413 GasAcu1.0_LG01_07,955,458    A/C         20 rs119103373  chrI:7955458
# 4 120258414 GasAcu1.0_LG01_09,345,491    A/T         20 rs119103374  chrI:9345491
# 5 120258415 GasAcu1.0_LG01_11,963,492    A/C         20 rs119103375 chrI:11963492
# 6 120258416 GasAcu1.0_LG01_12,038,660    A/C         20 rs119103376 chrI:12038660

# Create a chr variable from the snp id name
x$chr <- sapply( strsplit(x$submitter.snp.id, split = "_"), function(x){x[2]})
# x$chr <- gsub("Pitx1", "VIIpitx1", x$chr)  # use this if including the OLD pitx1 BAC
x$chr <- gsub("Pitx1", "VIIpitx1new", x$chr) # use this if including Felicity's NEW pitx1 BAC
x$chr <- gsub("LG", "", x$chr)

# Drop SNPs on the Y (all have Y in name)
z <- grepl("Y", x$chr)
x$chr[z]
# [1] "KitY" "KitY" "IdhY" "IdhY"
x <- x[!z, ]

# Convert numbers to roman
z <- grep("[0-9]{2}", x$chr)
x$chr[z] <- as.character(as.roman(x$chr[z]))
x$chr <- paste("chr", x$chr, sep = "")
unique(x$chr)
 # [1] "chrI"           "chrII"          "chrIII"         "chrIV"          "chrV"          
 # [6] "chrVI"          "chrVII"         "chrVIII"        "chrIX"          "chrX"          
# [11] "chrXI"          "chrXII"         "chrXIII"        "chrXIV"         "chrXV"         
# [16] "chrXVI"         "chrXVII"        "chrXVIII"       "chrXIX"         "chrXX"         
# [21] "chrXXI"         "chrUn"          "chrVIIpitx1new" "chrM"          

# Here we obtain the POS from the snp id names, but there's a 0-base vs 1-base issue that the
# results in the actual POS being one base pair more for nuclear chromosomes (not mtDNA!)
# Instead, use BLAT to identify actual POS and which allele is REF
# x$POS <- gsub( ",", "", sapply( strsplit(x$submitter.snp.id, split = "_"), function(x){x[3]}) )
# x$POS <- as.integer(x$POS)
# x$submitter.snp.id <- NULL
         # ss allele samplesize          rs snpname.jones  chr      POS
# 1 120258411    A/G         20 rs119103371   chrI:913033 chrI   913033
# 2 120258412    A/G         20 rs119103372  chrI:3494580 chrI  3494580
# 3 120258413    A/C         20 rs119103373  chrI:7955458 chrI  7955458
# 4 120258414    A/T         20 rs119103374  chrI:9345491 chrI  9345491
# 5 120258415    A/C         20 rs119103375 chrI:11963492 chrI 11963492
# 6 120258416    A/C         20 rs119103376 chrI:12038660 chrI 12038660

# -----------------
# Get flanking sequence for every snp

# Read SNP primer sequences from "rs_chNotOn (flanking sequence).fas"
# downloaded from ftp://ftp.ncbi.nih.gov/snp/organisms/stickleback_69293/rs_fasta/

y <- scan(file = "rs_chNotOn (flanking sequence).fas", what = character(), 
	comment.char = "#", sep = "\n", strip.white = TRUE)
y <- matrix(y, ncol = 4, byrow=TRUE)
z <- strsplit(y[,1], split = "[| ]+")
y[,1] <- sapply(z, function(z){z[3]})
colnames(y) <- c("rs","leftFlank","snp","rightFlank")
# nrow(y)
# [1] 1644
# nrow(x)
# [1] 1640

x$leftFlank <- matchcol(x$rs, y[,"rs"], y[,"leftFlank"])
x$snp <- matchcol(x$rs, y[,"rs"], y[,"snp"])
x$rightFlank <- matchcol(x$rs, y[,"rs"], y[,"rightFlank"])
head(x)
         # ss allele samplesize          rs snpname.jones  chr      POS
# 1 120258411    A/G         20 rs119103371   chrI:913033 chrI   913033
# 2 120258412    A/G         20 rs119103372  chrI:3494580 chrI  3494580
# 3 120258413    A/C         20 rs119103373  chrI:7955458 chrI  7955458
# 4 120258414    A/T         20 rs119103374  chrI:9345491 chrI  9345491
# 5 120258415    A/C         20 rs119103375 chrI:11963492 chrI 11963492
# 6 120258416    A/C         20 rs119103376 chrI:12038660 chrI 12038660
                                                          # leftFlank snp
# 1 TAGTGTTCAT ACATCTGCAG TATTATGGAC CGGAGGACAT TACTGACGAA AACATGCAGA   R
# 2 TGCATCCGAT CTGCACATAG AACGACTCCT CTTACCGAAG GCTCCTGAAC TCCCCCACAC   R
# 3 CTGTCTTTTG ACTGCAAGTG GACATTTGAT GCCTGAAACA TTTTCCTTTA ACTATGCATA   M
# 4 AAAAGAGAAT TGTGGAGCTG AGTATCCTGA TCAATTTTAT TACTCCTAGG ATGGATCTAC   W
# 5 AGCCTTGGGA ACGCAAGTAG CCTGCCAAAC TTAACATTAA CAACATCTTG AAATACTTTC   M
# 6 CTTATTTTTT TTATGTACTA AGATCTTTAT TTAATTGGGA TTTTTTTATT GAAGCCGCAG   M
                                                         # rightFlank
# 1 TTAGTTTTCC TCCCTTACTA CGCAGTGTTT ACGGTATTAA TACAGTAAAT GTTGGGAAAA
# 2 GGTGTCTACC GTTGTCTAAT ACTGGTGGTA TACGCAGAAA AAGCACAAGT TTCACAATAT
# 3 ATCGGTCTGT GGTGCTTATG AAATTGGAGT GATTTAAATT GTTTTAGTAC AACAAGCCTA
# 4 CACAATTGGG TATAAGATAC CTTATCATGG TGCACATCAC ACAAGTAATG TATACGTAAA
# 5 GAGCGTATTA ATCCATCTTC TAAATTGGAA CTTGTACTTT TTAAATTGTT GTCAAGTTTT
# 6 AGATGAGATT GTTTGATTGC ACAGGGAGCG TCCTCCGTAC TCTTCATGGT CCATCTGTAT

# Save interim results
# write.csv(x, "kingsley_snps.csv", row.names = FALSE) # if using the original pitx1 BAC
write.csv(x, "kingsley_snpspitx1new.csv", row.names = FALSE) # if using Felicity's NEW pitx1 BAC

# ---------------
# Use flanking sequence to get actual POS and REF allele for every snp using BLAT

# Generate fasta files containing flanking sequence for every snp
for(i in unique(x$chr)){
	# i <- "chrXXI"
	z <- x[x$chr == i , c("rs", "leftFlank", "snp", "rightFlank")]
	z$rs <- paste(">", z$rs)
	z <- as.matrix( z )
	z <- matrix(t(z), ncol = 1)
	outfile = paste("temp", i, "fa", sep = ".")
	write(z, file = outfile)
	}

# Copy the "temp" files to hermes (westgrid) and run "blat" on each chromosome
# Run R on hermes:
module load R/3.1.2
R
filenames <- list.files(pattern=glob2rx("temp.chr*.fa"))
filenames <- filenames[grep("chr", filenames)]
for(i in filenames){
	# i <- "temp.chrVIIpitx1new.fa"
	chrname <- sapply( strsplit(i, split = "[.]"), function(z){z[2]})
	database <- paste(chrname, "fa", sep = ".")
	outfile <- paste(chrname, "blat.out", sep = ".")
	# This BLAT command produced the same default output as in Ensembl, according to the Ensembl help
	cmd <- paste("/global/software/blat-3.4/bin/blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0", 
					database, i, outfile, "&", sep = " ")
	system(cmd)
	}
# Then copy *blat.out files back to computer

# Extract snp information from blat.out files
genomeDirName <- "../reference genomes/"
for(i in unique(x$chr)){
	# i <- "chrVIIpitx1new"
	blatfile <- paste(i, "blat.out", sep = ".")
	y <- read.table(blatfile, header = FALSE, stringsAsFactors = FALSE, skip = 5)
	y <- y[ , c(1,4,9,10,16,17)]
	names(y) <- c("match", "Ns", "strand", "rs", "start", "end")
	
	# Drop anything that is not a match of 120 (out of 121 bases blatted)
	y <- y[y$match == 120, ]
	
	head(y)
	  # match Ns strand          rs    start      end
	# 1   120  1      + rs119103522  3772198  3772319
	# 2   120  1      - rs119103523  8268391  8268512
	# 3   120  1      + rs119103524 10007823 10007944
	# 4   120  1      - rs119103525 10511222 10511343
	# 5   120  1      + rs119103526 11060149 11060270
	# 6   120  1      + rs119103527 11414323 11414444
	
	# Checked the results against Ensembl BLAT results, which has nice visuals
	# - rs119103522: Ensembl got "groupXXI:3772199-3772319" forward strand, middle snp "C" at  pos 3772259
	#	(which is 3772198 + 61 = 3772259, where 60 is the length of the left flanking region)
	# - rs119103523: Ensembl got "groupXXI:8268392-8268512", reverse strand, middle snp "A" at pos 8268452
	#	(which is 8268391 + 61 = 8268452). Complement of "A" on reverse strand is "T" on the forward strand.

	y$POS <- y$start + 61
	
	# Grab the reference genome sequence
	fastaName <- paste(genomeDirName, i, ".fa", sep = "")
	chr <- scan(fastaName, what=character()) 
	z <- grepl(">", chr) # is there a ">chr??" at the top
	chr <- chr[!z]  # drops it
	chr <- paste(chr, collapse="") # combined all the lines into a single word
	chrvec <- strsplit(chr, split="")[[1]] # break apart genome into individual bases
	rm(chr)
	# length(chrvec)
	# [1] 11717487
	
	# Obtain the corresponding REF base
	y$REF <- chrvec[y$POS]
	
	head(y)
	  # match Ns strand          rs    start      end      POS REF
	# 1   120  1      + rs119103522  3772198  3772319  3772259   C
	# 2   120  1      - rs119103523  8268391  8268512  8268452   T
	# 3   120  1      + rs119103524 10007823 10007944 10007884   G
	# 4   120  1      - rs119103525 10511222 10511343 10511283   G
	# 5   120  1      + rs119103526 11060149 11060270 11060210   C
	# 6   120  1      + rs119103527 11414323 11414444 11414384   G

	# Put results back into "x"
	x$POS[x$chr==i] 	<- matchcol(x$rs[x$chr==i], y$rs, y$POS)
	x$REF[x$chr==i] 	<- matchcol(x$rs[x$chr==i], y$rs, y$REF)
	x$strand[x$chr==i] 	<- matchcol(x$rs[x$chr==i], y$rs, y$strand)
	
	# x[x$chr==i, c(3,5,7,11,12,13)]
	     # allele          rs    chr      POS REF strand
	# 152     C/G rs119103522 chrXXI  3772259   C      +
	# 153     A/G rs119103523 chrXXI  8268452   T      -
	# 154     A/G rs119103524 chrXXI 10007884   G      +
	# 155     C/G rs119103525 chrXXI 10511283   G      -
	# 156     C/G rs119103526 chrXXI 11060210   C      +
	# 157     A/G rs119103527 chrXXI 11414384   G      +
	# 553     A/C rs119103923 chrXXI  5737466   A      +
	# 554     A/G rs119103924 chrXXI  5791520   A      +
	# 555     A/G rs119103925 chrXXI  5793104   T      -
	# 556     A/G rs119103926 chrXXI  6037993   A      +
	# 557     A/T rs119103927 chrXXI  7002179   T      +
	# 558     A/G rs119103928 chrXXI  7544042   G      +
	# 559     A/C rs119103929 chrXXI 10054379   T      -
	# 560     A/C rs119103930 chrXXI 10236130   C      +
	# 561     A/G rs119103931 chrXXI 10969153   G      +
	# 562     A/G rs119103932 chrXXI 11179444   A      +
	# 653     G/T rs119104113 chrXXI  7904440   T      +
	# 689     A/G rs119104149 chrXXI 10156752   A      +
	# 703     A/G rs119104163 chrXXI  1893295   G      +
	# 715     A/G rs119104175 chrXXI   774193   A      +
	# 720     A/C rs119104180 chrXXI  5716516   C      +
	# 767     G/T rs193648464 chrXXI   242381   G      +
	# 1415    A/T rs193649112 chrXXI  1146722   T      +
	# 1416    G/T rs193649113 chrXXI  1895966   G      +
	# 1417    A/G rs193649114 chrXXI  3082228   A      +
	# 1418    A/T rs193649115 chrXXI  6356755   A      +
	# 1419    A/C rs193649116 chrXXI  7466745   C      +
	# 1420    C/T rs193649117 chrXXI  9373718   C      +
	# 1421    A/C rs193649118 chrXXI  9820535   A      +
	# 1422    G/T rs193649119 chrXXI 11532914   T      +
	# 1519    C/T rs193649216 chrXXI   772528   T      +
	# 1520    G/T rs193649217 chrXXI  5648650   T      +
	# 1521    A/C rs193649218 chrXXI  5717549   A      +
		
	}

# Any missing?
table(x$REF, useNA = "always")
   # A    C    G    T <NA> 
 # 379  438  416  406    1 
 
# Drop the one missing case (on pitx1)
x <- x[!is.na(x$REF), ]

# When strand is "-" the SNP bases are the complement of the REF (and ALT) bases
# To get the ALT allele, need to convert, get the alt base, and convert back

# First, convert REF to the complementary base if strand is reverse (-)
snpref <- x$REF
snpref[x$strand == "-"] <- recode( snpref[x$strand == "-"], c("A","C","G","T"), c("T","G","C","A") )
 
# Next, figure out corresponding alt allele
z <- strsplit(x$allele, split = "/")
snpalt <- mapply(z, snpref, FUN = function(z,x){
	ref <- z[!is.element(z,x)]
	})

table(snpalt, useNA="always")
# snpalt
   # A    C    G    T <NA> 
 # 547  334  454  304    0 # results if using new pitx1 BAC (same as for original BAC, I think)
 
# Finally, convert the alt allele to the ALT allele on the complementary strand when strand is "-"
x$ALT <- snpalt
x$ALT[x$strand == "-"] <- recode( x$ALT[x$strand == "-"], c("A","C","G","T"), c("T","G","C","A") )

# Save results
names(x)[names(x) == "POSnew"] <- "POS"
# write.csv(x, "kingsley_snps.csv", row.names = FALSE)  # Use this is using original pitx1 BAC
write.csv(x, "kingsley_snpspitx1new.csv", row.names = FALSE) # Use this is using Felicity's new pitx1 BAC

# -----------------
# Generate the .vcf files for known snps
# Assuming that we are running on hermes (directory is /g01/home/schluter/)

# According to GATK, this is the minimum cvf file:
	
	##fileformat=VCFv4.1
	##contig=<ID=chrXXI,length=11717487>
	##reference=file:///g01/home/schluter/chrXXI.fa
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
	chrXXI	1	.	A	T	.	.	.

chrnames <- unique(x$chr)
genomeDirName <- "../reference genomes/"

for(i in chrnames){
	# i <- "chrXXI"
	fastaName <- paste(i, ".fa", sep = "")
	vcfname <- paste("knownSnp", i, "vcf", sep = ".")
	chr <- scan(paste(genomeDirName, fastaName, sep = ""), what=character()) 
	z <- grepl(">", chr) # is there a ">chr??" at the top
	chr <- chr[!z]  # drops it
	chr <- paste(chr, collapse="") # combined all the lines into a single word
	chrvec <- strsplit(chr, split="")[[1]] # break apart genome into individual bases
	chrlength <- length(chrvec)
	rm(chr)
	rm(chrvec)

	header <- c("##fileformat=VCFv4.1",
				paste("##contig=<ID=", i, ",length=", chrlength, ">", sep = ""),
				paste("##reference=file:///g01/home/schluter/", fastaName, sep = ""),
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	write(header, file = vcfname)
	
	y <- x[x$chr == i, c("chr", "POS", "REF", "ALT")]
	y <- y[order(y$POS),]
	z <- apply(y, 1, function(y){
		paste(y[1], y[2], ".", y[3], y[4], ".", ".", ".", sep = "\t")
		})
	write(z, file = vcfname, append = TRUE)
	
	}


