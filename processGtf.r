# Make stickleback transcript data base from gtf file

# First need to convert Ensembl gtf file to gasAcu1 coordinates (scaffolds converted to chrUn, etc)
# Then make the data base

# ==================
# 1) 
# Convert Ensembl gtf file to gasAcu1 coordinates 
# gtf file obtained at ftp://ftp.ensembl.org/pub/release-88/gtf/gasterosteus_aculeatus/

# I used Method B: Make "gasAcu1.ensembl.88.cleaned.gtf.gz" from "Gasterosteus_aculeatus.BROADS1.88.gtf.gz"
# Code is below.

# ==================
# 2)
# Create TranscriptDb objects from modified gtf files
# based on "/genomes/r files/makeSticklebackTranscriptDb_object_from_GTF.v88.r"

library(GenomicFeatures)

	# Using Biomart args. Has class "TranscriptDb"
gtfDb <- makeTxDbFromGFF(file="gasAcu1.ensembl.88.cleaned.gtf.gz", 
			format="gtf", organism = "Gasterosteus aculeatus",
			taxonomyId = 69293) 
	# Import genomic features from the file as a GRanges object ... OK
	# Prepare the 'metadata' data frame ... OK
	# Make the TxDb object ... OK
	# Warning messages:
	# 1: Named parameters not used in query: internal_chrom_id, chrom, length, is_circular 
	# 2: Named parameters not used in query: internal_id, name, type, chrom, strand, start, end 
	# 3: Named parameters not used in query: internal_id, name, chrom, strand, start, end 
	# 4: Named parameters not used in query: internal_id, name, chrom, strand, start, end 
	# 5: Named parameters not used in query: internal_tx_id, exon_rank, internal_exon_id, internal_cds_id 
	# 6: Named parameters not used in query: gene_id, internal_tx_id 

# The above warnings means that most parameters in gtf file are not included in the gtfDb object, such as gene names, etc.

# Maybe just make a dataframe with all this missing information for later determination of gene id?
# Can construct using (Not run)
	x <- read.table("gasAcu1.ensembl.88.cleaned.gtf.gz", sep = "\t", stringsAsFactors = FALSE, quote = "")
	names(x) <- c("seqid","source","feature","start","end","score","strand","frame","attributes")
	x1 <- strsplit(x$attributes, split = ";")


saveDb(gtfDb, file="gasAcuEnsemblGtf88.sqlite")
# gtfDb <- loadDb("gasAcuEnsemblGtf88.sqlite")
	# TxDb object:
	# Db type: TxDb
	# Supporting package: GenomicFeatures
	# Data source: gasAcu1.ensembl.88.cleaned.gtf.gz
	# Organism: Gasterosteus aculeatus
	# Taxonomy ID: 69293
	# miRBase build ID: NA
	# Genome: NA
	# transcript_nrow: 29245
	# exon_nrow: 245522
	# cds_nrow: 234262
	# Db created by: GenomicFeatures package from Bioconductor
	# Creation time: 2017-04-21 12:22:29 -0700 (Fri, 21 Apr 2017)
	# GenomicFeatures version at creation time: 1.26.4
	# RSQLite version at creation time: 1.1-2
	# DBSCHEMAVERSION: 1.1

# You are done. But here are some useful accessor functions
# See help for "select-methods {GenomicFeatures}" to interfqace with TxDb objects

gtfDb
metadata(gtfDb)
keytypes(gtfDb)
keys(gtfDb)
tx_chrom(gtfDb)
isActiveSeq(gtfDb) # which chromosomes are currently active
	# If you then wanted to only set Chromosome 1 to be active you could do it like this: (not run)
isActiveSeq(gtfDb)[seqlevels(gasaculEnsembl)] <- FALSE
isActiveSeq(gtfDb) <- c("chrXXI" = TRUE)
seqlevels(gtfDb) # can use: seqlevels(x) <- value
	 # [1] "chrI"     "chrII"    "chrIII"   "chrIV"    "chrV"     "chrVI"    "chrVII"   "chrVIII"  "chrIX"    "chrX"     "chrXI"   
	# [12] "chrXII"   "chrXIII"  "chrXIV"   "chrXV"    "chrXVI"   "chrXVII"  "chrXVIII" "chrXIX"   "chrXX"    "chrXXI"   "chrM"    
	# [23] "chrUn"   
seqinfo(gtfDb) # can use: seqinfo(x) <- value

transcripts(gtfDb)
transcriptsBy(...)

# Method A - edit the database after importing
----------------------------------------------
gtfDb <- makeTxDbFromGFF(file="Gasterosteus_aculeatus.BROADS1.88.gtf.gz", 
			format="gtf", organism = "Gasterosteus aculeatus",
			taxonomyId = 69293) # Using Biomart args. Has class "TranscriptDb"

metadata(gtfDb)
saveDb(gtfDb, file="gasAcuEnsemblGtf88.sqlite")
# gtfDb <- loadDb("gasAcuEnsemblGtf88.sqlite")

# To edit, dump to a list (taken from help page for TxDb objects)
	# txdb_dump <- as.list(txdb)
	# txdb_dump
	# txdb1 <- do.call(makeTxDb, txdb_dump)
	# stopifnot(identical(as.list(txdb1), txdb_dump))
txdb_dump <- as.list(gtfDb)
names(txdb_dump)
# [1] "transcripts" "splicings"   "genes"       "chrominfo"  

* Got to here* 
* Now edit the list items (shown below) and then go back and recreate the data base


# Method B - edit the gtf file directly - makeTxDbFromGFF failed
----------------------------------------------
# i. Download "Gasterosteus_aculeatus.BROADS1.88.gtf.gz" ftp://ftp.ensembl.org/pub/release-88/gtf/gasterosteus_aculeatus/
# First few lines are: (original included the #! at start of line)
	#!genome-build BROAD S1
	#!genome-version BROADS1
	#!genome-date 2006-02
	#!genebuild-last-updated 2010-05
	
# ii. Read gtf file into data frame

	# Details (from "https://www.rdocumentation.org/packages/refGenome/versions/1.7.0/topics/read.gtf")
	# GTF is an extension of the GFF file format. GTF contains tabled data: Nine columns separated by a tab 
	# delimiter. The last column expands into a list of attributes, separated by a semicolon an exactly one 
	# space. Each attribute consists of a type - value pair which are separated by one empty space. 
	# Enclosing quotation marks (") around attribute values are marks are skipped during import.
	# The first eight columns of the gtf table are fixed. The content is described in the following:
	# seqid: Chromosome identifier. Character.
	# source: Program which generated data.
	# feature: Feature type (e.g. 'exon', 'CDS'). Character.
	# start: Start position of feature (1-based). Integer.
	# end: End position of feature (inclusive). Integer.
	# score: Value between 0 and 1000 ("." for no score). Character.
	# strand: '+', '-' or '.'. Character.
	# frame: 0-2 for coding exons. '.' otherwise. Character.
	# (9th column is a list of attributes)

setwd("~/Desktop")

# Read the gtf file as a data frame
x <- read.table("Gasterosteus_aculeatus.BROADS1.88.gtf.gz", sep = "\t", stringsAsFactors = FALSE, quote = "")
names(x) <- c("seqid","source","feature","start","end","score","strand","frame","attributes")
head(x)
	    # seqid  source     feature start  end score strand frame
	# 1 groupIV ensembl        gene  3396 9380     .      -     .
	# 2 groupIV ensembl  transcript  3396 9380     .      -     .
	# 3 groupIV ensembl        exon  9231 9380     .      -     .
	# 4 groupIV ensembl        exon  7580 7601     .      -     .
	# 5 groupIV ensembl         CDS  7580 7591     .      -     0
	# 6 groupIV ensembl start_codon  7589 7591     .      -     0
	                                                                                                                                                                                                                                                                                                                                                # attributes
	# 1                                                                                                                                                                                                                                  gene_id "ENSGACG00000016217"; gene_version "1"; gene_name "rnf4"; gene_source "ensembl"; gene_biotype "protein_coding";
	# 2                                                                        gene_id "ENSGACG00000016217"; gene_version "1"; transcript_id "ENSGACT00000021436"; transcript_version "1"; gene_name "rnf4"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "rnf4-201"; transcript_source "ensembl"; transcript_biotype "protein_coding";
	# 3       gene_id "ENSGACG00000016217"; gene_version "1"; transcript_id "ENSGACT00000021436"; transcript_version "1"; exon_number "1"; gene_name "rnf4"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "rnf4-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSGACE00000190266"; exon_version "1";
	# 4       gene_id "ENSGACG00000016217"; gene_version "1"; transcript_id "ENSGACT00000021436"; transcript_version "1"; exon_number "2"; gene_name "rnf4"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "rnf4-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSGACE00000190269"; exon_version "1";
	# 5 gene_id "ENSGACG00000016217"; gene_version "1"; transcript_id "ENSGACT00000021436"; transcript_version "1"; exon_number "2"; gene_name "rnf4"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "rnf4-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSGACP00000021395"; protein_version "1";
	# 6                                                       gene_id "ENSGACG00000016217"; gene_version "1"; transcript_id "ENSGACT00000021436"; transcript_version "1"; exon_number "2"; gene_name "rnf4"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "rnf4-201"; transcript_source "ensembl"; transcript_biotype "protein_coding";

nrow(x)
	# [1] 683135

# iii. change "group" to "chr" in seqid
x$seqid <- sub("group", "chr", x$seqid)
grep("group", x$attributes) # 'group' not used anywhere else?
	# integer(0) 

# iv. Convert scaffold coordinates to chrUn coordinates (old assembly)

# -- First hive off the scaffolds from the gtf file
x1 <- x[grep("scaffold", x$seqid), ]
head(x1[,1:8],20)
	             # seqid  source         feature start   end score strand frame
	# 611945 scaffold_27 ensembl            gene 36048 39906     .      +     .
	# 611946 scaffold_27 ensembl      transcript 36048 39906     .      +     .
	# 611947 scaffold_27 ensembl            exon 36048 36615     .      +     .
	# 611948 scaffold_27 ensembl             CDS 36158 36615     .      +     0
	# 611949 scaffold_27 ensembl     start_codon 36158 36160     .      +     0
	# 611950 scaffold_27 ensembl            exon 38741 39035     .      +     .
	# 611951 scaffold_27 ensembl             CDS 38741 39035     .      +     1
	# 611952 scaffold_27 ensembl            exon 39114 39906     .      +     .
	# 611953 scaffold_27 ensembl             CDS 39114 39491     .      +     0
	# 611954 scaffold_27 ensembl      stop_codon 39492 39494     .      +     0
	# 611955 scaffold_27 ensembl  five_prime_utr 36048 36157     .      +     .
	# 611956 scaffold_27 ensembl three_prime_utr 39495 39906     .      +     .
	# 611957 scaffold_27 ensembl            gene 48405 52111     .      -     .
	# 611958 scaffold_27 ensembl      transcript 48405 52111     .      -     .
	# 611959 scaffold_27 ensembl            exon 51946 52111     .      -     .
	# 611960 scaffold_27 ensembl             CDS 51946 52078     .      -     0
	# 611961 scaffold_27 ensembl     start_codon 52076 52078     .      -     0
	# 611962 scaffold_27 ensembl            exon 51074 51192     .      -     .
	# 611963 scaffold_27 ensembl             CDS 51074 51192     .      -     2
	# 611964 scaffold_27 ensembl            exon 50304 50525     .      -     .

table(x1$strand)
	    # -     + 
	# 35422 35769 
	
table(x1$start <= x1$end) # start is always less than or equal to end, no matter the strand
	 # TRUE 
	# 71191 

# there are only up to scaffold_1877 in gtf file, 
#	though "glazerFileS4 NewScaffoldOrder.csv" lists up to 1934
z <- as.integer(sub("scaffold_", "", x1$seqid))
range(z)
	# [1]   27 1877 


# -- Next, read the scaffold positions from "glazerFileS4 NewScaffoldOrder.csv"
chrUnScaffolds <- read.csv("glazerFileS4 NewScaffoldOrder.csv")
chrUnScaffolds <- chrUnScaffolds[chrUnScaffolds$OldChr == "Un", ]
chrUnScaffolds <- chrUnScaffolds[order(chrUnScaffolds$Scaffold), ]
head(chrUnScaffolds)
	    # Scaffold  Length NewChr NewStart   NewEnd NewOrientation OldChr OldStart   OldEnd
	# 156       27 5105261     17  6993892 12099152        reverse     Un        1  5105261
	# 179       37 2648413     21  1792806  4441218        forward     Un  5106262  7754674
	# 66        47 1748180      7 11315283 13063462        reverse     Un  7755675  9503854
	# 42        48 1738986      5  1213557  2952542        forward     Un  9504855 11243840
	# 41        54 1212556      5        1  1212556        forward     Un 11244841 12457396
	# 70        56 1168921      8        1  1168921        forward     Un 12458397 13627317

# -- Loop through the scaffolds and covert the coordinates

Unstart <- rep(0, nrow(x1)) # initialize
Unend	<- rep(0, nrow(x1)) # initialize
for(i in 1:nrow(chrUnScaffolds)){
	# i <- 2
	# Grab data on Un positions from corresponding row of chrUnScaffolds
	scaffold <- paste("scaffold_", chrUnScaffolds$Scaffold[i], sep = "")
	oldstart <- chrUnScaffolds$OldStart[i]
	# oldend   <- chrUnScaffolds$OldEnd[i]
	Unstart[x1$seqid == scaffold] <- x1$start[x1$seqid == scaffold] + oldstart - 1
	Unend[x1$seqid == scaffold] <- x1$end[x1$seqid == scaffold] + oldstart - 1
	}

# Check a few
head(data.frame(x1$start, x1$end, Unstart, Unend)[x1$seqid == paste("scaffold_", "27", sep = ""), ])
head(data.frame(x1$start, x1$end, Unstart, Unend)[x1$seqid == paste("scaffold_", "56", sep = ""), ])
z <- data.frame(x1$start, x1$end, Unstart, Unend)
all(z$x1.end-z$x1.start == z$Unend - z$Unstart)
table(Unstart == 0) # check that all have been modified
table(Unend == 0)

# -- Change the values in the gtf data frame to the new values
x1$start <- Unstart
x1$end <- Unend
x[grep("scaffold", x$seqid), ] <- x1
x$seqid[grep("scaffold", x$seqid)] <- "chrUn"

head(x[grep("Un", x$seqid), 1:8])

# -- Check the seqid values
table(x$seqid)
    chrI    chrII   chrIII    chrIV    chrIX    chrUn     chrV    chrVI   chrVII  chrVIII 
   39435    29925    30418    39279    33981    71191    23608    25774    40386    29264 
    chrX    chrXI   chrXII  chrXIII   chrXIV   chrXIX    chrXV   chrXVI  chrXVII chrXVIII 
   22776    31207    29896    32737    24967    33670    25899    26053    22652    24741 
   chrXX   chrXXI       MT 
   29957    15172      147 

# -- change "MT" to "chrM"
x$seqid[x$seqid == "MT"] <- "chrM"
table(x$seqid)
    # chrI    chrII   chrIII    chrIV    chrIX     chrM    chrUn     chrV    chrVI   chrVII 
   # 39435    29925    30418    39279    33981      147    71191    23608    25774    40386 
 # chrVIII     chrX    chrXI   chrXII  chrXIII   chrXIV   chrXIX    chrXV   chrXVI  chrXVII 
   # 29264    22776    31207    29896    32737    24967    33670    25899    26053    22652 
# chrXVIII    chrXX   chrXXI 
   # 24741    29957    15172 

head(x)
tail(x)

x[39771:39775,1:8] # check whether "200000" isbeing represented as "2e+5"
      seqid  source         feature  start    end score strand frame
39771  chrI ensembl three_prime_utr 206203 206224     .      +     .
39772  chrI ensembl three_prime_utr 206304 206656     .      +     .
39773  chrI ensembl      transcript 200000 206664     .      +     .
39774  chrI ensembl            exon 200000 200079     .      +     .
39775  chrI ensembl             CDS 200009 200079     .      +     0


# v. Write results to new gtf file
	# *to prevent it writing 200000 as 2e+05, change integers to character - DIDN'T WORK
	# x$start <- as.character(x$start)
	# x$end <- as.character(x$end)
# This worked:
options(scipen=999)
write.table(x, file = "gasAcu1.ensembl.88.cleaned.gtf", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# vi. Paste these 4 lines back into "gasAcu1.ensembl.88.cleaned.gtf" by hand
	#!genome-build BROAD S1
	#!genome-version BROADS1
	#!genome-date 2006-02
	#!genebuild-last-updated 2010-05

# vii. Gzip new gtf file
system("gzip gasAcu1.ensembl.88.cleaned.gtf")

