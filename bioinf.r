git <- function (githubfile){
    library(RCurl, quietly = TRUE)
    script <- getURL(paste("https://raw.githubusercontent.com/dschluter/RScripts/master/", 
        githubfile, sep = ""), .opts = list(ssl.verifypeer = FALSE))
    eval(parse(text = script), envir = .GlobalEnv)
	}

bestClust <- function(x = pca.scores[, 1:2], j1 = -1.2, j2 = 1.2, by = 0.2, plot = FALSE, minProb = 0.95){
	# Function to return a list with results from an Mclust cluster analysis using the (hopefully single) best projection
	# Find the best projection on which to classify fish into three genotypes MM, MF, FF
	# Evaluates several different projection vectors and fits to a univariate Gaussian mixture model using Mclust
	# Projections are vectors starting at the origin and going to (1,-1.2), (1, -0.8), (1, -0.6), (1, -0.4), (1, -0.2), or (1, 0) by default
	# The best projection to be used for classification is the one having highest mean posterior probability of classification
	# To understand projection calculation, refer to figure at https://en.wikipedia.org/wiki/Scalar_projection
	# 	b is the vector onto which points are being projected
	# 	x contains the data (usually pc1 and pc2)
	# 	proj is a1 in Wiki figure, the shadow of the vector x on b
	# 	distance is a2 in the Wiki figure, the vector "rejection", i.e., the perpendicular distance between x and proj.

	library(mclust)

	results <- list()

	pca.scores <- as.data.frame(pca.scores[, 1:2])
	colnames(pca.scores) <- c("PC1", "PC2")

	if(plot) pdf()
	# Convert both axes to have range [0,1]
	x$PC1 <- (x$PC1 - min(x$PC1)) / (max(x$PC1) - min(x$PC1)) # convert to (0,1)
	x$PC2 <- (x$PC2 - min(x$PC2)) / (max(x$PC2) - min(x$PC2)) # convert to (0,1)
	
	# Fit the univariate mclust model to data on each projection in turn
	k <- 1
	for(j in seq(j1, j2, by = by)){
		# j <- -0.6
		b <- c(1, j)
		blength <- sqrt(sum(b^2))
		proj <- as.matrix(x) %*% cbind(b) / blength
		xlength <- sqrt(x$PC1^2 + x$PC2^2)
		distance <- unname(sqrt(xlength^2 - proj^2))
		# z <- densityMclust(proj, modelNames = "E", G = 3, plot = FALSE)
		results[[k]] <- densityMclust(proj, modelNames = "E", G = 3, plot = FALSE)
		names(results)[k] <- paste0("proj", j)
		k <- k + 1
		}
	
	# Grab the average posterior probability for each model
	AvgPostProb <- sapply(results, function(z){
		# maxposterior calculates the posterior probability for the class into which the fish was placed
		# avgPosterior averages this across all fish to measure overall confidence of assignments
		maxposterior <- apply(z$z, 1, max)
		return(mean(maxposterior))
		})
	return(results[which(AvgPostProb == max(AvgPostProb))])
	}


ecopeak2vcf <- function(vcf.in, vcf.out, chrname, start, end, name = NULL, 
	fastaname = "stickleback_v5_assembly.fa", samples = NULL, snptable = FALSE){
	# Extracts genotype data for the SNP range from main vcf file and creates a new vcf file and index file
	# vcf.in is the main vcf file containing snps
	# vcf.out is the new vcf file to be created, is bgz compressed
	# fastaname is the name of the genome fasta file
	# ecopeak.ranges is a GRanges object with seqname, start, end information
	# If samples is provided as a vector of ID's 
	#	e.g., c("Marine-Pac-LittleCampbellRiver_BC.LCR-45.RV", "Solitary-KirkLake.KL-5.RV"), 
	#	then the sample list is subsetted
	# If snptable = TRUE, then a snp table file is also generated
	library(VariantAnnotation)
	
	# Make GRanges object
	ecopeak.ranges <- makeGRangesFromDataFrame(data.frame(seqnames = chrname, start = start, end = end))
	if(!is.null(name)) ecopeak.ranges$name <- name
	
	if(is.null(samples)) 
		param <- ScanVcfParam(fixed = "ALT", geno = c("GT", "PL"), info = NA, which = ecopeak.ranges) else
		param <- ScanVcfParam(fixed = "ALT", geno = c("GT", "PL"), info = NA, which = ecopeak.ranges, samples = samples)
	vcf <- readVcf(vcf.in, genome = fastaname, param = param)
	
	# Remove spanning deletions (*), or null alleles ( ) which are included still.
	# no spanning deletions here, just blank alleles
	ALT <- CharacterList(alt(vcf))
	has.null.allele <- sapply(ALT, function(x){"" %in% x})
	has.span.deletion <- sapply(ALT, function(x){"*" %in% x})
	z <- sapply(ALT, length)
	table(has.null.allele)
		# FALSE  TRUE 
		 # 4737    72
 	table(has.span.deletion)
	
	# Delete spanning deletions (*), or null alleles ( )
	vcf <- vcf[!has.span.deletion & !has.null.allele]
	
	print("Dimensions of genotype array:")
	print(dim(geno(vcf)$GT))
		# [1] 4737  1158
	
	vcf.out <- sub(".gz$", "", vcf.out) # remove trailing "gz" if present
	writeVcf(vcf, file = vcf.out, index = TRUE) # makes a bgz file
	
	if(snptable){
		snps <- genotypeToSnpMatrix(vcf) # rows = samples, columns = snps
		z <- apply(snps$genotypes@.Data, 2, function(x) as.integer(x))
		rownames(z) <- colnames(geno(vcf)$GT)
		snps <- as.data.frame(z)
		snps.out <- sub(".vcf$", ".snptable.csv", vcf.out)
		write.csv(snps, file = snps.out, row.names = TRUE)
		}
	}

