sizeadjust.lm <- function(trait, size, group){
	# My size-correction function using linear model
    # Size corrections assume a common slope among levels of group variable.
    # Traits are adjusted to a fish of size msize = mean(size)
    if(!is.factor(group)) group <- factor(group)
    
    msize <- mean(size, na.rm = TRUE)
    z <- lm(trait ~ size + group, na.action = na.exclude)
    
    # Trait value predicted at the mean size
    y1 <- predict(z, newdata = 
                data.frame(size = rep(msize, length(group)), group) )
                
    # residuals
    y2 <- trait - predict(z, newdata = data.frame(size, group))
    
    # adjusted values
    y3 <- y1 + y2
    return(y3)
    }

sizeadjust.cpc <- function(trait, size, group){
	# My size-correction function using common principal components
	# Size corrections assume a common slope among levels of group variable.
	# Trait is adjusted to have the same mean after adjustment as before
	# trait <- x$grl1.mm; size = x$stdl_new
	library(cpcbp)
	f <- factor(group)
	mtrait <- mean(trait, na.rm = TRUE)
	dat <- data.frame(trait, size)
	    
	nomissing <- which(apply(dat, 1, function(x){all(!is.na(x))}))
	x1 <- dat[nomissing, ]
	f1 <- f[nomissing]
	xlist <- split(x1, f1)
	n <- sapply(xlist, nrow)
	covlist <- lapply(xlist, cov)
	z <- new_cpc(matList = covlist, n = n)
	xadj <- bpfun(dat, f, z[,1], center = FALSE)
	colnames(xadj) <- names(dat)
	trait.adj <- xadj[,1] - mean(xadj[,1], na.rm = TRUE) + mtrait
		return(trait.adj)
	}

sizeadjust.rma <- function(trait, size, group){
    library(dplyr)
	# My size-correction function for reduced major axis regression RMA
	# (equivalent to CPC if traits are standardized so using correlation matrix)
    # Size corrections assume a common slope among levels of group variable.
    # Traits are adjusted to a fish of size msize = mean(size)
    # trait <- x$grl1.mm; size = x$stdl_new; group = x$sexpoplake
    f <- factor(group)

    # Mean size
    msize <- mean(size, na.rm = TRUE)
    
    # Assemble traits and groups
    dat <- data.frame(trait, size)
    nomissing <- which(apply(dat, 1, function(x){all(!is.na(x))}))
    x1 <- dat[nomissing, ]
    f1 <- f[nomissing]
    xlist <- split(x1, f1)
    n <- sapply(xlist, nrow)
    covlist <- lapply(xlist, cov)
    
    # Calculate variances and pool them across groups
    vars <- do.call("rbind", lapply(covlist, diag))
    var_trait <- weighted.mean(vars[, "trait"], w = n - 1)
    var_size <- weighted.mean(vars[, "size"], w = n - 1)
    b <- sqrt(var_trait/var_size)

    # Calculate predicted values and residuals
    # meansize and meantrait are calculated separately for each group
    dat <- data.frame(dat, f)
    dat <- mutate(group_by(dat, f), 
        meansize = mean(size, na.rm = TRUE), 
        meantrait = mean(trait, na.rm = TRUE))
    dat$pred.rma <- b * (dat$size - dat$meansize) + dat$meantrait
    dat$resid.rma <- dat$trait - dat$pred.rma
    
    # calculate size-adjusted values
    trait.adj <- b * (msize - dat$meansize) + dat$meantrait + dat$resid.rma
    return(trait.adj)
    }
    
time2hr <- function(x){
	# x must be a character
	# library(chron) - doesn't accept hours > 24, eg "28:26:12"
	# time2hr("3-03:03:26")
	# time2hr("07:10:51")
	# time2hr("00:43:48")
	# time2hr("04:31.056")
	# time2hr("01:29:36")
	# if(!is.character(x)) x <- deparse(substitute(x))
	hrs <- 0
	ncolon <- sapply(regmatches(x, gregexpr(":", x)), length) # No. of ":" in character
	x1 <- unlist(strsplit(x, split="-"))
	if(length(x1) == 2){
		hrs <- hrs + 24 * as.integer(x1[1])
		x1 <- x1[-1]
		}
	x2 <- unlist(strsplit(x1, split="[.]"))	
	if(length(x2) == 2) x2 <- x2[1] # drop the fractions of seconds
	
	if(ncolon == 1) x2 <- paste("00:", x2, sep = "")
	x3 <- as.integer(unlist(strsplit(x2, split = ":")))
	hrs <- hrs + x3[1] + x3[2]/60 + x3[3]/3600
	return(hrs)
	}
		
is_number <- function(A){ !grepl("[^0-9]+", A) & nchar(A) > 0}

curl <- function(Rfile){
	system(paste("curl -o " , Rfile, " https://raw.githubusercontent.com/dschluter/RScripts/master/", Rfile, sep = ""))
	}

showq <- function(){system("squeue -u schluter")}
squeue <- function(){system("squeue -u schluter")}
scontrol <- function(jobid){system( paste("scontrol show jobid -dd", jobid) )}
scancel <- function(jobid){system( paste("scancel", jobid) )}
scancelAllJobs <- function(){system("scancel -u schluter")}
scancelAllPendingJobs <- function(){system("scancel -t PENDING -u schluter")}

chunpair<-function(a,sep=","){
# Returns c("a") if "a" does not contain sep, otherwise returns c("a1","a2")
# ***This routine fails if a is a vector***
	a<-as.character(a)
	z<-unlist(strsplit(a,split=sep))
	if(length(z)!=2) return(a)
	else {return(z)}
}

month2integer <- function(x){
	z1 <- casefold(substr(x, 1, 3)) # round to 3 digits
	months <- c("january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december")
	months.abb <- casefold(substr(months, 1, 3))
	z2 <- recode(z1, months.abb, 1:12)
	}

nice <- function(x){
	if(is.data.frame(x)) write.table(x, sep = " ", row.names = FALSE, quote = FALSE)
	else if(is.list(x)) for(i in 1:length(x)){
		cat("$",names(x)[i],"\n", sep="")
		print(x[[i]])
		}
	else if(is.vector(x)) sapply(1:length(x), function(i){ cat(x[i], "\n") })
	else stop("Not a data frame, list, or vector")
	return(invisible())
	}

nicelist <- function(z){
	# Function to print out list output (e.g., bio300$anova.oneway) for easy copy and paste to .r file
	for(i in 1:length(z)){
		cat("$",names(z)[i],"\n", sep="")
		print(z[[i]])
		}
	}
	
nicevector <- function(x){
	if(!is.vector(x)) stop("not a vector\n")
	z <- sapply(1:length(x), function(i){
		cat(x[i], "\n")
		})
	return(invisible())
	}
	
niceframe <- function(x){
	if(!is.data.frame(x)) stop("not a data frame\n")
	write.table(x, sep = " ", row.names = FALSE, quote = FALSE)
	return(invisible())
	}

names.quotes <- function(...){
	# use as names.quotes(c(a,b,c))
	#names.quotes <- function(x = c(a,b,c)){  # alternative
	x <- deparse(substitute(...), width.cutoff=500)
	x <- gsub("c(", "", x, fixed=TRUE)
	x <- gsub(")", "", x, fixed=TRUE)
	x <- unlist(strsplit(x, split="[ ,]+"))
	return(x)
	}
addquotes <- names.quotes # synonym

integer2binary <- function(x){
        xi <- as.integer(x)
        if(any(is.na(xi) | ((x-xi)!=0)))
                print(list(ERROR="x not integer", x=x))
        N <- length(x)
        xMax <- max(x)
        ndigits <- (floor(logb(xMax, base=2))+1)
        Base2 <- array(NA, dim=c(N, ndigits))
        for(i in 1:ndigits){#i <- 1
                Base2[, ndigits-i+1] <- (x %% 2)
                x <- (x %/% 2)
        }
        if(N ==1) Base2[1, ] else Base2
} 

logit<-function(x){
	log(x/(1. - x))
}

antilogit<-function(x){
	exp(x)/(exp(x) + 1)
}

tally<-function(x, weights=NULL){
	# Counts up cases in different categories as a prelude to log-linear analysis
	# "x" is a data frame of categorical variables ONLY (or discrete variables)
	# The function counts the number of cases in each combination of categories
	# "weights" refers to a frequency variable if the data are already partly tallied by category,
	# otherwise the weight is assumed to be 1 for each row.
	x<-na.omit(as.data.frame(x))
	if(is.null(weights)){
		weights<-rep(1,nrow(x))
		}

	y<-as.list(x)
	# adjust the levels of y in case some are missing
	for(i in 1:length(y)){
		y[[i]]<-as.factor(as.character(y[[i]]))
		}

	z<-tapply(weights,y,sum)
	z<-as.vector(z)
	z[is.na(z)]<-0
	
	vars<-list(length(y))
	n<-length(z)
	levs<-unlist(lapply(y,function(x){length(levels(x))}))
	levs<-c(1,levs[-length(levs)])
	for(i in 1:length(y)){
		k<-prod(levs[1:i])
		vars[[i]]<-rep(levels(y[[i]]),each=k,length.out=n)
		}
	do.call("cbind.data.frame",vars)
	z<-cbind.data.frame(vars,z)
	names(z)<-c(names(y),"frequency")
	return(z)
}

matchcol<-function(matchto, matchfrom, fromvar){
	# Gets the position in "matchfrom" of the elements in the single vector variable "matchto"
	# Then "fromvar" is the variable/array in "matchfrom" corresponding to those positions
	# At this point matchto and fromvar will have the same length
	# normally after running this macro you will want to one of the following:
	#  1) bind "matchto" to result
	#  2) replace the variable "matchfrom" in result with the elements of "matchto" to eliminate NA's
	position <- match(matchto, matchfrom)
	if(is.null(ncol(fromvar))){
		result <- fromvar[position]
		}
	else{
		result <- fromvar[position,]
		}
	return(result)
}

recode <- function(..., ret = c("numeric", "factor"), none = if (ret == "numeric") 0 else "none", na){
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

