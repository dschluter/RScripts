git <- function(githubfile){
	# Sources the contents of a script file located at github
	# Use as git("genome.r")
	library(RCurl)
	script <- getURL( paste("https://raw.githubusercontent.com/dschluter/genomeScripts/master/", githubfile, sep="") )
	eval(parse(text = script), envir = .GlobalEnv)
	}

month2integer <- function(x){
	z1 <- casefold(substr(x, 1, 3)) # round to 3 digits
	months <- c("january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december")
	months.abb <- casefold(substr(months, 1, 3))
	z2 <- recode(z1, months.abb, 1:12)
	}

sizeadjust <- function(x, size, group){
	# Size-adjust traits contained in data frame "x" (no other traits allowed in "x")
	# Size corrections assume a common slope among levels of group variable and no need to take logs
	# Traits are adjusted to a fish of size msize = mean(size)
	# nb: Can't just apply the function "residuals" to the model object because of missing values
	# i <- 1
	if(!is.factor(group)) group <- factor(group)
	msize <- mean(size, na.rm=TRUE)
	results <- x # initialize
	for(i in 1:ncol(x)){
		z <- lm(x[,i] ~ size + group)
		y1 <- predict( z, newdata = data.frame(size = rep(msize, length(group)), group) ) # predicted values
		y2 <- x[,i] - predict(z, newdata = data.frame(size, group)) # residuals
		y3 <- y1 + y2 # adjusted values
		results[,i] <- y3
		}
	return(results)
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

digitize<-function(x,xlim=c(0,1),ylim=c(0,1),se=FALSE){
# function to convert digitized points into actual scale
# x is assumed to be in data table format
# first three points are assumed to be the limits of the scale, 
# digitized in order as max(Y), [min(X),min(Y)], max(X)
# xlim and ylim are the values on the graph corresponding to the limits of the scale
# if se=TRUE then it is assumed that the last half of the cases are upper or lower limits
# corresponding to the means in the first half
	xlen<-xlim[2]-xlim[1]
	ylen<-ylim[2]-ylim[1]
	x[,1]<-x[,1]-x[2,1]
	x[,2]<-x[,2]-x[2,2]
	x[,1]<-x[,1]/x[3,1]
	x[,2]<-x[,2]/x[1,2]
	x<-x[-c(1:3),]
	x[,1]<-x[,1]*xlen + xlim[1]
	x[,2]<-x[,2]*ylen + ylim[1]
	if(se){
		x[,3]<-x[(1+nrow(x)/2):nrow(x),2]
		x<-x[1:(nrow(x)/2),]
		x[,3]<-abs(x[,3]-x[,2])
		names(x)<-c("x","y","se.y")
		}
	row.names(x)<-as.character(seq(1:nrow(x)))
	return(x)
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

writeclip <- function(text = "http://ezproxy.library.ubc.ca/login?url=http://", ...){
	dat <- pipe("pbcopy", "w")
	write(text, dat)
	on.exit(close(dat))
	return(invisible())
	}

scanclip <- function(...){
	# use "scan" to scan the clipboard
	# ... are options passed to scan, eg
	#   what = character()
	#   sep=","
	#   etc
	dat <- pipe("pbpaste")
	x <- scan(dat, ...)
	on.exit(close(dat))
	return(x)
	}

readclip <- function(...){
	# use "read.table" to scan the clipboard
	# ... are options passed to read.table
	dat <- pipe("pbpaste")
	x <- read.table(dat, ...)
	#on.exit(close(dat)) # causes an error
	return(x)
	}

clip <- function(primary = TRUE, ...) {
	# Modified from "clipsource" in the svMisc package
	# Source data from the clipboard, manage clipboard correctly depending
	# on the OS
	if (.Platform$OS.type == "windows") { # Windows OS
		data <- file("clipboard")
	} else if (grepl("mac.binary",.Platform$pkgType)) {	# Mac OS
		data <- pipe("pbpaste")
	} else {	# Must be Linux/Unix
		if (primary) {
			data <- file("X11_clipboard")
		} else {
			data <- file("X11_secondary")
		}
	}
	on.exit(close(data))
	# Invoke source() with the data from the clipboard
	res <- source(data, ...)
	return(invisible(res))	
} 

recode<-function (..., ret = c("numeric", "factor"), none = if (ret == "numeric") 0 else "none", na){
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



ls2<-function(type=NULL,envir=.GlobalEnv){
	# Selectively lists objects of specified type
	# To remove all of any type, use "rm(list=ls2(type="character"))"
	if(is.null(type)){
		type<-select.list(c("character","data.frame","function","list",
			"matrix","numeric","vector"))
		}
	if(type=="character")
		return(ls(envir=envir)[sapply(ls(envir=envir), function(x) is.character(get(x)))])
	if(type=="data.frame")
		return(ls(envir=envir)[sapply(ls(envir=envir), function(x) is.data.frame(get(x)))])
	else if(type=="function")
		return(ls(envir=envir)[sapply(ls(envir=envir), function(x) is.function(get(x)))])
	else if(type=="list")
		return(ls(envir=envir)[sapply(ls(envir=envir), function(x) is.list(get(x)))])
	else if(type=="matrix")
		return(ls(envir=envir)[sapply(ls(envir=envir), function(x) is.matrix(get(x)))])
	else if(type=="numeric")
		return(ls(envir=envir)[sapply(ls(envir=envir), function(x) is.numeric(get(x)))])
	else if(type=="vector")
		return(ls(envir=envir)[sapply(ls(envir=envir), function(x) is.vector(get(x)))])
}

logit<-function(x)	{
	log(x/(1. - x))
}

antilogit<-function(x){
	exp(x)/(exp(x) + 1)
}

compute<-function(x, strata, FUN,...){
	# To compute statistics on columns of a data frame or matrix x.
	# strata is the category variable defining the strata.
	# FUN is the function to be carried out (not in quotes).
	# ... is used to pass additional info to tapply, such as "na.rm=TRUE"
	if(is.null(ncol(x))) {
		result <- cbind(tapply(x, strata, FUN,...))
	}
	else {
		result <- vector()
		for(j in 1.:ncol(x)) {
			result <- cbind(result, tapply(x[, j], strata, FUN,...))
		}
	}
	if(!is.null(dimnames(x)[2.]))
		dimnames(result)[2.] <- dimnames(x)[2.]
	return(as.data.frame(result))
}

#view(x)<-function(x){x<-edit(x)}

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

view<-function(x){
	library(tcltk)
	library(relimp)
	if(!is.data.frame(x) & !is.matrix(x)) stop("Object must be matrix or data frame")
	showData(x,font='Courier 10',title=deparse(substitute(x)))
	invisible()
}

graph<-function(square=FALSE, cex = 1., las = 1., pch = 1., mkh = 0.05, pty1 = "m",...){
	if(square) pty1 <- "s"
	dev.set(which=1)
	par(pty = pty1, cex = cex, las = las, pch = pch, mkh = mkh,...)
	invisible()
}

#graph<-function(square=FALSE, cex = 1., mar = c(3., 3., 1., 1.), 
#		mgp = c(2., 1., 0.), las = 1., pch = 1., mkh = 0.05){
#	pty1 <- "m"
#	if(square) pty1 <- "s"
#	dev.set(which=1)
#	par(pty = pty1, cex = cex, mar = mar, mgp = mgp, las = las, pch = pch, mkh = mkh)
#	invisible()
#}

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

#recode<-function # See character.s

# plot.regression<-function(x,y,treatment,residuals=FALSE,ancova=FALSE,title=NULL){
	# # Plots data separately for each treatment
	# # Adds regression lines, optionally separately for each treatment
	# # residuals=T plots residuals from the single regression of y on x (regardless of ancova)
	# # ancova=T plots separate regression lines for each treatment, otherwise a single regression line is plotted
	# # example: 
	# # plot.regression(x$standlength,x$bodydepth,x$treatment)
# #	if(newplot) windows(record=TRUE)
	
	# if(residuals){
		# z<-lm(y~x)
		# plot(x,residuals(z),type="n",xlab=deparse(substitute(x)),ylab="Residuals",las=1)
		# for(i in 1:length(groups)){
			# points(x[treatment==groups[i]],residuals(z)[treatment==groups[i]],pch=i)
			# }
		# lines(c(x1,x2),c(0,0))
		# title(title)
		# }

	# else{
		# plot(x,y,type="n",xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),las=1)
		# title(title)
		# groups<-unique(treatment)
		# for(i in 1:length(groups)){
			# points(x[treatment==groups[i]],y[treatment==groups[i]],pch=i)
			# }
		# if(ancova){
			# for(i in 1:length(groups)){
				# z<-lm(y[treatment==groups[i]]~x[treatment==groups[i]])
				# x1<-min(x[treatment==groups[i]])
				# x2<-max(x[treatment==groups[i]])
				# intercept<-coef(z)[1]
				# slope<-coef(z)[2]
				# lines(c(x1,x2),intercept+slope*c(x1,x2))
				# }
			# }
		# else{
			# z<-lm(y~x)
			# x1<-min(x)
			# x2<-max(x)
			# intercept<-coef(z)[1]
			# slope<-coef(z)[2]
			# lines(c(x1,x2),intercept+slope*c(x1,x2))
			# }
		# }
# }

plot.robust<-function(x,y,treatment,residuals=FALSE,ancova=FALSE,title=NULL,method="lqs",...){
	# residuals=T plots residuals from the single regression of y on x (regardless of ancova)
	# ancova=T plots separate regression lines for each treatment, otherwise a single regression line is plotted
	# The dots (...) pass extra info to the "lqs" function
	# example: 
	# plot.regression(x$standlength[x$family=="ipe164"],x$bodydepth[x$family=="ipe164"],x$treatment[x$family=="ipe164"])
	library(MASS)
	if(residuals){
#		z<-ltsReg(y~x,alpha=0.50)
		z<-lqs(y~x,method=method,...)
		plot(x,residuals(z),type="n",xlab=deparse(substitute(x)),ylab="Residuals",las=1)
		for(i in 1:length(groups)){
			points(x[treatment==groups[i]],residuals(z)[treatment==groups[i]],pch=i)
			}
		lines(c(x1,x2),c(0,0))
		title(title)
		}

	else{
		plot(x,y,type="n",xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),las=1)
		title(title)
		groups<-unique(treatment)
		for(i in 1:length(groups)){
			points(x[treatment==groups[i]],y[treatment==groups[i]],pch=i)
			}
		if(ancova){
			for(i in 1:length(groups)){
#				z<-ltsReg(y[treatment==groups[i]]~x[treatment==groups[i]],alpha=0.50)
				z<-lqs(y[treatment==groups[i]]~x[treatment==groups[i]],method=method,...)
				x1<-min(x[treatment==groups[i]])
				x2<-max(x[treatment==groups[i]])
				intercept<-coef(z)[1]
				slope<-coef(z)[2]
				lines(c(x1,x2),intercept+slope*c(x1,x2))
				}
			}
		else{
#			z<-ltsReg(y~x,alpha=0.50)
			z<-lqs(y~x,method=method,...)
			x1<-min(x)
			x2<-max(x)
			intercept<-coef(z)[1]
			slope<-coef(z)[2]
			lines(c(x1,x2),intercept+slope*c(x1,x2))
			}
		}
}

