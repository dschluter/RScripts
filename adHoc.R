#!/usr/bin/Rscript

# Edited directly to execute a very specific command
# Takes no arguments.

project <- "BenlimAllMarine"
chrname <- c("chrXXI")
interestingPairs <- list(
	c("paxl", "paxb"), c("pril", "prib"), c("qryl", "qryb"), c("ensl", "ensb"),
	c("marine-pac", "paxb"), c("marine-pac", "paxl"), 
	c("marine-pac", "solitary"),
	c("paxb", "solitary"), c("prib", "solitary"), c("qryb", "solitary"), c("ensb", "solitary")
	)

stepsize <- 500
nsteps.per.window <- 5 	# window size is (nsteps.per.window)*(stepsize), e.g., 5*500 = 2500
windowNmin = 100 # minimum number of good bases in window

g$plotSlidewinInterestingPairsByChr(method = "vara", project, chrname, interestingPairs,
			stepsize = 500, nsteps.per.window = 5, windowNmin = 100, orderChr = TRUE,
			Glazerize = TRUE, scafFile = "glazerFileS4 NewScaffoldOrder.csv")
