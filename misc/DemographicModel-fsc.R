.fsc.profileFile <- paste(tempdir(),"/fsc_profile.tpl",sep="")
.fsc.paramFile 	 <- paste(tempdir(),"/fsc_pars.def",sep="")
.fsc.workdir 	 <- tempdir()
.fsc.jsfsfile 	 <- paste(tempdir(),"/fsc_profile_jointDAFpop1_0.obs",sep="")
.fsc.cleanDir 	 <- paste(tempdir(),"/fsc_*",sep="")

.fsc.cmd <- paste("cd",.fsc.workdir,";","fastsimcoal -dsqx -n1 -t",
		  .fsc.profileFile,"-f",.fsc.paramFile,"> /dev/null")

.fsc.generateCmd <- function(dm,relParamValue,nLoci){
	if (dm@profiled == F) .fsc.createProfile(dm)
	.fsc.createParFile(dm,relParamValue)
	return(.fsc.cmd)
}

.fsc.getJSFS <- function(dm){
	jsfs  <-  as.matrix(read.table(.fsc.jsfsfile,skip=1))
	unlink(.fsc.cleanDir)
	return(jsfs)
}

.fsc.createProfile <- function(dm){
	.fsc.appendToProfile("//Number of populations",append=F)
	.fsc.appendToProfile("2")

	#Effective population sizes
	.fsc.appendToProfile("//Effective population sizes" )
	.fsc.appendToProfile( "Ne1" )
	.fsc.appendToProfile( "Ne2" )
	
	#Sample sizes
	.fsc.appendToProfile( "//Sample Sizes" )
	.fsc.appendToProfile( dm@sampleSizes[1] )
	.fsc.appendToProfile( dm@sampleSizes[2] )

	#Growth rates
	.fsc.appendToProfile( "//Growth rates" )
	.fsc.appendToProfile( "0" )
	.fsc.appendToProfile( "0" )

	#Migration
	.fsc.appendToProfile( "//Number of migration matrices" )
	if (.getParameter(dm,"migration") == 0) {
		.fsc.appendToProfile("0")
	} else {
		.fsc.appendToProfile("2")
		.fsc.appendToProfile("//Matrix 0\n0 0\n0 0")
		.fsc.appendToProfile("//Matrix 1\n0 ",.getParameter(dm,"migration"),"\n",
		    		 	.getParameter(dm,"migration")," 0")
	}

	#Historical event
	.fsc.appendToProfile("//Historical event")
	.fsc.appendToProfile(.getNumOfHistoricalEvent(dm)," historical event")

	for (i in 1:dim(dm@features)[1]){
		type <- dm@features$type[i]

		if (type == "split") 
			.fsc.appendToProfile(.getParameter(dm,"split",1),
					     " ",1," ",0,
					     " ",1," ",1,
					     " ",0," ",0)
		
		if (type == "migration")
			.fsc.appendToProfile(0,
					     " ",1," ",0,
					     " ",1," ",1,
					     " ",0," ",0)
	}

	#Number of loci
	.fsc.appendToProfile("//Number of loci")
	.fsc.appendToProfile(dm@nLoci," ",0)

	#Blocks
	.fsc.appendToProfile("//Number of Block per Loci")
	.fsc.appendToProfile("1")
	.fsc.appendToProfile("//Blocks")
	#.fsc.appendToProfile("DNA"," ",dm@seqLength,
	#	  		   " ",.getParameter(dm,"recombination"),
	#    	  		   " ",.getParameter(dm,"mutation"),
	#	  		   " ",dm@tsTvRatio)
	.fsc.appendToProfile("SNP"," ",dm@seqLength,
		  		   " ",.getParameter(dm,"recombination"),
	    	  		   " ",.getParameter(dm,"mutation"))

	#dm@profiled <<- T
}

.fsc.createParFile <- function(dm,absValues){
	if (!is.numeric(absValues)) stop("absValues must be numeric!")
	if (!is.matrix(absValues)) absValues <- matrix(pars,1)
	colnames(absValues) <- dm@parameters

	# Get size of population 1
	Ne1Par <- .getParameter(dm,"effPopSize",1)
	if (Ne1Par == 0) {
		Ne1 <- 10000
		absValues <- cbind(absValues,Ne1)
	} else if ( is.numeric(Ne1Par) ) {
		Ne1 <- Ne1Par
		absValues <- cbind(absValues,Ne1)
	} else { 
		absValues[,'Ne1'] <- round(absValues[,'Ne1'],0)
	}
	
	# Get size of population 2
	Ne2Par <- .getParameter(dm,"effPopSize",2)
	if (Ne2Par == 0) {
		Ne2 <- 10000
		absValues <- cbind(absValues,Ne2)
	} else if ( is.numeric(Ne2Par) ) {
		Ne2 <- Ne2Par
		absValues <- cbind(absValues,Ne2)
	} else { 
		absValues[,'Ne2'] <- round(absValues[,'Ne2'],0)
	}

	#Scale Parameters:
	# * Parameters that are scaled with 4*Ne1
	pars <- dm@features$parameter[dm@features$type == "migration"]
	absValues[,pars] <- absValues[,pars] / (4 * absValues[,'Ne1'])

	# * Parameters that are scaled with N*Ne1*seqLength
	pars <- dm@features$parameter[
			  dm@features$type == "mutation"
			| dm@features$type == "recombination" 
		]
	absValues[,pars] <- absValues[,pars] / (4 * (absValues[,'Ne1']+absValues[,'Ne2']) * dm@seqLength)
	
	# * times are scaled with 1/(4*Ne1) 
       	pars <- dm@features$parameter[dm@features$type == "split"]
	absValues[,pars] <- round( absValues[,pars] * 4 * absValues[,'Ne1'], 0 )

	write.table(absValues,file=.fsc.paramFile,row.names=F,quote=F)
	return(absValues)
}

.fsc.appendToProfile <- function(...,append=T){
	cat(...,"\n",file=.fsc.profileFile,sep="",append=append)
}

.fsc.simSumStats <- function(dm,parameters,sumStatFunc=.fsc.calcSumStats){
	system(.fsc.generateCmd(dm,parameters))
	lines <- scan(.fsc.jsfsfile, what="character", sep="\n", quiet=T)

	nSumStats <- length(sumStatFunc(matrix(0,dm@sampleSizes[1]+1,dm@sampleSizes[2]+1)))
        nSims	  <- max(dim(parameters)[1],1)
	print(c(nSims,nSumStats))
	sumStats  <- matrix(0,nSims,nSumStats)

	for ( n in 1:nSims ) {
		i <- n+1 + (n-1) * (dm@sampleSizes[2] + 2)
		tC <- textConnection(lines[i:(i+dm@sampleSizes[2]+1)])
		jsfs <- t(as.matrix(read.table(tC,header=T)))
		close(tC)
		sumStats[n,] <- sumStatFunc(jsfs)
		unlink(.fsc.cleanDir)
	}

	return(sumStats)
}

.fsc.calcSumStats <- function(jsfs) {
  n <- nrow(jsfs)
  m <- ncol(jsfs)
  c(sum(jsfs[1,2:3]),
    sum(jsfs[2:3,1]),
    sum(jsfs[1,4:(m-3)]),
    sum(jsfs[4:(n-3),1]),
    sum(jsfs[1,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),1]),
    sum(jsfs[2:3,2:3]),
    sum(jsfs[2:3,4:(m-3)]),
    sum(jsfs[4:(n-3),2:3]),
    sum(jsfs[(n-2):(n-1),4:(m-3)]),
    sum(jsfs[4:(n-3),(m-2):(m-1)]),
    sum(jsfs[2:3,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),2:3]),
    sum(jsfs[4:(n-3),4:(m-3)]),
    sum(jsfs[(n-2):(n-1),(m-2):(m-1)]),
    jsfs[1,m],
    jsfs[n,1],
    sum(jsfs[n,2:3]),
    sum(jsfs[2:3,m]),
    sum(jsfs[n,4:(m-3)]),
    sum(jsfs[4:(n-3),m]),
    sum(jsfs[n,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),m]) )
}
