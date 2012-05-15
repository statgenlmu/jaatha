.ms.outfile <- function() {
	return( paste(tempdir(),"/ms.out",sep="") )
}

.ms.generateCmd <- function(dm,parameters,outfile=.ms.outfile()){
	.dm.log(dm,"Called .ms.generateCmd()")

	nSample <- dm@sampleSizes
	cmd <- c("ms",sum(nSample),dm@nLoci)                     #repetitions and number of loci
	cmd <- c(cmd,"-I 2",nSample[1],nSample[2])            	 #two populations
	
	for (i in 1:dim(dm@features)[1] ){
		type <- as.character(dm@features[i,"type"])
		#cat(i,":",type,"\n")

		if (type == "mutation"){ 
			if (dm@externalTheta) cmd <- c(cmd,"-t 5") 
			else cmd <- c(cmd,"-t",.ms.getParameter(dm,parameters,type))
		}
		
		if (type == "split") 
			cmd <- c(cmd,"-ej",.ms.getParameter(dm,parameters,type,1),"2 1")
		
		if (type == "migration")
			cmd <- c(cmd,"-m 1 2",.ms.getParameter(dm,parameters,type),
				     "-m 2 1",.ms.getParameter(dm,parameters,type))
		
		if (type == "recombination") 
			cmd <- c(cmd,"-r",.ms.getParameter(dm,parameters,type),dm@seqLength)

		#if (type == "size")
		#	cmd <- c(cmd,"-n",dm@features[i,"population"],
		#		     parameters[i])
		
		#if (type == "sizeChange")
		#	cmd <- c(cmd,"-g",dm@features[i,"population"],
		#		     parameters[i])
	}

	if (outfile != F) cmd <- c(cmd,">",outfile)
	cmd <- paste(cmd,collapse=" ")
	.dm.log(dm,"Generated simulation cmd:",cmd)
	.dm.log(dm,"Finished .ms.generateCmd()")
	return(cmd)
}

.ms.getParameter <- function(dm,parameters,type,population=NA){
	feature <- .getFeature(dm,type,population)
	if ( dim(feature)[1] != 1 ) stop("Error creating ms command")
	if ( feature$lowerRange[1] == feature$upperRange[1] )
		return(feature$lowerRange[1]) 
	else
		return(parameters[feature$parameter])
}

#dyn.load(paste("src/msFile2jsfs",.Platform$dynlib.ext, sep=""))
.ms.getJSFS <- function(dm){
	.dm.log(dm,"Called .ms.getJSFS()")
	resultSize <- (dm@sampleSizes[1]+1)*(dm@sampleSizes[2]+1)
	jsfs <- ( .C("msFile2jsfs",
		as.character(.ms.outfile()),
		as.integer(dm@sampleSizes[1]),
		as.integer(dm@sampleSizes[2]),
		as.integer(dm@nLoci),
		res=integer(resultSize),
	        PACKAGE="jaatha")$res )
	unlink(.ms.outfile())
	.dm.log(dm,"Finished .ms.getJSFS()")
	return(matrix(jsfs,dm@sampleSizes[1]+1,dm@sampleSizes[2]+1))
}

.ms.simSumStats <- function(dm,parameters,sumStatFunc){
	.dm.log(dm,"Called .ms.simSumStats()")

	nSumStats <- length(sumStatFunc(matrix(0,dm@sampleSizes[1],dm@sampleSizes[2])))
        nSims	  <- max(dim(parameters)[1],1)
	sumStats  <- matrix(0,nSims,nSumStats)

	.dm.log(dm,"Simulating",nSumStats,"summary statistics for",nSims,"parameter combination(s)")

	for ( n in 1:nSims ) {
		system(.ms.generateCmd(dm,parameters[n,]))
		sumStats[n,] <- sumStatFunc(.ms.getJSFS(dm))
	}

	.dm.log(dm,"Finished .ms.simSumStats()")
	return(sumStats)
}

