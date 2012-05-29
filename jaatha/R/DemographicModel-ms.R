dm.mscmd <- "ms"

.ms.outfile <- function() {
	return( paste(tempdir(),"/ms.out",sep="") )
}

.ms.generateCmd <- function(dm,parameters,outfile=.ms.outfile()){
	.dm.log(dm,"Called .ms.generateCmd()")

	nSample <- dm@sampleSizes
	cmd <- c(sum(nSample),dm@nLoci)                     #repetitions and number of loci
	cmd <- c(cmd,"-I 2",nSample[1],nSample[2])            	 #two populations
	
	tau <-  .ms.getParameter(dm,parameters,"split")

	for (i in 1:dim(dm@features)[1] ){
		type <- as.character(dm@features[i,"type"])
		pop <- dm@features[i,"population"]
		param <- .ms.getParameter(dm,parameters,type,pop)

		if (type == "mutation"){ 
			if (dm@externalTheta) cmd <- c(cmd,"-t 5") 
			else cmd <- c(cmd,"-t",param)
		}
		
		if (type == "split") 
			cmd <- c(cmd,"-ej",param,"2 1")
		
		if (type == "migration")
			cmd <- c(cmd,"-m 1 2",param,"-m 2 1",param)
		
		if (type == "recombination") 
			cmd <- c(cmd,"-r",param,dm@seqLength)

		if (type == "splitSize"){
		        cmd <- c(cmd,"-g",pop,log(.ms.getPresentSize(dm,parameters,pop)/param)/tau)  
			cmd <- c(cmd,"-eN",tau,1+param)               
		}

		if (type == "presentSize"){
			cmd <- c(cmd,"-n",pop,param)
		}
	}

	#if (outfile != F) cmd <- c(cmd,">",outfile)
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

.ms.getPresentSize <- function(dm,parameters,population){
	feature <- .getFeature(dm,"presentSize",population)
	if ( dim(feature)[1] != 1 ) return(1)
	if ( feature$lowerRange[1] == feature$upperRange[1] )
		return(feature$lowerRange[1]) 
	else
		return(parameters[feature$parameter])
}

.ms.getJSFS <- function(dm,msOut){
	.dm.log(dm,"Called .ms.getJSFS()")
	resultSize <- (dm@sampleSizes[1]+1)*(dm@sampleSizes[2]+1)
	jsfs <- ( .C("ms2jsfs",
		as.character(msOut),
		as.integer(dm@sampleSizes[1]),
		as.integer(dm@sampleSizes[2]),
		as.integer(dm@nLoci),
		res=integer(resultSize),
	        PACKAGE="jaatha")$res )
	#unlink(.ms.outfile())
	.dm.log(dm,"Finished .ms.getJSFS()")
	return(matrix(jsfs,dm@sampleSizes[1]+1,dm@sampleSizes[2]+1))
}

.ms.simSumStats <- function(dm,parameters,sumStatFunc){
	.dm.log(dm,"Called .ms.simSumStats()")

	nSumStats <- length(sumStatFunc(matrix(0,dm@sampleSizes[1],dm@sampleSizes[2])))
        nSims	  <- max(dim(parameters)[1],1)
	sumStats  <- matrix(0,nSims,nSumStats)

	.dm.log(dm,"Simulating",nSumStats,"summary statistics for",nSims,"parameter combination(s)")
	
	wd <- getwd()
	setwd(tempdir())
	.ms.setSeed()

	for ( n in 1:nSims ) {
		suppressWarnings(ms <- system2(dm.mscmd,.ms.generateCmd(dm,parameters[n,]),stdout=T))
		sumStats[n,] <- sumStatFunc(.ms.getJSFS(dm,ms))
		.dm.log(dm,"SumStats:",sumStats[n,])
	}

	setwd(wd)
	.dm.log(dm,"Finished .ms.simSumStats()")
	return(sumStats)
}

.ms.setSeed <- function(){
	cat(sample(1:65525,3),"\n",file="seedms")
	return(T)
}
