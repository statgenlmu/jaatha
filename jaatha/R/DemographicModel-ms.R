executable <-
  c("ms","hudsons_ms","/bin/ms","/usr/local/bin/ms","/usr/local/bin/hudsons_ms","~/bin/ms")
features <- c("mutation","migration","split","recombination","splitSize","presentSize")

# From phyclust package
callMs <- function(nsam = NULL, nreps = 1, opts = NULL){
  temp.file.ms <- tempfile("ms.")

  if (is.null(nsam)) stop("No sample size given")
  if (is.null(opts)) stop("No options for calling ms")
    
  nsam <- as.character(nsam)
  nreps <- as.character(nreps)
  argv <- c("ms", nsam, nreps, unlist(strsplit(opts, " ")))

  .Call("R_ms_main", argv, temp.file.ms, PACKAGE = "jaatha")
  ret <- scan(file = temp.file.ms,
              what = "character", sep = "\n", quiet = TRUE)

  unlink(temp.file.ms)
  return(ret)
}


.ms.generateCmd <- function(dm,parameters){
	.log3("Called .ms.generateCmd()")

	nSample <- dm@sampleSizes
	cmd <- c(sum(nSample),dm@nLoci)                     #repetitions and number of loci
	cmd <- c(cmd,"-I 2",nSample[1],nSample[2])          #two populations
	
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
	.log3("Generated simulation cmd:",cmd)
	.log3("Finished .ms.generateCmd()")

    return(cmd)
}

.ms.getParameter <- function(dm,parameters,type,population=NA){
	feature <- .getFeature(dm,type,population)
	if ( dim(feature)[1] != 1 ) stop("Error creating ms command")
	if ( dm@externalTheta & type == "mutation") return(0)
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
	.log3("Called .ms.getJSFS()")
	resultSize <- (dm@sampleSizes[1]+1)*(dm@sampleSizes[2]+1)
	jsfs <- ( .C("ms2jsfs",
		as.character(msOut),
		as.integer(dm@sampleSizes[1]),
		as.integer(dm@sampleSizes[2]),
		as.integer(dm@nLoci),
		res=integer(resultSize),
	        PACKAGE="jaatha")$res )
	#unlink(.ms.outfile())
	.log3("Finished .ms.getJSFS()")
	return(matrix(jsfs,dm@sampleSizes[1]+1,dm@sampleSizes[2]+1))
}


#' These are the default summary statistics for Jaatha
#' 
#' @param jsfs        The joint site frequency spectrum of two populations
#' @return        A vector with sums over different areas of the JSFS
##' @export
#'
#' @examples
#' jsfs <- matrix(rpois(26*26,5),26,26)
#' dm.defaultSumStats(jsfs)
.ms.defaultSumStats <- function(jsfs) {
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

.ms.sumStatFunc <- function(dm,jsfs,simOutput){
  return(.ms.defaultSumStats(jsfs))
}

.ms.seedFunc <- function(){
  seedfile <- paste(tempdir(),"/","seedms",sep="")
  cat(sample(1:65525,3), "\n", file=seedfile)
}

ms <- new("SimProgram","ms",executable,features,.ms.generateCmd,
          .ms.sumStatFunc,.ms.getJSFS,.ms.seedFunc)

if(!exists("dm.defaultSimProgs")) {
  dm.defaultSimProgs <- list(ms=ms)
} else {
  dm.defaultSimProgs[['ms']] <- ms
}
