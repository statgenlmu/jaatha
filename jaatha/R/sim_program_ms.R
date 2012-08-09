possible.features  <- c("mutation","migration","split",
                        "recombination","splitSize","presentSize")
possible.sum.stats <- c("jsfs")

#' Function to perform simulation using ms 
#' 
#' Copyright note: Originally taken from 'phyclust' package and modified.
#'
#' @param opts The options to pass to ms. Must either be a character or character
#' vector.
callMs <- function(opts){
  if (missing(opts)) stop("No options given!")
  opts <- unlist(strsplit(opts, " "))
  .log3("Called callMs")
  .log3("Options:", opts)

  ms.file <- tempfile("ms.")
  argv <- c("ms", opts)
  
  .log3("Calling ms...")
  .Call("R_ms_main", argv, ms.file, PACKAGE = "jaatha")
  .log3("ms finished. Finished callMs()")
  return(ms.file)
}

generateMsOptions <- function(dm, parameters) {
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

msOut2Jsfs <- function(dm, ms.out) {
  .log3("Called .ms.getJSFS()")
  jsfs <- rep(0,(dm@sampleSizes[1]+1)*(dm@sampleSizes[2]+1))
  jsfs <- matrix( .C("msFile2jsfs",
                     as.character(ms.out),
                     as.integer(dm@sampleSizes[1]),
                     as.integer(dm@sampleSizes[2]),
                     res=as.integer(jsfs),
                     PACKAGE="jaatha")$res,
                 dm@sampleSizes[1] + 1 ,
                 dm@sampleSizes[2] + 1 ,
                 byrow=T)
  .log3("Finished .ms.getJSFS()")
  return(jsfs)
}

msSingleSimFunc <- function(dm, parameters) {
  .log3("Called msSumSunc()")
  .check.dm(dm)
  checkType(parameters, "num")
  if (length(parameters) != dm.getNPar(dm)) 
    stop("Wrong number of parameters!")

  .log3("Calling ms.")
  ms.options <- generateMsOptions(dm, parameters)
  ms.out  <- callMs(ms.options)
  .log3("Simulation output in file", ms.out)
  
  .log3("Calculation jsfs...")
  jsfs  <- msOut2Jsfs(dm, ms.out) 
  unlink(ms.out)
  return(jsfs)
}

#' These are the default summary statistics for Jaatha
#' 
#' @param jsfs    The joint site frequency spectrum of two populations
#' @param dm      The corresponding demographic model. Not needed at the moment.
#' @return        A vector with sums over different areas of the JSFS
#' @export
#'
#' @examples
#' jsfs <- matrix(rpois(26*26,5),26,26)
#' ms.defaultSumStats(jsfs = jsfs)
ms.defaultSumStats <- function(dm = NULL, jsfs) {
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

createSimProgram("ms", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=msSingleSimFunc,
                 defaultSumStatFunc=ms.defaultSumStats)
