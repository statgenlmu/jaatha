possible.features  <- c("mutation","migration","split",
                        "recombination","splitSize","presentSize",
                        "time.point")
possible.sum.stats <- c("jsfs")

#' Function to perform simulation using ms 
#' 
#' @param opts The options to pass to ms. Must either be a character or character
#' vector.
callMs <- function(opts){
  if (missing(opts)) stop("No options given!")
  opts <- unlist(strsplit(opts, " "))
  .log3("Called callMs")
  .log3("Options:", opts)

  ms.file <- tempfile("ms.")
  
  .log3("Calling ms...")
  .Call("R_ms_main", opts, ms.file, PACKAGE = "jaatha")
  .log3("ms finished. Finished callMs()")
  return(ms.file)
}

generateMsOptions <- function(dm, parameters) {
	.log3("Called .ms.generateCmd()")

    par.names <- unique(dm.getParameters(dm, include.fixed=F))
    for (i in seq(along = par.names)){
      eval(parse(text=paste(par.names[i],"<-",parameters[i])))
    }

    fixed.features <- dm@features[dm@features$fixed & !is.na(dm@features$parameter), ]
    if (nrow(fixed.features) > 0) {
      for (i in 1:nrow(fixed.features)){
        eval(parse(text=paste(fixed.features$parameter[i],"<-",fixed.features$lower.range[i])))
      }
    }

    print(ls())

    nSample <- dm@sampleSizes
	cmd <- c('c(', '"ms"', ",", sum(nSample), ",", dm@nLoci , ",")
	cmd <- c(cmd,'"-I 2"', ",", nSample[1], ",", nSample[2], ",")

    par.names <- dm@features$parameter

	for (i in 1:dim(dm@features)[1] ) {
		type <- as.character(dm@features[i,"type"])
        feat <- unlist(dm@features[i, ])
        
        if (type == "time.point") next

		if (type == "mutation") {
			if (dm@externalTheta) cmd <- c(cmd,'"-t 5"', ",") 
			else cmd <- c(cmd,'"-t"', ',', par.names[i], ',')
		}
		
		if (type == "split") {
			cmd <- c(cmd, '"-ej"', ',', feat["time.point"], ',',
                     feat["pop.sink"], ',', feat["pop.source"], ',')
        }
		
		if (type == "migration")
			cmd <- c(cmd, '"-m 1 2"', ',', feat['parameter'], ',',
                          '"-m 2 1"', ',', feat['parameter'], ',')
		
		if (type == "recombination") 
			cmd <- c(cmd,"-r",param,dm@seqLength)

		if (type == "splitSize"){
		        cmd <- c(cmd,"-g",feat["pop.source"],
                         log(.ms.getPresentSize(dm,parameters,feat["pop.source"])/param)/tau)  
			cmd <- c(cmd,"-eN",tau,1+param)               
		}

		if (type == "presentSize"){
			cmd <- c(cmd,"-n",feat["pop.source"],param)
		}
	}

    cmd <- c(cmd, '"")')
    
    cmd <-  eval(parse(text=cmd))

	.log3("Finished .ms.generateCmd()")

    return(cmd)
}

.ms.getParameter <- function(dm, parameters, row.nr, par.name=NA) {
    # Calling with row.nr instead of par name
    if (is.na(par.name)) par.name <- dm@features$parameter[row.nr]
    # Calling with row.nr and the par is fixed
    if (is.na(par.name)) return(dm@features$lower.range[row.nr])
    #par.name <- unlist(par.name)
    par.mask <- dm.getParameters(dm) == par.name
    par.nr <- seq(along = par.mask)[par.mask]
    return(parameters[par.nr])
}

.ms.getPresentSize <- function(dm, parameters, pop.source){
	feature <- getFeature(dm, "presentSize", pop.source)
	if ( dim(feature)[1] != 1 ) return(1)
	if ( feature$lower.range[1] == feature$upper.range[1] )
		return(feature$lower.range[1]) 
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
  checkType(dm, "dm")
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
