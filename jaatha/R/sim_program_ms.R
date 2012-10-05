# --------------------------------------------------------------
# sim_prog_ms.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

possible.features  <- c("mutation","migration","split",
                        "recombination","size.change","growth")
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

  ms.file <- getTempFile("ms")
  
  .log3("Calling ms...")
  .Call("R_ms_main", opts, ms.file, PACKAGE = "jaatha")
  .log3("ms finished. Finished callMs()")
  return(ms.file)
}

# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsOptionsCommand <- function(dm) {
    nSample <- dm@sampleSizes
	cmd <- c('c(', '"ms"', ",", sum(nSample), ",", dm@nLoci , ",")
	cmd <- c(cmd,'"-I 2"', ",", nSample[1], ",", nSample[2], ",")

	for (i in 1:dim(dm@features)[1] ) {
		type <- as.character(dm@features[i,"type"])
        feat <- unlist(dm@features[i, ])
        
		if (type == "mutation") {
			if (dm@externalTheta) cmd <- c(cmd,'"-t 5"', ",")
            else 
              cmd <- c(cmd,'"-t"', ',', feat["parameter"], ',')
		}
		
		if (type == "split") {
			cmd <- c(cmd, '"-ej"', ',', feat["time.point"], ',',
                     feat["pop.sink"], ',', feat["pop.source"], ',')
        }
		
		if (type == "migration")
			cmd <- c(cmd, '"-em"', ',', feat['time.point'], ',',
                     feat['pop.sink'], ',', feat['pop.source']  , ',',
                     feat['parameter'], ',')
		
		if (type == "recombination") 
			cmd <- c(cmd, '"-r"', ',', feat['parameter'], ',', dm@seqLength, ',')

		if (type == "size.change"){
		        cmd <- c(cmd, '"-en"', ',', feat['time.point'], ',',
                         feat["pop.source"], ',', feat['parameter'], ',')
		}

		if (type == "growth"){
			cmd <- c(cmd, '"-eg"', ',' , feat["time.point"], ',',
                     feat["pop.source"], ',', feat["parameter"], ',')
		}
	}

    cmd <- c(cmd, '"")')
}

generateMsOptions <- function(dm, parameters) {
	.log3("Called .ms.generateCmd()")

    par.names <- dm.getParameters(dm)
    for (i in seq(along = par.names)){
      eval(parse(text=paste(par.names[i],"<-",parameters[i])))
    }

    fixed.pars <- dm@parameters[dm@parameters$fixed, ]
    if (nrow(fixed.pars) > 0) {
      for (i in 1:nrow(fixed.pars)){
        eval(parse(text=paste(fixed.pars$name[i],"<-",fixed.pars$lower.range[i])))
      }
    }

    cmd <- generateMsOptionsCommand(dm)
    cmd <- eval(parse(text=cmd))

	.log3("Finished .ms.generateCmd()")

    return(cmd)
}

printMsCommand <- function(dm) {
  for (i in 1:nrow(dm@parameters)) {
    eval(parse(text=paste(dm@parameters[i, "name"],
                          "<-",'\"',dm@parameters[i, "name"],'\"',sep="")))
  }

  cmd <- generateMsOptionsCommand(dm)
  cmd <- eval(parse(text=cmd))
  cmd <- paste(cmd, collapse=" ")

  return(cmd)
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
  .log3("Called msSingleSimFunc()")
  .log3("parameter:",parameters)
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

createSimProgram("ms", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=msSingleSimFunc)
