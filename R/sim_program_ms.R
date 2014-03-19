# --------------------------------------------------------------
# Translates an demographic model to an ms command and 
# executes the simulation.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Licence:  GPLv3 or later
# --------------------------------------------------------------

possible.features  <- c("sample", "loci.number", "loci.length",
                        "mutation", "migration", "split",
                        "recombination", "size.change", "growth")
possible.sum.stats <- c("jsfs", "4pc", "tree", "seg.sites", "file")

#' Function to perform simulation using ms 
#' 
#' @param opts The options to pass to ms. Must either be a character or character
#' vector.
#' @param dm The demographic model we are using
#' @return The file containing the output of ms
callMs <- function(opts, dm){
  if (missing(opts)) stop("No options given!")
  opts <- unlist(strsplit(opts, " "))

  ms.file <- getTempFile("ms")

  ms(sum(dm.getSampleSize(dm)), dm.getLociNumber(dm), opts, ms.file)
  return(ms.file)
}

# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsOptionsCommand <- function(dm) {
  nSample <- dm.getSampleSize(dm)
  cmd <- c('c(')
  cmd <- c(cmd,'"-I"', ",", length(nSample), ',', 
           paste(nSample, collapse=","), ',')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if (type == "mutation") {
      cmd <- c(cmd,'"-t"', ',', feat["parameter"], ',')
    }

    else if (type == "split") {
      cmd <- c(cmd, '"-ej"', ',', feat["time.point"], ',',
               feat["pop.sink"], ',', feat["pop.source"], ',')
    }

    else if (type == "migration")
      cmd <- c(cmd, '"-em"', ',', feat['time.point'], ',',
               feat['pop.sink'], ',', feat['pop.source']  , ',',
               feat['parameter'], ',')

    else if (type == "recombination") 
      cmd <- c(cmd, '"-r"', ',', feat['parameter'], ',', dm.getLociLength(dm), ',')

    else if (type == "size.change"){
      cmd <- c(cmd, '"-en"', ',', feat['time.point'], ',',
               feat["pop.source"], ',', feat['parameter'], ',')
    }

    else if (type == "growth"){
      cmd <- c(cmd, '"-eg"', ',' , feat["time.point"], ',',
               feat["pop.source"], ',', feat["parameter"], ',')
      }

    else if (type %in% c("sample", "loci.number", "loci.length")) {}
    else stop("Unknown feature:", type)
  }

  if ('trees' %in% dm@sum.stats) cmd <- c(cmd, '"-T",')
  cmd <- c(cmd, '" ")')
}

generateMsOptions <- function(dm, parameters) {
  .log3("Called .ms.generateCmd()")
  ms.tmp <- new.env()

  par.names <- dm.getParameters(dm)
  for (i in seq(along = par.names)){
    ms.tmp[[ par.names[i] ]] <- parameters[i]
  }

  fixed.pars <- dm@parameters[dm@parameters$fixed, ]
  if (nrow(fixed.pars) > 0) {
    for (i in 1:nrow(fixed.pars)){
      ms.tmp[[ fixed.pars$name[i] ]] <- fixed.pars$lower.range[i]
    }
  }

  if ( !is.null( dm@options[['ms.cmd']] ) )
    cmd <- dm@options[['ms.cmd']]
  else
    cmd <- generateMsOptionsCommand(dm)
  cmd <- eval(parse(text=cmd), envir=ms.tmp)

  .log3("Finished .ms.generateCmd()")

  return(cmd)
}

printMsCommand <- function(dm) {
  cmd <- generateMsOptionsCommand(dm)

  cmd <- cmd[cmd != ","]
  cmd <- cmd[-c(1, length(cmd))]

  cmd <- paste(cmd, collapse=" ")

  cmd <- gsub(",", " ", cmd)
  cmd <- gsub('\"', "", cmd)
  cmd <- gsub('"', " ", cmd)

  return(cmd)
}

msOut2Jsfs <- function(dm, ms.out) {
  .log3("Called .ms.getJSFS()")
  sample.size <- dm.getSampleSize(dm)
  jsfs <- matrix(.Call("msFile2jsfs", ms.out,sample.size[1], 
                       sample.size[2]),
                 sample.size[1] + 1 ,
                 sample.size[2] + 1,
                 byrow=T)
  .log3("Finished .ms.getJSFS()")
  return(jsfs)
}

msSingleSimFunc <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")

  if (length(parameters) != dm.getNPar(dm)) stop("Wrong number of parameters!")

  ms.options <- generateMsOptions(dm, parameters)
  sim.time <- system.time(ms.out <- callMs(ms.options, dm))

  sum.stats <- list(pars=parameters)

  if ("jsfs" %in% dm@sum.stats) {
    sum.stats[['jsfs']] <- msOut2Jsfs(dm, ms.out)
  }

  if ("file" %in% dm@sum.stats) {
    sum.stats[['file']] <- ms.out
  }

  if (any(c('seg.sites', 'tree', '4pc') %in% dm@sum.stats)) {
    output <- scan(ms.out, character(), sep="\n", quiet=TRUE)

    if ("seg.sites" %in% dm@sum.stats) {
      sum.stats[['seg.sites']] <- readSegSitesFromOutput(output,
                                                         dm.getSampleSize(dm))
    }

    if ("4pc" %in% dm@sum.stats) {
      if (!is.null(sum.stats$seg.sites)) seg.sites <- sum.stats$seg.sites
      else seg.sites <- readSegSitesFromOutput(output, dm.getSampleSize(dm))

      sum.stats[['4pc']] <- calcFpcSumStat(seg.sites, dm)
    }
  }

  if (!'file' %in% dm@sum.stats) unlink(ms.out)
  return(sum.stats)
}

finalizeMs <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}

readSegSitesFromOutput <- function(output, pop.sizes) {
  seg.sites.begin <- which(grepl('^segsites: [0-9]+$', output))

  lapply(seg.sites.begin, function(begin) {
    positions <- strsplit(output[begin+1], ' ')[[1]][-1]
    positions <- positions[positions != ""]
    if (length(positions) == 0) {
      return(matrix(0, sum(pop.sizes), 0))
    }
    stopifnot( all(!is.na(positions)) )
    seg.sites.char <- output[1:sum(pop.sizes)+begin+1]
    seg.sites <- matrix(as.integer(unlist(strsplit(seg.sites.char, split= ''))),
                        length(seg.sites.char), byrow=TRUE)
    colnames(seg.sites) <- positions
    seg.sites
  })
}

calcFpcSumStat <- function(seg.sites, dm) {
  breaks.near <- dm@options[['4pc.breaks.near']]
  breaks.far  <- dm@options[['4pc.breaks.far']]
  breaks.theta <- dm@options[['4pc.breaks.theta']]

  fpc <- array(0, 
               list(length(breaks.near), length(breaks.far), breaks.theta),
               list(c(2:length(breaks.near)-1,'NA'),
                    c(2:length(breaks.far)-1,'NA'),
                    1:breaks.theta)) 

  loci.class <- sapply(seg.sites, calcPercentFpcViolations)
  loci.class['near',] <- cut(loci.class['near',], breaks.near, include.lowest=TRUE, labels=FALSE)
  loci.class['far',] <- cut(loci.class['far',], breaks.far, include.lowest=TRUE, labels=FALSE)
  loci.class['theta',] <- cut(loci.class['theta',], breaks.theta, include.lowest=TRUE, labels=FALSE)

  for (j in 1:ncol(loci.class)) {
    class.near <- ifelse(is.na(loci.class['near', j]), 'NA', loci.class['near', j])
    class.far <- ifelse(is.na(loci.class['far', j]), 'NA', loci.class['far', j])
    class.theta <- ifelse(is.na(loci.class['theta', j]), 'NA', loci.class['theta', j])
    fpc[class.near, class.far, class.theta] <- fpc[class.near, class.far, class.theta] + 1
  }

  return(fpc)
}

calcPercentFpcViolations <- function(snp.matrix) {
  snp.matrix <- snp.matrix[, colSums(snp.matrix)>1, drop=FALSE]
  if (ncol(snp.matrix) <= 1) return(c(near=NaN, far=NaN, theta=0))
  snp.state <- apply(combn(1:ncol(snp.matrix), 2), 2, violatesFpc, snp.matrix)
  return(c(near=sum(snp.state[2, snp.state[1, ]])/sum(snp.state[1, ]),
           far=sum(snp.state[2, !snp.state[1, ]])/sum(!snp.state[1, ]),
           theta=ncol(snp.matrix)/sum(1/1:(nrow(snp.matrix)-1)) ))
}

violatesFpc <- function(sites, snp.matrix, near=.1) {
  is.near <- diff(as.numeric(colnames(snp.matrix)[sites])) < near
  status <- snp.matrix[ ,sites[1]] * 2 + snp.matrix[ ,sites[2]] 
  if (all(0:3 %in% status)) return(c(near=is.near, violates=TRUE))
  return(c(near=is.near, violates=FALSE))
}

calcFpcBreaks <- function(dm, seg.sites, number=5) {
  props <- seq(0, 1, length.out = number + 2)[-c(1, number+2)]
  fpc.percent <- t(sapply(seg.sites, calcPercentFpcViolations))
  dm@options[['4pc.breaks.near']] <- unique(c(0, quantile(fpc.percent[ ,'near'], props, na.rm=TRUE), 1))
  dm@options[['4pc.breaks.far']] <- unique(c(0, quantile(fpc.percent[ ,'far'], props, na.rm=TRUE), 1))
  dm@options[['4pc.breaks.theta']] <- number
  dm
}

createSimProgram("ms", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=msSingleSimFunc,
                 finalizationFunc=finalizeMs)
