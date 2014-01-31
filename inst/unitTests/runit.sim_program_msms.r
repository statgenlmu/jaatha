# ---------------------------------------------------
# Unit tests for functions in sim_prog_msms.r
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Licence:  GPLv3 or later
# ---------------------------------------------------

example.msms.output <- 
  c('ms -oAFS jAFS 5 2 -r 10 100 -t 5 -I 2 3 2 1 -seed 446666  [3.2rc Build:147]',
    '0x6d0ca',
    '',
    '//',
    'segsites: 25',
    'positions: 0.00198 0.12622 0.22100 0.28564 0.29314 0.33176 0.39234 0.40854 0.46026 0.46430 0.47525 0.50528 0.57724 0.60490 0.60827 0.61665 0.65229 0.66336 0.67762 0.74818 0.86222 0.87351 0.96671 0.97816 0.98457 ',
    '0001100100010011001101000',
    '0101100110000010001101000',
    '0010011000001000010010111',
    '1001100100110011001101000',
    '0001100011000100100010111',
    '',
    'AFS: 12 6 5 2 ',
    'jAFS 0 vrs 1',
    '0 5 0',
    '7 6 0',
    '0 5 2',
    '0 0 0',
    '',
    '//',
    'segsites: 24',
    'positions: 0.02793 0.05862 0.09954 0.12986 0.15803 0.17430 0.29267 0.31987 0.41199 0.47040 0.47904 0.50675 0.55857 0.57181 0.58251 0.59448 0.60219 0.65597 0.68429 0.71842 0.72497 0.72744 0.78606 0.91129 ',
    '100100010000000000011010',
    '101011001000011000000000',
    '101001000000011000000010',
    '010000100011111101000001',
    '000000100100100010100111',
    '',
    'AFS: 16 4 4 0 ',
    'jAFS 0 vrs 1',
    '0 10 2',
    '6 0 0',
    '2 3 0',
    '1 0 0',
    '',
    'Summary AFS:28 10 9 2',
    'Summary jAFS',
    'jAFS 0 vrs 1',
    '0 15 2',
    '13 6 0',
    '2 8 2',
    '1 0 0',
    '')

test.readJsfsFromOutput <- function() {
  jsfs <- readJsfsFromOutput(example.msms.output)
  checkEquals(jsfs, matrix(c(0,13,2,1,15,6,8,0,2,0,2,0), 4, 3))
  checkException(readJsfsFromOutput(example.msms.output[-(38:42)]))
}

test.readSegSitesFromOutput <- function() {
  seg.sites <- readSegSitesFromOutput(example.msms.output, c(3,2))
  checkTrue(is.list(seg.sites))
  checkEquals(length(seg.sites), 2)
  checkTrue(all(seg.sites[[1]][1,1:5] == c(0, 0, 0, 1, 1)))
  checkTrue(all(seg.sites[[2]][1,1:5] == c(1, 0, 0, 1, 0)))
  checkTrue(all(colnames(seg.sites[[2]])[1:2] == c(0.02793, 0.05862)))
}

test.generateFourPointStat <- function() {
  seg.sites <- readSegSitesFromOutput(example.msms.output, c(3,2))
  fps <- generateFourPointStat(seg.sites, c(3,2))
  checkTrue(is.vector(fps))
  checkTrue(is.numeric(fps))
  checkTrue(length(fps) == 6)
  checkTrue(all(fps[c(1,2,4,5)] == 0))
  checkEquals(fps[3], 3)
}

test.msms <- function() {
  jar.path <- '~/bin/msms.jar'
  if (!file.exists(jar.path)) {
    warning('Can not test msms: jar not found')
    return()
  }

  ms.args <- '5 1 -r 10 100 -t 5 -I 2 3 2 1'
  msms.args <- ''
  pop.sizes <- c(3,2)
  
  set.seed(17)
  sum.stats <- msms(jar.path, ms.args, msms.args, pop.sizes) 
  checkTrue( is.list(sum.stats) )
  checkTrue( length(sum.stats) == 1 )
  checkTrue( !is.null(sum.stats$raw) )

  possible.sum.stats <- c('raw', 'jsfs', 'seg.sites', 'four.point.condition')
  sum.stats <- msms(jar.path, ms.args, msms.args, pop.sizes, possible.sum.stats) 
  checkTrue( is.list(sum.stats) )
  checkEquals( length(sum.stats), length(possible.sum.stats) )

  checkTrue( is.vector(sum.stats$raw) )
  checkTrue( is.matrix(sum.stats$jsfs) )
  checkEquals( dim(sum.stats$jsfs), 4:3 )
  checkTrue( sum(sum.stats$jsfs) > 0 )

  checkTrue( is.list(sum.stats$seg.sites) )
  checkTrue( is.vector(sum.stats$four.point.condition) )

  checkException( msms(jar.path, ms.args, msms.args, pop.sizes, 'blub') ) 
}

test.msms.singleSimFunc <- function() {
  jar.path <- '~/bin/msms.jar'
  if (!file.exists(jar.path)) {
    warning('Can not test msms: jar not found')
    return()
  }

  checkForMsms()
  sum.stats <- msmsSingleSimFunc(dm.sel, c(1, 1.5, 1500, 5))
  print(sum.stats)
  checkTrue( is.list(sum.stats) )
  checkTrue( is.matrix(sum.stats$jsfs) )
}
