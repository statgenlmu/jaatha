# --------------------------------------------------------------
# Uses scrm to simulate demographic models
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Licence:  GPLv3 or later
# --------------------------------------------------------------

scrm_features  <- c("sample", "loci.number", "loci.length",
                    "mutation", "migration", "split",
                    "recombination", "size.change", "growth")

scrm_sum_stats <- c("jsfs", "fpc", "seg.sites", "file", "pmc")

scrm_simulate <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")

  if (length(parameters) != dm.getNPar(dm)) stop("Wrong number of parameters!")

  args <- paste(sum(dm.getSampleSize(dm)),
                dm.getLociNumber(dm),
                paste(generateMsOptions(dm, parameters), collapse = ' '))
  
  sum_stats <- scrm(args)
  generateSumStats(sum_stats, 2, parameters, dm)
}

scrm_finalize <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}

createSimProgram("scrm", scrm_features, scrm_sum_stats,
                  scrm_simulate, scrm_finalize, printMsCommand, 100)
rm(scrm_features, scrm_sum_stats, scrm_simulate, scrm_finalize)