runit.generateMsmsOptionsCommand <- function() {
  dm <- dm.addPositiveSelection(dm.tt, 1, 2, population=1, at.time="2") 
  opts <- generateMsmsOptionsCommand(dm)
  checkTrue( "-SI" %in% eval(parse(text=opts)) )
}
