test.generateMsmsOptionsCommand <- function() {
  dm <- dm.addPositiveSelection(dm.tt, 100, 500, population=1, at.time="2") 
  opts <- generateMsmsOptionsCommand(dm)
  s <- 5
  checkTrue( "-SI" %in% eval(parse(text=opts)) )
  checkTrue( "-SAA" %in% eval(parse(text=opts)) )
  checkTrue( "-SAa" %in% eval(parse(text=opts)) )
}
