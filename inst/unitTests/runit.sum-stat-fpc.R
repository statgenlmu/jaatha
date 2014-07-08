test.calcFpcBreaks = function() {
  dm = calcFpcBreaks(dm.fpc, seg.sites)
  checkTrue( !is.null(dm@options[['fpc.breaks.near']]) )
  checkTrue( !is.null(dm@options[['fpc.breaks.far']]) )
  
  dm = calcFpcBreaks(dm.fpc, seg.sites, group=1)
  checkTrue( !is.null(dm@options[['group.1']][['fpc.breaks.near']]) )
  checkTrue( !is.null(dm@options[['group.1']][['fpc.breaks.far']]) )
  
  dm = calcFpcBreaks(dm, seg.sites, group=2)  
  checkTrue( !is.null(dm@options[['group.1']][['fpc.breaks.near']]) )
  checkTrue( !is.null(dm@options[['group.1']][['fpc.breaks.far']]) )
  checkTrue( !is.null(dm@options[['group.2']][['fpc.breaks.near']]) )
  checkTrue( !is.null(dm@options[['group.2']][['fpc.breaks.far']]) )
}