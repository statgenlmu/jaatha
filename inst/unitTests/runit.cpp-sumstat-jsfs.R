test.addSegSitesToJsfs <- function() {
  seg.sites <- matrix(c(1,0,0,0, 1,1,0,1, 1,0,0,1), 4, 3)  
  jsfs <- matrix(0, 3, 3)
  jsfs.new <- addSegSitesToJsfs(seg.sites, c(2,2), jsfs)
  checkTrue( is.matrix(jsfs.new) )
  checkEquals( c(3,3), dim(jsfs.new) )
  checkEquals( 0, sum(jsfs) )
  checkEquals( 3, sum(jsfs.new) )
  checkEquals( 1, jsfs.new[2,1] )
  checkEquals( 1, jsfs.new[3,2] )
  checkEquals( 1, jsfs.new[2,2] )
}
