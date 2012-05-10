#Loads the package binaries on package load
.First.lib <- function(lib, pkg)
{
  library.dynam("jaatha", pkg, lib)
}

#Unloads the package binaries on package unload
.Last.lib <- function(libpath)
{
  library.dynam.unload("jaatha", libpath)
}
