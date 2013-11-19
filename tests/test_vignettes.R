jaatha.howto <- system.file("doc/jaatha_howto.pdf", package="jaatha", mustWork=TRUE)
csi.howto <- system.file("doc/custom_simulator_howto.pdf", package="jaatha", mustWork=TRUE)

file.info(jaatha.howto)[1, "size"] > 75000 || stop("Seems we are shipping a dummy vignette")
file.info(csi.howto)[1, "size"] > 75000 || stop("Seems we are shipping a dummy vignette")
