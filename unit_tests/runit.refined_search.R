test.calcBorders <- function() {
  point <- c(.5, .2)
  borders <- calcBorders(point, 0.1)
  checkEquals( borders[1, ], c(lower=0.4, upper=0.6) )
  checkEquals( borders[2, ], c(lower=0.1, upper=0.3) )

  borders <- calcBorders(point, 0.3)
  checkEquals( borders[1, ], c(lower=0.2, upper=0.8) )
  checkEquals( borders[2, ], c(lower=0.0, upper=0.6) )

  borders <- calcBorders(point+.4, 0.3)
  checkEquals( borders[1, ], c(lower=0.4, upper=1) )
  checkEquals( borders[2, ], c(lower=0.3, upper=.9) )

  checkException(calcBorder(c(-.5, .4), .1))
  checkException(calcBorder(c(.5, 1), .1))
  checkException(calcBorder(c(.5, 1), -.1))
}
