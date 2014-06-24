test.logfac <- function() {
  setJaathaVariable("logfacs", c(0))
  checkTrue(sum(abs(calcLogFactorial(0:10)-log(factorial(0:10)))) < 1e-10)
  checkTrue(sum(abs(calcLogFactorial(0:20)-log(factorial(0:20)))) < 1e-10)
  checkTrue(sum(abs(calcLogFactorial(25:30)-log(factorial(25:30)))) < 1e-10)
  checkTrue(sum(abs(calcLogFactorial(0:10)-log(factorial(0:10)))) < 1e-10)
}
