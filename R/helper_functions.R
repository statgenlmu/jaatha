# --- Global Variables ---------------------------------------------------------
if (!exists(".jaatha")) jaatha_env <- new.env()

setJaathaVariable <- function(name, value) {
  jaatha_env[[name]] <- value
}

getJaathaVariable <- function(name) {
  return(get(name, envir = jaatha_env))
}

isJaathaVariable <- function(name) {
  exists(name, envir = jaatha_env)
}


# --- Argument checking --------------------------------------------------------
is_single_numeric <- function(value) is.numeric(value) && length(value) == 1


# --- Printing -----------------------------------------------------------------
.print <- function(...) { cat(..., "\n") }


# --- Seeds --------------------------------------------------------------------
sampleSeed <- function(n = 1) sample.int(2 ^ 20, n)
