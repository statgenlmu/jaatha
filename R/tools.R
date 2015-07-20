# --- Global Variables ---------------------------------------------------------
if (!exists("jaatha_env")) jaatha_env <- new.env()

set_jaatha_var <- function(name, value) {
  jaatha_env[[name]] <- value
}

get_jaatha_var <- function(name) {
  return(get(name, envir = jaatha_env))
}

is_jaatha_var <- function(name) {
  exists(name, envir = jaatha_env)
}


# --- Argument checking --------------------------------------------------------
is_single_numeric <- function(value) is.numeric(value) && length(value) == 1
is_single_logical <- function(value) is.logical(value) && length(value) == 1

# --- Printing -----------------------------------------------------------------
.print <- function(...) { cat(..., "\n") }


# --- Seeds --------------------------------------------------------------------
sample_seed <- function(n = 1) sample.int(2 ^ 20, n)
