#' @importFrom R6 R6Class
block_class <- R6Class("Block", 
  private = list(border = NULL),
  public = list(
    initialize = function(border) {
      assert_that(is.matrix(border))
      assert_that(ncol(border) == 2)
      assert_that(nrow(border) >= 1)
      assert_that(all(border[ , 1] < border[ , 2]))
      private$border <- border
    },
    get_border = function() private$border,
    print = function() print(private$border),
    print_border = function(jaatha) {
      lower <- denormalize(private$border[ , 1], jaatha)
      upper <- denormalize(private$border[ , 2], jaatha)
      paste0(round(lower, 3), "-", round(upper, 3), collapse=" x ")
    },
    includes = function(point) {
      assert_that(length(point) == nrow(private$border))
      all(private$border[ , 1] - 1e-15 <= point & 
            point <= private$border[ , 2] + 1e-15)
    },
    get_middle = function() {
      m <- (private$border[ , 2] - private$border[ , 1]) / 2 + 
        private$border[ , 1]
      names(m) <- rownames(private$border)
      m
    },
    get_corners = function() {
      corners <- expand.grid(lapply(1:nrow(private$border), function(i) {
        private$border[i , , drop = FALSE]
      }), KEEP.OUT.ATTRS = FALSE)
      colnames(corners) <- rownames(private$border)
      as.matrix(corners)
    }
  )
)
