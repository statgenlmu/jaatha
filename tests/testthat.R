if (require("testthat")) {
  test_check("jaatha")
} else {
  warning("testthat not available. Skipping unittests!")
}
