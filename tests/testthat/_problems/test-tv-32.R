# Extracted from test-tv.R:32

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "MTVGARCH", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
st <- (1:10) / 10
obj <- tv(st, shape = tvshape$single)
expect_equal(obj@nr.pars, as.integer(1))
