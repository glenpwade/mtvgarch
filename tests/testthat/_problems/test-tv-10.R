# Extracted from test-tv.R:10

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "MTVGARCH", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
Tobs <- 100
st <- (1:Tobs) / Tobs
obj_delta0 <- tv(st, shape = tvshape$delta0only)
expect_s4_class(obj_delta0, "tv_class")
expect_equal(obj_delta0$shape, tvshape$delta0only)
expect_equal(obj@nr.pars, as.integer(1))
