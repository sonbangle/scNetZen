test_that("sonic apply mpi works", {
  x = list(1,2,3,4)
  fun = function(x)
    {
    return(x*2)
  }

  out = sonicApply(x, fun, method = "mpi")
  expect_equal(out, list(2,4,6,8))

})

test_that("sonic apply multicore works", {
  x = list(1,2,3,4)
  fun = function(x) x*2
  out = sonicApply(x, fun, method = "multicore")
  expect_equal(out, list(2,4,6,8))
})

test_that("sonic apply single works", {
  x = list(1,2,3,4)
  fun = function(x) x*2
  out = sonicApply(x, fun, method = "singlecore")
  expect_equal(out, list(2,4,6,8))
})
