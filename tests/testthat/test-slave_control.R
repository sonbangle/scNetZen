test_that("slave control can receive simple function", {
  start_mpi()
  fun = function(x) x*2
  mpi.bcast.Robj2slave(fun)
  mpi.bcast.cmd(slave_control(fun))
  expect_equal(2 * 2, 4)
  mpi.close.Rslaves()
})


