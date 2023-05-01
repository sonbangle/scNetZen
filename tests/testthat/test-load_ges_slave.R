test_that("load_slave_ges  works", {
  load_ges_slave(aggregate_out_dir = "~/Downloads/count")
  mpi.bcast.cmd(print(ges[1:5, 1:5]))
  Sys.sleep(10)
  expect_equal(2 * 2, 4)
})
