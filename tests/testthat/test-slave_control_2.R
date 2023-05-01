
test_that("slave control can receive function", {
  aggregate_out_dir = "~/Downloads/count"
  #barcode.names = get_data(aggregate_out_dir = aggregate_out_dir)
  #aggregation_samples_file = "~/Downloads/count/aggregation.csv"

  #load_ges_slave(aggregate_out_dir, n_replicates = 5)
  start_mpi()
  fun= extract_task_with_ges_predefine
  mpi.bcast.Robj2slave(fun)
  mpi.bcast.cmd(slave_control(fun))
  Sys.sleep(30)
  mpi.close.Rslaves()
  expect_equal(2 * 2, 4)

})

