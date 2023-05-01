test_that("get task from aggregation sample  works", {
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  aggregate_out_dir = aggregate_out_dir = "~/Downloads/count"
  ges_list = read_ges(aggregate_out_dir)
  barcode.names = get_barcode_names_from_ges_list(ges_list)
  print(barcode.names[1:5,])
  tasks = get_task_from_aggregation_samples(aggregation_samples_file, barcode.names = barcode.names)
  print(names(tasks[[1]]))
  expect_true(is.list(tasks))
  expect_equal(names(tasks[[1]]), c("barcodes", "barcode_group_name", "n_replicates") )
})
