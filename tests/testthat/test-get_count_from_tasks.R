test_that("get_count_from_task works", {


  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  load_ges_slave(aggregate_out_dir)


  counts =get_count_from_tasks(tasks,
                        method= "mpi",
                        close_slaves = TRUE,
                        n_replicates = 5)
  cell_group_counts = counts[[1]]
  cell_group_counts_replicates =  counts[[2]]
  print(cell_group_counts[1:5,])
  print(cell_group_counts_replicates[1:5,])
  expect_equal(ncol(cell_group_counts), length(tasks) +1)
  expect_equal(ncol(cell_group_counts_replicates), length(tasks) * 5 + 1)

  expect_equal(2 * 2, 4)
})
