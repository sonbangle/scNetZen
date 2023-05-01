test_that("get task from aggregation sample  works", {

  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_data(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  print(tasks[[1]])
  print(paste("total number of tasks:", length(tasks)))

  expect(length(tasks), length(barcode.names))

})
