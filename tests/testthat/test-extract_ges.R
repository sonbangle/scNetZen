test_that("extract ges works", {
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  aggregate_out_dir = aggregate_out_dir = "~/Downloads/count"
  ges_list = read_ges(aggregate_out_dir)
  barcode.names = get_barcode_names_from_ges_list(ges_list)
  tasks = get_task_from_aggregation_samples(aggregation_samples_file, barcode.names = barcode.names)
  task = tasks[[1]]
  barcodes = task[["barcodes"]]
  barcode_group_name = task[["barcode_group_name"]]
  ges = ges_list[[1]]
  n_replicates = 5
  out = extract_ges(barcodes, ges,barcode_group_name,n_replicates = 5)
  print(out[["barcodes_ges_sum"]][1:10,, drop=FALSE])
  expect_equal(colnames(out[["barcodes_ges_sum"]]),barcode_group_name  )
})

test_that("extract_task_ges works", {
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  aggregate_out_dir = aggregate_out_dir = "~/Downloads/count"
  ges_list = read_ges(aggregate_out_dir)
  barcode.names = get_barcode_names_from_ges_list(ges_list)
  tasks = get_task_from_aggregation_samples(aggregation_samples_file, barcode.names = barcode.names)
  task = tasks[[1]]
  out = extract_task_ges(task, ges = ges_list[[1]])
  print(out[["barcodes_ges_sum"]][1:10,, drop=FALSE])
  expect_equal(colnames(out[["barcodes_ges_sum"]]), "mF4-2nd")
})

test_that("extract_task_with_ges_predefine works", {
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  aggregate_out_dir = aggregate_out_dir = "~/Downloads/count"
  ges_list = read_ges(aggregate_out_dir)
  barcode.names = get_barcode_names_from_ges_list(ges_list)
  tasks = get_task_from_aggregation_samples(aggregation_samples_file, barcode.names = barcode.names)
  task = tasks[[1]]
  ges = ges_list[[1]]
  out = extract_task_with_ges_predefine(task, ges_name = "ges")
  print(out[["barcodes_ges_sum"]][1:10,, drop=FALSE])
  expect_equal(colnames(out[["barcodes_ges_sum"]]), "mF4-2nd")
})
