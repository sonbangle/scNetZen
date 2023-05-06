test_that("extract ges  by lapply works", {

  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  #load_ges_slave(aggregate_out_dir, n_replicates = 5)
  ges<-read_ges(aggregate_out_dir)[[1]]
  out = extract_ges_parallel(tasks, method = "singlecore")
  res = out[[1]]
  print(res[["barcodes_ges_sum"]][1:5,, drop=FALSE])
  print(res[["barcodes_ges_sum_replicates"]][1:5,1:5])
  expect_equal(sum(names(res) == c("barcodes_ges_sum", "barcodes_ges_sum_replicates", "is_null_result")),3)
})

test_that("extract ges  by mc lapply works", {

  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  ges<-read_ges(aggregate_out_dir)[[1]]
  out = extract_ges_parallel(tasks, method = "multicore")
  res = out[[1]]
  print(res[["barcodes_ges_sum"]][1:5,, drop=FALSE])
  print(res[["barcodes_ges_sum_replicates"]][1:5,1:5])
  expect_equal(sum(names(res) == c("barcodes_ges_sum", "barcodes_ges_sum_replicates", "is_null_result")),3)
})

test_that("extract ges mpi works", {

  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  load_ges_slave(aggregate_out_dir)
  out = extract_ges_parallel(tasks, method = "mpi")
  res = out[[1]]
  print(res[["barcodes_ges_sum"]][1:5,, drop=FALSE])
  print(res[["barcodes_ges_sum_replicates"]][1:5,1:5])
  Sys.sleep(30)
  expect_equal(sum(names(res) == c("barcodes_ges_sum", "barcodes_ges_sum_replicates", "is_null_result")),3)
})

test_that("extract ges  by lapply works for n_replicates =2 ", {

  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  #load_ges_slave(aggregate_out_dir, n_replicates = 5)
  ges<-read_ges(aggregate_out_dir)[[1]]
  out = extract_ges_parallel(tasks, method = "singlecore", n_replicates = 2)
  res = out[[1]]
  print(res[["barcodes_ges_sum"]][1:5,, drop=FALSE])
  print(res[["barcodes_ges_sum_replicates"]][1:5,])
  expect_equal(sum(names(res) == c("barcodes_ges_sum", "barcodes_ges_sum_replicates", "is_null_result")),3)
  expect_equal(ncol(res[["barcodes_ges_sum_replicates"]]), 2)
})

test_that("extract ges  by mc lapply works for n_replicates =2  ", {

  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  ges<-read_ges(aggregate_out_dir)[[1]]
  out = extract_ges_parallel(tasks, method = "multicore", n_replicates = 2)
  res = out[[1]]
  print(res[["barcodes_ges_sum"]][1:5,, drop=FALSE])
  print(res[["barcodes_ges_sum_replicates"]][1:5,])
  expect_equal(ncol(res[["barcodes_ges_sum_replicates"]]), 2)
  expect_equal(sum(names(res) == c("barcodes_ges_sum", "barcodes_ges_sum_replicates", "is_null_result")),3)
})

test_that("extract ges mpi works for n_replicates =2 ", {

  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  load_ges_slave(aggregate_out_dir)
  out = extract_ges_parallel(tasks, method = "mpi",  n_replicates = 2)
  res = out[[1]]
  print(res[["barcodes_ges_sum"]][1:5,, drop=FALSE])
  print(res[["barcodes_ges_sum_replicates"]][1:5,])
  Sys.sleep(30)
  expect_equal(ncol(res[["barcodes_ges_sum_replicates"]]), 2)
  expect_equal(sum(names(res) == c("barcodes_ges_sum", "barcodes_ges_sum_replicates", "is_null_result")),3)
})


