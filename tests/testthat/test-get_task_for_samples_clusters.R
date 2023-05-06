test_that("get_task_for_samples_clusters works", {
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  clusters_file = "~/Downloads/count/analysis/clustering/graphclust/clusters.csv"
  res = get_task_for_samples_clusters(aggregation_samples_file , clusters_file,  barcode.names)
  print(res)
  tasks = res[["tasks"]]
  sample_cluster_matrix = res[["sample_cluster_matrix"]]
  expect_equal(sum(names(tasks[[1]])  == c("barcodes", "barcode_group_name")) , 2)
})
