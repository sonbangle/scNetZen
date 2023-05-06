test_that("get_sample_cluster_gene_expressions works", {

  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  clusters_file = "~/Downloads/count/analysis/clustering/graphclust/clusters.csv"
  organism = "mouse"
  outdir =  "~/Downloads/count/GES"
  get_sample_cluster_gene_expressions(aggregate_out_dir,
                                                 aggregation_samples_file,
                                                 clusters_file,
                                                 outdir,
                                                 organism,
                                                 method = "mpi",
                                                 close_slaves = TRUE,
                                                 n_replicates = 5,
                                                 min_n_reads_per_cell_group = 20000)


  check_file = paste0(outdir, "/SAMPLE_CLUSTER/sample_cluster_ges.csv")
  expect_true(file.exists(check_file))
})
