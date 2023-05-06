test_that("get_aggregation_sample_gene_expressions using multicore works", {
  aggregate_out_dir = "~/Downloads/count"
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  outdir = paste0(aggregate_out_dir, "/GES_MULTICORE")
  organism = "mouse"
  get_aggregation_sample_gene_expressions(aggregate_out_dir,
                                          aggregation_samples_file,
                                          outdir,
                                          organism,
                                          method = "multicore",
                                          close_slaves = TRUE,
                                          n_replicates = 5,
                                          min_n_reads_per_cell_group = 20000)

  prefix = "aggregation_sample"
  check_file =  paste0(outdir, "/",toupper(prefix), "/", prefix, "_ges_cpm.csv")
  print(check_file)
  expect_true(file.exists(check_file))
})



test_that("get_aggregation_sample_gene_expressions using mpi works", {
  aggregate_out_dir = "~/Downloads/count"
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  outdir = paste0(aggregate_out_dir, "/GES_MPI")
  organism = "mouse"
  get_aggregation_sample_gene_expressions(aggregate_out_dir,
                                                     aggregation_samples_file,
                                                     outdir,
                                                     organism,
                                                     method = "mpi",
                                                     close_slaves = TRUE,
                                                     n_replicates = 5,
                                                     min_n_reads_per_cell_group = 20000)

  prefix = "aggregation_sample"
  check_file =  paste0(outdir, "/",toupper(prefix), "/", prefix, "_ges_cpm.csv")
  print(check_file)
  expect_true(file.exists(check_file))
})


test_that("get_aggregation_sample_gene_expressions using single core works", {
  aggregate_out_dir = "~/Downloads/count"
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  outdir = paste0(aggregate_out_dir, "/GES_SINGLE")
  organism = "mouse"
  get_aggregation_sample_gene_expressions(aggregate_out_dir,
                                          aggregation_samples_file,
                                          outdir,
                                          organism,
                                          method = "singlecore",
                                          close_slaves = TRUE,
                                          n_replicates = 5,
                                          min_n_reads_per_cell_group = 20000)

  prefix = "aggregation_sample"
  check_file =  paste0(outdir, "/",toupper(prefix), "/", prefix, "_ges_cpm.csv")
  print(check_file)
  expect_true(file.exists(check_file))
})



