test_that("get_cell_group_gene_expression_from_tasks works", {
  aggregate_out_dir = "~/Downloads/count"
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  aggregation_samples_file = "~/Downloads/count/aggregation.csv"
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  outdir = paste0(aggregate_out_dir, "/GES")
  prefix = "Aggregation_sample"
  organism = "mouse"
  out = get_cell_group_gene_expression_from_tasks(tasks,
                                            aggregate_out_dir,
                                            outdir,
                                            prefix,
                                            organism,
                                            method = "mpi",
                                            close_slaves = TRUE,
                                            n_replicates = 5,
                                            min_n_reads_per_cell_group = 20000)

  check_file =  paste0(outdir, "/",toupper(prefix), "/", prefix, "_ges_cpm.csv")
  print(check_file)
  expect_true(file.exists(check_file))
})
