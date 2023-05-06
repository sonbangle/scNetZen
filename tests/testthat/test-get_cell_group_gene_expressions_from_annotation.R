test_that("get_cell_group_gene_expressions_from_annotation works", {
  aggregate_out_dir = "~/Downloads/count"
  annotation_file ="~/Downloads/count/Annotation.csv"
  outdir = paste0(aggregate_out_dir, "/GES")
  organism = "mouse"
  get_cell_group_gene_expressions_from_annotation(aggregate_out_dir,
                                                             annotation_file,
                                                             outdir,
                                                             organism,
                                                             method = "mpi",
                                                             close_slaves = TRUE,
                                                             n_replicates = 5,
                                                             min_n_reads_per_cell_group = 20000,
                                                             prefix = "cell_group")
  check_file= "/home/sonle/Downloads/count/GES/CELL_GROUP/CELL_GROUP_COUNTS/CONSOLIDATED_COUNTS/tpm/gene_level/consolidated_tpm_table.csv"
  expect_true(file.exists(check_file))
})
