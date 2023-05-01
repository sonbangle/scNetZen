test_that("write_cell_group_gene_expression_from_counts works", {
  extract_ges_obj = readRDS("/home/sonle/Downloads/scNetZen/tests/extract_ges_obj.RDS")
  counts = make_count_data_frame_from_extract_ges_obj(extract_ges_obj)
  cell_group_counts = counts[["counts"]]

  aggregate_out_dir = "~/Downloads/count"
  gene_names <-read_ges(aggregate_out_dir)[[2]][2]
  organism = "mouse"
  outdir = paste0(aggregate_out_dir, "/GES")
  prefix = "sample"
  min_n_reads_per_cell_group =20000
  write_cell_group_gene_expression_from_counts(cell_group_counts,
                                                          gene_names,
                                                          outdir,
                                                          prefix,
                                                          organism,
                                                          min_n_reads_per_cell_group)
  check_file =  paste0(outdir, "/",toupper(prefix), "/", prefix, "_ges_cpm.csv")
  print(check_file)
  expect_true(file.exists(check_file))
})
