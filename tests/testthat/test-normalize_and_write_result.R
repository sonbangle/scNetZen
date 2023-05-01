test_that("normalize_and_write_result works", {

  extract_ges_obj = readRDS("/home/sonle/Downloads/scNetZen/tests/extract_ges_obj.RDS")
  counts = make_count_data_frame_from_extract_ges_obj(extract_ges_obj)
  cell_group_ges = counts[["counts"]]

  aggregate_out_dir = "~/Downloads/count"
  gene_names <-read_ges(aggregate_out_dir)[[2]][2]

  out = get_formated_ges_from_counts(cell_group_ges,  gene_names, min_n_reads_per_cell_group= 20000)
  outdir = paste0(aggregate_out_dir, "/GES")
  prefix = "sample"
  normalize_and_write_result(out , outdir = outdir , prefix  = prefix, organism = "mouse")
  check_file =  paste0(outdir, "/", prefix, "_ges_cpm.csv")
  print(check_file)
  expect_true(file.exists(check_file))
})
