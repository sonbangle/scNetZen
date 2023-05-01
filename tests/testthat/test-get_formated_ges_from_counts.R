test_that("get_formated_ges_from_count works", {

  extract_ges_obj = readRDS("/home/sonle/Downloads/scNetZen/tests/extract_ges_obj.RDS")
  counts = make_count_data_frame_from_extract_ges_obj(extract_ges_obj)
  cell_group_ges = counts[["counts"]]
  cell_group_ges_replicates = counts[["count_replicates"]]

  aggregate_out_dir = "~/Downloads/count"
  gene_names <-read_ges(aggregate_out_dir)[[2]][2]

  out = get_formated_ges_from_counts(cell_group_ges,  gene_names, min_n_reads_per_cell_group= 20000)
  for (df in out)
  {
    print(df[1:5,1:2])
  }
  print(colnames(out[["cell_group_ges"]]))


  expect_equal(sum(colnames(out[["cell_group_ges"]]) == c("SYMBOL", "mF4-2nd", "dP-2nd", "Dan_EV")), 4  )
})
