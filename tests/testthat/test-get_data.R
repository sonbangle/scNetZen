test_that("get data work", {
  aggregate_out_dir = "~/Downloads/count"
  out = get_data(aggregate_out_dir = aggregate_out_dir)
  print(out[1:5,])
  expect_equal(2 * 2, 4)
})
