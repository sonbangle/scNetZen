test_that("get_task_from_annotation works", {
  tasks = get_task_from_annotation(annotation_file ="~/Downloads/count/Annotation.csv")
  print(tasks)
  task = tasks[[1]]
  barcode_group_name = task[["barcode_group_name"]]
  print(task)
  expect_true(is.character(barcode_group_name))
})
