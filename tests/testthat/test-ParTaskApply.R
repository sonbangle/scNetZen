test_that("TaskApply", {

  start_mpi()
  tasks = c(1:5)
  get_double = function (x)
  {
    print("running function")
    out = x *2
    print(paste("in:", x, "out:", out))
    return(out)
    }
  out = TaskApply(tasks = tasks, get_double)
  print(out)
  mpi.close.Rslaves(dellog = FALSE)
  mpi.finalize()

  expect_equal(unlist(out), tasks *2)
})
