test_that("TaskApply works with simple function function with optional arguments", {
  tasks = list(1,2,3,4)
  print(length(tasks))
  start_mpi()
  fun = function(x, a)
  {
    print(x)
    return(x*a)
  }
  a = 3
  out = TaskApply(tasks, fun, a=a)
  print(out)

  expect_equal(out,  list(3,6,9,12))
})

test_that("TaskApply works with simple function function without optional orguments", {
  tasks = list(1,2,3,4)
  print(length(tasks))
  start_mpi()
  fun = function(x)
  {
    print(x)
    return(x*2)
  }
  out = TaskApply(tasks, fun)
  print(out)

  expect_equal(out,  list(2,4,6,8))
})




