#' Title
#'
#' @param FUN
#'
#' @return
#'
#' @examples
slave_control = function(FUN, ...)
  #Manage the slave processing
{
  # sent:  1=ready_for_task, 2=done_task, 3=exiting received:  1=task, 2=done_tasks

  junk <- 0
  done <- 0 # Flag if all tasks are done

  while (done != 1) {
    # Signal being ready to receive a new task
    mpi.send.Robj(junk, 0, 1) # Send to Master: ready for task, request task
    task <-
      mpi.recv.Robj(mpi.any.source(), mpi.any.tag()) # Receive a task
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]
    if (tag == 1) {
      #task exits
      #results = extract_task_ges(task, ges)
      task_id = task[["id"]]
      task_data = task[["data"]]
      result = FUN(task_data,...)
      result = list(result=result, id =task_id)
      # Send a results message back to the master, tag: done_task
      mpi.send.Robj(result, 0, 2) # send task result
    }
    else if (tag == 2) {
      done <- 1
    } # We'll just ignore any unknown messages
  }
  mpi.send.Robj(junk, 0, 3)
  # Exiting received

}


#' Title
#' Parallelly apply the function to list, similar to lapply but with mpi version
#' Initialization and closing mpi should be done manually so that each slave can be prepared(loading the dataset) before processing the data
#'
#' @param tasks the list of data, each element is a data to apply function
#' @param fun The name of function on the slave. Should be exist in slave node either by creating on the node or transfer from the master manually
#'
#' @return the list of concatentated output, result of applying function into each input element
#' @export
#'
#' @examples
TaskApply = function(tasks, fun)
{

  #ns <- mpi.universe.size() - 1
  ns = mpi.comm.size() -1
  print(paste("Number of slaves:",  ns))
  closed_slaves <- 0
  n_tasks = length(tasks)
  results = as.list(integer(n_tasks))
  current_task = 1
  mpi.bcast.Robj2slave(fun)
  mpi.bcast.cmd(slave_control(fun))

  while (closed_slaves < ns) {
    # Receive a message from a slave
    junk <- 0
    message <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())

    message_info <- mpi.get.sourcetag()
    slave_id <- message_info[1]
    tag <- message_info[2]
    # print(paste("get mssg from :", slave_id))
    # print(paste("tag received:", tag))
    # print("message:")
    # print(message)

    if (tag == 1) {
      # slave ready for a task. Give it the next task, or tell it tasks are done.
      if (current_task <= length(tasks)) {
        # Send a task
        task_to_send = list(data =tasks[[current_task]], id = current_task )
        mpi.send.Robj(task_to_send, slave_id, 1)
        current_task = current_task + 1

      }
      else {
        mpi.send.Robj(junk, slave_id, 2)
      }
    } else if (tag == 2) {
      # Do something with the results. Store in the data structure

    result = message[["result"]]
    task_id_received = message[["id"]]
    results[[task_id_received]] =result

    } else if (tag == 3) {
      closed_slaves <- closed_slaves + 1
    }
  }



  return(results)

}
