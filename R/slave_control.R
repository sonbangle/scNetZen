slave_control = function()
  #Manage the slave processing
{
  # sent:  1=ready_for_task, 2=done_task, 3=exiting received:  1=task, 2=done_tasks
  ges_list = read_ges(aggregate_out_dir) # Each slave read the data from the hard disk as rmpi cannot pass huge dataset over the network. Error Long Vector not supported
  ges = ges_list[[1]]

  junk <- 0
  done <- 0 # Flag if all tasks are done
  while (done != 1) {
    # Signal being ready to receive a new task
    mpi.send.Robj(junk, 0, 1) # Send to Master: ready for task
    task <-
      mpi.recv.Robj(mpi.any.source(), mpi.any.tag()) # Receive a task
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]
    if (tag == 1) {
      #task exits
      print("processing task")
      results = extract_cluster_ges(task, ges)
      print("done for this task, sending back result")
      # Send a results message back to the master, tag: done_task
      mpi.send.Robj(results, 0, 2)
    }
    else if (tag == 2) {
      done <- 1
    } # We'll just ignore any unknown messages
  }
  mpi.send.Robj(junk, 0, 3)
  # Exiting received

}
