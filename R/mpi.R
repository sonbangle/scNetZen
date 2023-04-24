
#' Title
#'
#' @return ns number of slave
#'
#' @examples
#'
#'
start_mpi = function()
{

  # run mpi

  print("beginning mpi")
  ns <- mpi.universe.size() - 1

  if (ns == 0)
    #desktop
  {
    ns = 2
  }

  print(paste("Number of slaves:", ns))

  n_current_slaves =  mpi.comm.size() -1
  if (n_current_slaves <0)
    {
    n_current_slaves = 0
  }
  print(paste("Current number of unclosed slaves:",  n_current_slaves))
  n_slaves_to_spawn = ns - n_current_slaves
  if (n_slaves_to_spawn > 0)
  {
    mpi.spawn.Rslaves(nslaves = n_slaves_to_spawn)

  }
  n_current_slaves =  mpi.comm.size() -1
  print(paste("Current number of unclosed slaves:",  n_current_slaves))

  # Tell all slaves to return a message identifying themselves
  mpi.bcast.cmd(id <- mpi.comm.rank())
  mpi.bcast.cmd(print(paste("rank:", id)))
  mpi.bcast.cmd(ns <- mpi.comm.size())
  mpi.bcast.cmd(host <- mpi.get.processor.name())
  mpi.bcast.Robj2slave(slave_control)
  #mpi.bcast.Robj2slave(read_ges)
  #mpi.bcast.Robj2slave(extract_task_ges)
  print("Done starting mpi")

  return(ns)
}



