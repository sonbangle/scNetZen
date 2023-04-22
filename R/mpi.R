
#' Title
#'
#' @return ns number of slave
#'
#' @examples
#'
#'
start_mpi = function()
{
  print("github updated ")

  print("Starting extraction single cell gene expresion profiles")

  # run mpi

  print("beginning mpi")
  closed_slaves <- 0
  ns <- mpi.universe.size() - 1
  print(paste("Number of slaves:", ns))
  if (ns == 0)
    #desktop
  {
    ns = 1
  }

  n_current_slaves =  mpi.comm.size()
  print(paste("Current number of unclosed slaves:",  n_current_slaves))
  n_slaves_to_spawn = ns - n_current_slaves
  if (n_slaves_to_spawn > 0)
  {
    mpi.spawn.Rslaves(nslaves = n_slaves_to_spawn)

  }

  # Tell all slaves to return a message identifying themselves
  mpi.bcast.cmd(id <- mpi.comm.rank())
  mpi.bcast.cmd(ns <- mpi.comm.size())
  mpi.bcast.cmd(host <- mpi.get.processor.name())
  mpi.bcast.cmd(library("Matrix"))
  return(ns)
}
