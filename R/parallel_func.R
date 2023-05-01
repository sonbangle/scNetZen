
#' Title
#' A parallel version of lapply using either mpi on clusters, multicore or simple lapply
#' Advatange comparing with Rmpi mpi.applyBL: More flexible, can do some preprocessing in the slave node first, for example reading gene expression in each
#' slave node, not passing the big share data from master to slaves.
#' @param x a list of elements to process
#' @param fun function to apply
#' @param method method to do parallel processing, either mpi, multicore or single core
#' @param close_slaves close slave after doing parallel apply. default is TRUE, set to FALSE when you want to reuse slave.
#' @return
#' @export
#'
#' @examples
sonicApply = function(x, fun, method, close_slaves=TRUE, ...)
{
out = switch(method,
              "mpi"={
                start_mpi()
                out = TaskApply(tasks = x, fun = fun, ...)
                if (close_slaves) mpi.close.Rslaves()
                return(out)
              },

              "multicore"={
                out = parallel::mclapply(x, fun, ...)
                return(out)

              },
              "singlecore" ={
                out = lapply(x, fun, ...)
                return(out)
              }

)
return(out)

}
