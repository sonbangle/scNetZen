

#' Title
#'
#' @param task the task to calculate expression, a list of barcodes, n_replicates, barcode_group_name
#' @param ges
#' @param min_n_cells_per_cluster
#'
#' @return a list ( result: list of (barcodes_ges_sum dataframe, barcodes_ges_sum_replicates dataframe, is_null_result), task_id)
#'
#' @examples
extract_task_ges = function(task,
                            ges, n_replicates=5)
  # Slave function
  #Input: a list (sample, cluster, barcodes)
  #Output: list(cluster_ges_sum, cluster_ges_sum_replicates)

{
  barcodes = task[["barcodes"]]
  barcode_group_name = task[["barcode_group_name"]]
  out = extract_ges(barcodes = barcodes,
                    ges = ges,
                    barcode_group_name = barcode_group_name,
                    n_replicates = n_replicates
  )

  return(out)


}





extract_task_with_ges_predefine <- function(task, n_replicates =5,  ges_name = "ges") {
  # Attempt to retrieve the gene expression matrix (ges) using dynamic scoping
  tryCatch({
    ges <- dynGet(ges_name)
  }, error = function(e) {
    stop("Gene expression matrix does not exist")
  })

  # Continue with the execution if ges is found
  return(extract_task_ges(task, ges, n_replicates))
}

#' Title
#'Extract gene expression from group of cells, such as a cluster of cells, the whole sample or any group of cells
#' @param barcodes the barcodes of single cells in the group of interest.
#' @param ges gene expression profile, a sparse matrix of class dgTMatrix or data frame, row are genes, columns are cell barcodes
#' @param barcode_group_name barcode group name
#' @param task_id tak_id
#' @param min_n_cells minimal number of cells in the group to calculate
#' @param n_replicates number of replicates when do sampling of cells in the group
#'
#' @return
#' @export
#'
#' @examples
#'

extract_ges = function(barcodes,
                        ges,
                       barcode_group_name,
                       n_replicates = 5
                       )
  # Slave function
  #Input: a list (sample, cell group, barcodes)
  #Output: list(cell group_ges_sum, cluster_ges_sum_replicates)

{

 # print(paste("I am", mpi.comm.rank(), "of", mpi.comm.size()))

    barcodes_ges_sum_replicates = data.frame(matrix(0., nrow = nrow(ges), ncol = n_replicates))

    if (length(barcodes) == 0)
    {
      barcodes_ges_sum = data.frame(barcode_group_name = rep(0,nrow(ges)))
      rownames(barcodes_ges_sum) = rownames(ges)
      colnames(barcodes_ges_sum_replicates) = sapply(c(1:n_replicates) , function(k) paste0(barcode_group_name, "^rep_", k))
      rownames(barcodes_ges_sum_replicates) = rownames(ges)

      return(list(barcodes_ges_sum = barcodes_ges_sum,
             barcodes_ges_sum_replicates = barcodes_ges_sum_replicates,
             is_null_result = 1))
    }
    barcodes_ges = ges[, barcodes, drop=FALSE]
    #barcodes_ges = as.data.frame(barcodes_ges)
    barcodes_ges_sum = Matrix::rowSums(barcodes_ges)
    barcodes_ges_sum = data.frame(barcodes_ges_sum)
    colnames(barcodes_ges_sum) = barcode_group_name

    #make replicate count matrix
    #barcodes_ges_sum_replicates = data.frame(matrix(NA, nrow = nrow(barcodes_ges), ncol = n_replicates))

    for (k in c(1:n_replicates))
    {
      barcodes_replicate = sample(barcodes,
                                          size = length(barcodes),
                                          replace = TRUE)
      barcodes_ges_replicate = barcodes_ges[, barcodes_replicate, drop = FALSE]
      barcodes_ges_sum_replicate = Matrix::rowSums(barcodes_ges_replicate)
      barcodes_ges_sum_replicates[, k] = barcodes_ges_sum_replicate
      colnames(barcodes_ges_sum_replicates)[k] = paste0(barcode_group_name, "^rep_", k)
    }
    rownames(barcodes_ges_sum_replicates) = rownames(barcodes_ges)
    return(
      list(
        barcodes_ges_sum = barcodes_ges_sum,
        barcodes_ges_sum_replicates = barcodes_ges_sum_replicates,
        is_null_result = 0
      )
    )

}
