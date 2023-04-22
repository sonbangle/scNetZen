
extract_cluster_ges = function(sample_cluster_barcode,
                               ges,
                               min_n_cells_per_cluster = 5)
  # Slave function
  #Input: a list (sample, cluster, barcodes)
  #Output: list(cluster_ges_sum, cluster_ges_sum_replicates)

{
  sample = sample_cluster_barcode[["sample"]]
  cluster = sample_cluster_barcode[["cluster"]]
  cluster_barcodes = sample_cluster_barcode[["cluster_barcodes"]]
  sample_index = sample_cluster_barcode[["i"]]
  cluster_index = sample_cluster_barcode[["j"]]
  n_replicates = sample_cluster_barcode[["n_replicates"]]

  out = extract_ges(barcodes = cluster_barcodes,
                     ges = ges,
                     barcode_group_name = paste0(sample, "^", cluster),
                     min_n_cells = min_n_cells_per_cluster,
                     n_replicates = n_replicates
  )

  cluster_ges_sum = out[["barcodes_ges_sum"]]
  cluster_ges_sum_replicates = out[["barcodes_ges_sum_replicates"]]
  is_null_result = out[["is_null_result"]]
    return(
      list(
        sample = sample,
        cluster = cluster,
        cluster_ges_sum = cluster_ges_sum,
        cluster_ges_sum_replicates = cluster_ges_sum_replicates,
        sample_index = sample_index,
        cluster_index = cluster_index,
        is_null_result = 0
      )
    )

}



#' Title
#'
#' @param task the task to calculate expression, a list of barcodes, n_replicates, task_id, barcode_group_name
#' @param ges
#' @param min_n_cells_per_cluster
#'
#' @return a list ( result: list of (barcodes_ges_sum dataframe, barcodes_ges_sum_replicates dataframe, is_null_result), task_id)
#'
#' @examples
extract_task_ges = function(task,
                            ges,
                            min_n_cells_per_cluster = 5)
  # Slave function
  #Input: a list (sample, cluster, barcodes)
  #Output: list(cluster_ges_sum, cluster_ges_sum_replicates)

{

  barcodes = task[["barcodes"]]
  n_replicates = task[["n_replicates"]]
  task_id = task[["task_id"]]
  barcode_group_name = task[["barcode_group_name"]]
  out = extract_ges(barcodes = barcodes,
                    ges = ges,
                    barcode_group_name = barcode_group_name,
                    min_n_cells = min_n_cells_per_cluster,
                    n_replicates = n_replicates
  )

  return(
    list(result = out,
         task_id = task_id
    )
  )


}



#' Title
#'Extract gene expression from group of cells, such as a cluster of cells, the whole sample or any group of cells
#' @param barcodes the barcodes of single cells in the group of interest.
#' @param ges gene expression profile
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
                        min_n_cells = 5,
                       n_replicates = 5,
                       )
  # Slave function
  #Input: a list (sample, cell group, barcodes)
  #Output: list(cell group_ges_sum, cluster_ges_sum_replicates)

{

  print(paste("I am", mpi.comm.rank(), "of", mpi.comm.size()))

  if (length(barcodes) >=  min_n_cells)
    # At least more than x cells in barcodes group
    # eliminate the cell group that has less than predefined number of  cells  in this cell group
  {
    barcodes_ges = ges[, barcodes]
    barcodes_ges_sum = rowSums(barcodes_ges)
    barcodes_ges_sum = data.frame(barcodes_ges_sum)
    colnames(barcodes_ges_sum) = barcode_group_name

    #make replicate count matrix
    barcodes_ges_sum_replicates = data.frame(matrix(NA, nrow = nrow(barcodes_ges), ncol = n_replicates))

    for (k in c(1:n_replicates))
    {
      print(paste("replicate", k))
      barcodes_replicate = sample(barcodes,
                                          size = length(barcodes),
                                          replace = TRUE)
      barcodes_ges_replicate = barcodes_ges[, barcodes_replicate]
      barcodes_ges_sum_replicate = rowSums(barcodes_ges_replicate)
      barcodes_ges_sum_replicates[, k] = barcodes_ges_sum_replicate
      colnames(barcodes_ges_sum_replicates)[k] = paste0(barcode_group_name, "^rep_", k)
    }
    print("done")
    return(
      list(
        barcodes_ges_sum = barcodes_ges_sum,
        barcodes_ges_sum_replicates = barcodes_ges_sum_replicates,
        is_null_result = 0
      )
    )

  } else
  {
    print(paste("too few cells in this cluster:", length(barcodes)))
    return(
      list(
        barcodes_ges_sum = NULL,
        barcodes_ges_sum_replicates = NULL,
        is_null_result = 1
      )
    )
  }

}