
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
  print(sample)
  print(cluster)

  print(paste("I am", mpi.comm.rank(), "of", mpi.comm.size()))
  if (length(cluster_barcodes) >=  min_n_cells_per_cluster)
    # At least more than x cells in clusters
    # eliminate the cluster that has less than predefined number of  cells  in this cluster
  {
    cluster_ges = ges[, cluster_barcodes]
    cluster_ges_sum = rowSums(cluster_ges)
    # sample_cluster_ges[, ges_col_index] = cluster_ges_sum
    # colnames(sample_cluster_ges)[ges_col_index] = paste0(sample, "^", cluster)
    # ges_col_index = ges_col_index + 1
    cluster_ges_sum = data.frame(cluster_ges_sum)
    colnames(cluster_ges_sum) = c(paste0(sample, "^", cluster))

    #make replicate count matrix
    cluster_ges_sum_replicates = data.frame(matrix(NA, nrow = nrow(cluster_ges), ncol = n_replicates))

    for (k in c(1:n_replicates))
    {
      print(paste("replicate", k))
      cluster_barcodes_replicate = sample(cluster_barcodes,
                                          size = length(cluster_barcodes),
                                          replace = TRUE)
      cluster_ges_replicate = cluster_ges[, cluster_barcodes_replicate]
      cluster_ges_sum_replicate = rowSums(cluster_ges_replicate)
      cluster_ges_sum_replicates[, k] = cluster_ges_sum_replicate
      colnames(cluster_ges_sum_replicates)[k] = paste0(sample, "^" , cluster, "^rep_", k)
    }
    print("done")
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

  } else
  {
    print(paste("too few cells in this cluster:", length(cluster_barcodes)))
    return(
      list(
        sample = sample,
        cluster = cluster,
        cluster_ges_sum = NULL,
        cluster_ges_sum_replicates = NULL,
        sample_index = sample_index,
        cluster_index = cluster_index,
        is_null_result = 1
      )
    )
  }

}
