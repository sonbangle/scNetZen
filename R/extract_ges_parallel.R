

#' Title
#'
#' @param aggregate_out_dir Aggregation outdir, parent folder of folder filtered_feature_bc_matrix
#' @param n_replicates  number of replicates to sample from each sing cell group
#'
#' @return
#' @export
#'
#' @examples
load_ges_slave = function(aggregate_out_dir)
{
  ns = start_mpi()
  mpi.bcast.Robj2slave(aggregate_out_dir)
  mpi.bcast.Robj2slave(read_ges)
  mpi.bcast.cmd(ges<-read_ges(aggregate_out_dir)[[1]], nonblock = FALSE) # should use ges<-read_ges(aggregate_out_dir), not ges=read_ges(aggregate_out_dir) otherwise it will generate error
}




#' Title
#'
#' @param aggregate_out_dir
#' @param min_n_reads_per_cluster
#' @param n_replicates
#' @param outdir
#' @param organism
#' @param clusters_file
#' @param extract_cell_group
#' @param method method of parallel processing, choose from :mpi, multicore, singlecore.
#'
#' @return
#' @export
#'
#' @examples
#'
extract_ges_parallel = function(tasks,
                           method= "mpi",
                           close_slaves = TRUE,
                           n_replicates = 5)


{
  # This function extracts the percentage and absolute counts of cells in each cluster and sammple, raw number of reads per cluster per sample  and normalized count per million(cpm) per each cluster per each sample
  # Input aggreated_out_dir: "outs" Folder of aggreated sample, parent of filtered_feature_bc_matrix
  # Clusters_file: file with cluster annotation, by default using 10X cluster but, can be used to provided customised clusters such as from Seurat. Table with 2 columns Barcode, Cluster with header.
  #extract_cell_group: TRUE: extract  cluster expression for each sample, cluster. FALSE: extract only sample expression)
  # Happen at the Master Only
  # This is function of master node, rank 0

  outs = sonicApply(x =tasks , fun =extract_task_with_ges_predefine, method =method, close_slaves= close_slaves, n_replicates = n_replicates)
  return(outs)

}

make_count_data_frame_from_extract_ges_obj = function(extract_ges_obj )
{
  # From the result of extract ges parallel, format the results as two dataframe, one of all cells in cell group, another for cell replicates)

  n_cell_groups = length(extract_ges_obj)
  res1 = extract_ges_obj[[1]]
  n_genes = nrow(res1[["barcodes_ges_sum"]])
  genes =  rownames(res1[["barcodes_ges_sum"]])
  n_replicates = ncol(res1[["barcodes_ges_sum_replicates"]])
  cell_group_ges = matrix(data = 0,
                          nrow = n_genes ,
                          ncol = n_cell_groups + 1)
  cell_group_ges = data.frame(cell_group_ges)
  rownames(cell_group_ges) = genes
  cell_group_ges[, 1] = genes
  colnames(cell_group_ges)[1] = "SYMBOL"
  cell_group_ges_replicates = matrix(
    data = 0,
    nrow = n_genes ,
    ncol = n_cell_groups * n_replicates + 1
  )
  cell_group_ges_replicates = data.frame(cell_group_ges_replicates)
  rownames(cell_group_ges_replicates) = genes
  cell_group_ges_replicates[, 1] = genes
  colnames(cell_group_ges_replicates)[1] = "SYMBOL"
  ges_col_index = 2
  ges_col_index_replicates = 2
  for (cell_group_profile in extract_ges_obj )
  {
    barcodes_ges_sum = cell_group_profile[["barcodes_ges_sum"]]
    barcodes_ges_sum_replicates = cell_group_profile[["barcodes_ges_sum_replicates"]]
    cell_group_ges[, ges_col_index] = barcodes_ges_sum
    colnames(cell_group_ges)[ges_col_index] = colnames(barcodes_ges_sum)
    ges_col_index_replicates_end = ges_col_index_replicates + n_replicates -1
    cell_group_ges_replicates[, ges_col_index_replicates:ges_col_index_replicates_end] = barcodes_ges_sum_replicates
    colnames(cell_group_ges_replicates)[ges_col_index_replicates:ges_col_index_replicates_end] = colnames(barcodes_ges_sum_replicates)
    ges_col_index = ges_col_index +1
    ges_col_index_replicates = ges_col_index_replicates + n_replicates

  }

  return(list(counts= cell_group_ges, count_replicates = cell_group_ges_replicates))

}


get_count_from_tasks = function(tasks,
                                method= "mpi",
                                close_slaves = TRUE,
                                n_replicates = 5)
{

  extract_ges_obj = extract_ges_parallel(tasks,
                                         method,
                                         close_slaves,
                                         n_replicates )
  counts = make_count_data_frame_from_extract_ges_obj(extract_ges_obj)
  return(counts)
}






