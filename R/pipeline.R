#' Title
#'
#' @param tasks
#' @param aggregate_out_dir
#' @param outdir
#' @param prefix
#' @param organism
#' @param method
#' @param close_slaves
#' @param n_replicates
#' @param min_n_reads_per_cell_group
#'
#' @return
#' @export
#'
#' @examples
get_cell_group_gene_expression_from_tasks = function(tasks,
                                                     aggregate_out_dir,
                                                     outdir,
                                                     prefix,
                                                     organism,
                                                     method = "mpi",
                                                     close_slaves = TRUE,
                                                     n_replicates = 5,
                                                     min_n_reads_per_cell_group = 20000)
{
  ges_list <- read_ges(aggregate_out_dir)[[1]]
  ges = ges_list[[1]]
  gene_names = ges_list[[2]][2]
  if (method == "mpi")
  {
    load_ges_slave(aggregate_out_dir)
  }

  counts = get_count_from_tasks(
    tasks,
    method = method,
    close_slaves = close_slaves,
    n_replicates = n_replicates
  )
  cell_group_counts = counts[[1]]
  cell_group_counts_replicates =  counts[[2]]
  to_write_data = list(list(prefix, cell_group_counts),
                       list(
                         paste0(prefix, "_replicates"),
                         cell_group_counts_replicates
                       ))
  lapply(to_write_data, function(x) {
    write_cell_group_gene_expression_from_counts(x[[2]],
                                                 gene_names,
                                                 outdir,
                                                 x[[1]],
                                                 organism,
                                                 min_n_reads_per_cell_group = 20000)
  })
}




#' Title
#' From cell group counts data frame, write different cell group gene_expression normalized data, and filtered for the cell groups with total numbers of reads above threshold
#' @param cell_group_counts count dataframe where first column is ENS id, other columns are cell group names, rows are genes, values are counts
#' @param gene_names the corresponding normal gene names
#' @param outdir output directory
#' @param prefix prefix to add to the file path, such as sample, sample_replicates, sample_cluster, sample_cluster_duplicates, etc.
#' @param organism either mouse or human . Use to translate gene ENSM id into normal gene names
#' @param min_n_reads_per_cell_group Minimal number of reads per cell group threshold to filter
#'
#' @return
#' @export
#'
#' @examples
write_cell_group_gene_expression_from_counts = function(cell_group_counts,
                                                        gene_names,
                                                        outdir,
                                                        prefix,
                                                        organism,
                                                        min_n_reads_per_cell_group = 20000)
{
  out = get_formated_ges_from_counts(cell_group_counts,  gene_names, min_n_reads_per_cell_group = 20000)
  normalize_and_write_result(out , outdir = outdir , prefix  = prefix, organism)
  return(out)
}
