

#' Title
#'from count data frame of cell group gene expression, filter cell groups that have number of read higher than threshold, and the total number of reads for each cell group
#' @param cell_group_ges a count dataframe where first column is ENS id, other columns are cell group names, rows are genes, values are counts
#' @param gene_names the corresponding normal gene names
#' @param min_n_reads_per_cell_group Minimal number of reads per cell group threshold to filter
#'
#' @return  a list of:
#'  cell_group_ges:original count data frame
#' cell_group_ges_cpm : cpm gene expression dataframe
#' cell_group_ges_hugo :count data frame with normal gene names instead of ENS id
#' cell_group_ges_hugo_cpm : cpm gene expression data frame with normal gene names instead of ENS id
#' filtered_cell_group_ges :  filtered count data frame, where cell groups have total number of reads more than min_n_reads_per_cell_group
#' filtered_cell_group_ges_hugo :  filtered count data frame with normal gene names instead of ENS id, where cell groups have total number of reads more than min_n_reads_per_cell_group
#' total_n_reads_per_cell_group : total number of reads per cell group data frame with two columns:  total_n_reads_per_cell_group,  cell_group.
#' total_n_reads_per_cell_group_more_than_threshold :  total number of reads per cell group that pass threshold data frame with two columns:  total_n_reads_per_cell_group,  cell_group.
#'
#' @export
#'
#' @examples
#'
get_formated_ges_from_counts = function(cell_group_ges,   gene_names, min_n_reads_per_cell_group= 20000)
{
  total_n_reads_per_cell_group = colSums(cell_group_ges[2:ncol(cell_group_ges)])
  total_n_reads_per_cell_group = as.data.frame(total_n_reads_per_cell_group)
  total_n_reads_per_cell_group$cell_group = rownames(total_n_reads_per_cell_group)
  total_n_reads_per_cell_group = total_n_reads_per_cell_group[, c(2,1)]
  total_n_reads_per_cell_group_more_than_threshold = as.data.frame(total_n_reads_per_cell_group[total_n_reads_per_cell_group$total_n_reads_per_cell_group > min_n_reads_per_cell_group,])
  cell_group_ges_cpm = cell_group_ges
  for (i in c(2:ncol(cell_group_ges_cpm)))
  {
    cell_group_ges_cpm[, i] =  cell_group_ges[, i] / (total_n_reads_per_cell_group[i -1, 2] + 1) * 1e6
  }
  cell_group_ges_hugo = cell_group_ges
  cell_group_ges_hugo[, 1] = gene_names
  cell_group_ges_hugo_cpm = cell_group_ges_cpm
  cell_group_ges_hugo_cpm[, 1] = gene_names

  filtered_cell_groups = total_n_reads_per_cell_group_more_than_threshold$cell_group  # Name of cell groups  that have total number of counts more than threshold.
  filtered_cell_group_ges = cell_group_ges[, c("SYMBOL", filtered_cell_groups)]
  filtered_cell_group_ges_hugo = cell_group_ges_hugo[, c("SYMBOL", filtered_cell_groups)]

    return(list(cell_group_ges = cell_group_ges,
              cell_group_ges_cpm = cell_group_ges_cpm,
              cell_group_ges_hugo = cell_group_ges_hugo,
              cell_group_ges_hugo_cpm = cell_group_ges_hugo_cpm,
              filtered_cell_group_ges = filtered_cell_group_ges,
              filtered_cell_group_ges_hugo = filtered_cell_group_ges_hugo,
              total_n_reads_per_cell_group = total_n_reads_per_cell_group,
              total_n_reads_per_cell_group_more_than_threshold = total_n_reads_per_cell_group_more_than_threshold))
}

#' Title
#' Save the result of get_formated_ges_from_counts function to disk , also get normalized gene expression by calling the normalized_and_export_counts_data function
#' @param formated_ges_list the result of get_formated_ges_from_counts
#' @param outdir output directory
#' @param prefix prefix to add to the file path, such as sample, sample_replicates, sample_cluster, sample_cluster_duplicates, etc.
#' @param organism either mouse or human . Use to translate gene ENSM id into normal gene names
#'
#' @return
#' @export
#'
#' @examples
normalize_and_write_result = function(formated_ges_list , outdir, prefix  = "sample", organism = "human")
{
  # Write result from  get_formated_ges_from_counts function to disk, normalize the filtered_cell_group_ges and filtered_cell_group_replicates_ges and write to disk

  dir.create(outdir)
  outdir = paste0(outdir, "/", toupper(prefix))
  dir.create(outdir)

  write.table(
    formated_ges_list[["cell_group_ges"]],
    file = paste0(outdir, "/", prefix, "_ges.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    formated_ges_list[["cell_group_ges_cpm"]],
    file = paste0(outdir, "/", prefix, "_ges_cpm.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  write.table(
    formated_ges_list[["cell_group_ges_hugo"]],
    file = paste0(outdir, "/", prefix, "_ges_translated.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    formated_ges_list[["cell_group_ges_hugo_cpm"]],
    file = paste0(outdir, "/", prefix, "_ges_cpm_translated.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  write.table(
    formated_ges_list[["filtered_cell_group_ges_hugo"]],
    file = paste0(outdir, "/", prefix, "_ges_translated.filtered.csv"),
    quote = FALSE,
    row.names = FALSE
  )


  write.table(
    formated_ges_list[["total_n_reads_per_cell_group_more_than_threshold"]],
    file = paste0(
      outdir,
      "/total_n_reads_per_", prefix, "_more_than_threshold.csv"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )


  prefix = toupper(prefix)
  count_outdir = paste0(outdir, "/", prefix, "_COUNTS")

  normalized_and_export_counts_data(
    counts =   formated_ges_list[["filtered_cell_group_ges"]],
    outdir = count_outdir,
    organism = organism
  )

}









