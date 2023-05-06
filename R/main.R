#' Title
#' Extract aggregated gene expression of samples in scRNA seq, convert from single cell RNAseq to bulk RNAseq.
#' @param aggregate_out_dir aggregate output directory from 10X run, parent directory of "filtered_feature_bc_matrix"
#' @param aggregation_samples_file aggregation sample file for 10X pipeline, containing columns: library_id, molecule_h5 and other metadata column
#' @param outdir output directory
#' @param organism  either mouse or human . Use to translate gene ENSM id into normal gene names
#' @param method method to run parallel processing. Choose from "mpi", "multicore", "singlecore"
#' @param close_slaves default = TRUE. if method is mpi, then close slaves node
#' @param n_replicates number of replicates to sample from
#' @param min_n_reads_per_cell_group minimal number of reads per cell group (sample)
#'
#' @return write in output directory the pseudo bulk RNA seq with following format: count, tpm, rpkm, cpm.
#' @export
#'
#' @examples
get_aggregation_sample_gene_expressions = function(aggregate_out_dir,
                                                  aggregation_samples_file,
                                                  outdir,
                                                  organism,
                                                  method = "mpi",
                                                  close_slaves = TRUE,
                                                  n_replicates = 5,
                                                  min_n_reads_per_cell_group = 20000)

{
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
  prefix = "aggregation_sample"
  out = get_cell_group_gene_expression_from_tasks(tasks,
                                                  aggregate_out_dir,
                                                  outdir,
                                                  prefix,
                                                  organism,
                                                  method = method,
                                                  close_slaves = close_slaves,
                                                  n_replicates = n_replicates,
                                                  min_n_reads_per_cell_group = min_n_reads_per_cell_group)

  }


#' Title
#'  Extract aggregated gene expression of cell group in scRNA seq, convert from single cell RNAseq to bulk RNAseq. Each cell group is a cluster of individual sample.
#' @param aggregate_out_dir aggregate output directory from 10X run, parent directory of "filtered_feature_bc_matrix"
#' @param aggregation_samples_file aggregation sample file for 10X pipeline, containing columns: library_id, molecule_h5 and other metadata column
#' @param clusters_file clusters file. A table with two columns: Barcode , Cluster
#' @param outdir  output directory
#' @param organism either mouse or human . Use to translate gene ENSM id into normal gene names
#' @param method method to run parallel processing. Choose from "mpi", "multicore", "singlecore"
#' @param close_slaves default = TRUE. if method is mpi, then close slaves node
#' @param n_replicates  number of replicates to sample from
#' @param min_n_reads_per_cell_group minimal number of reads per cell group (sample_cluster)
#'
#' @return
#' @export
#'
#' @examples
get_sample_cluster_gene_expressions = function(aggregate_out_dir,
                                                   aggregation_samples_file,
                                                    clusters_file,
                                                   outdir,
                                                   organism,
                                                   method = "mpi",
                                                   close_slaves = TRUE,
                                                   n_replicates = 5,
                                                   min_n_reads_per_cell_group = 20000)

{
  barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
  tasks = get_task_for_samples_clusters(aggregation_samples_file = aggregation_samples_file , clusters_file = clusters_file,  barcode.names = barcode.names)[["tasks"]]
  prefix = "sample_cluster"
  out = get_cell_group_gene_expression_from_tasks(tasks,
                                                  aggregate_out_dir,
                                                  outdir,
                                                  prefix,
                                                  organism,
                                                  method = method,
                                                  close_slaves = close_slaves,
                                                  n_replicates = n_replicates,
                                                  min_n_reads_per_cell_group = min_n_reads_per_cell_group)

}



#' Title
#' Extract aggregated gene expression of cell group in scRNA seq, convert from single cell RNAseq to bulk RNAseq. Each cell group is a cell group in annotated file.
#' @param aggregate_out_dir aggregate output directory from 10X run, parent directory of "filtered_feature_bc_matrix"
#' @param annotation_file aggregation sample file for 10X pipeline, containing columns: library_id, molecule_h5 and other metadata column
#' @param outdir  output directory
#' @param organism either mouse or human . Use to translate gene ENSM id into normal gene names
#' @param method  method to run parallel processing. Choose from "mpi", "multicore", "singlecore"
#' @param close_slaves  default = TRUE. if method is mpi, then close slaves node
#' @param n_replicates  number of replicates to sample from
#' @param min_n_reads_per_cell_group minimal number of reads per cell group (sample_cluster)
#' @param prefix prefix to add to the output subdirectories and files. For example: "immunocell_annotation"
#'
#' @return
#' @export
#'
#' @examples
get_cell_group_gene_expressions_from_annotation = function(aggregate_out_dir,
                                                           annotation_file,
                                               outdir,
                                               organism,

                                               method = "mpi",
                                               close_slaves = TRUE,
                                               n_replicates = 5,
                                               min_n_reads_per_cell_group = 20000,
                                               prefix = "cell_group")
{
  tasks = get_task_from_annotation(annotation_file)
  out = get_cell_group_gene_expression_from_tasks(tasks,
                                                  aggregate_out_dir,
                                                  outdir,
                                                  prefix,
                                                  organism,
                                                  method = method,
                                                  close_slaves = close_slaves,
                                                  n_replicates = n_replicates,
                                                  min_n_reads_per_cell_group = min_n_reads_per_cell_group)

}

