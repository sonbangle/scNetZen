% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/task_generator.R
\name{get_task_for_samples_clusters}
\alias{get_task_for_samples_clusters}
\title{Title}
\usage{
get_task_for_samples_clusters(
  aggregation_samples_file,
  clusters_file,
  barcode.names
)
}
\value{
list of tasks and sample_cluster_matrix where tasks is a list of barcodes, barcode_group_name, sample_cluster_matrix is matrix of number of cells in each sample cluster.
Rows are clusters, columsn are samples
}
\description{
Title
}
