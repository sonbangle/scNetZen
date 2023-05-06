

#' Title
#'
#' @param aggregate_out_dir
#'
#' @return
#'
#' @examples
get_barcode_names = function(aggregate_out_dir, barcode_length=16)
  {
  ges_list = read_ges(aggregate_out_dir)
  ges = ges_list[[1]]
  feature.names = ges_list[[2]]
  barcode.names = ges_list[[3]]
  barcode.names$samples = substr(barcode.names$V1, barcode_length + 2 , nchar(barcode.names$V1))
  return(barcode.names)
}

get_barcode_names_from_ges_list = function(ges_list)
{
  barcode.names = ges_list[[3]]
  barcode.names$samples = substr(barcode.names$V1, 18, nchar(barcode.names$V1))
  return(barcode.names)
}


#' Title
#'
#' @param aggregation_samples_file
#' @param barcode.names
#' @return
#' @examples
get_task_from_aggregation_samples = function(aggregation_samples_file ="aggregation.csv", barcode.names)
{
  aggregation_samples <-
    read.csv(aggregation_samples_file,
             stringsAsFactors = FALSE)
  colnames(aggregation_samples)[1] = "library_id" # Change the first column name to library_id as it can be sample_id initially.
  samples = aggregation_samples$library_id
  n_samples = length(samples)
  print(paste("Total number of samples", n_samples))
  tasks = as.list(integer(n_samples))
  for (i in c(1:n_samples))
  {
    sample = samples[i]
    sample_barcodes = barcode.names[barcode.names$samples == i, 1]

    tasks[[i]] = list(
      barcodes = sample_barcodes,
      barcode_group_name = sample
    )
  }
  return(tasks)
}



#' Title
#'
#' @param aggregation_samples_file
#' @param barcode.names
#' @return list of tasks and sample_cluster_matrix where tasks is a list of barcodes, barcode_group_name, sample_cluster_matrix is matrix of number of cells in each sample cluster.
#' Rows are clusters, columsn are samples
#' @examples
get_task_for_samples_clusters = function(aggregation_samples_file , clusters_file,  barcode.names)
{
  #samples
  aggregation_samples <-
    read.csv(aggregation_samples_file,
             stringsAsFactors = FALSE)
  colnames(aggregation_samples)[1] = "library_id" # Change the first column name to library_id as it can be sample_id initially.
  samples = aggregation_samples$library_id
  n_samples = length(samples)
  print(paste("Total number of samples", n_samples))

  #cluster
  clusters <- read.csv(clusters_file, stringsAsFactors = FALSE)
  clusters = as.data.frame(clusters)
  colnames(clusters) = c("Barcode", "Cluster")
  clusters = as.data.frame(clusters)
  clusters[clusters$Cluster == "", 2] = "Unassigned"
  # Remove space in the Cluster names
  clusters$Cluster = gsub(" ", "_", clusters$Cluster)
  uniqclusters = unique(clusters$Cluster)
  n_clusters =  length(uniqclusters)


  sample_cluster_matrix = matrix(data = 0,
                                 nrow = n_clusters,
                                 ncol = n_samples)
  sample_cluster_matrix = data.frame(sample_cluster_matrix)
  colnames(sample_cluster_matrix) = samples
  uniqclusters = as.vector(uniqclusters)
  uniqclusters = sort(as.numeric(uniqclusters))
  rownames(sample_cluster_matrix) = uniqclusters


  tasks = as.list(integer(n_samples * n_clusters))
  task_id = 1
  for (i in c(1:n_samples))
  {
    sample = samples[i]
    sample_barcodes = barcode.names[barcode.names$samples == i, 1]

      cluster_assignments  = clusters[clusters$Barcode %in% sample_barcodes,] # cluster in this particular sample
      for (j  in c(1:n_clusters))
      {
        cluster = uniqclusters[j]
        cluster_barcodes = as.vector(cluster_assignments[cluster_assignments$Cluster == cluster, "Barcode"])
        n_barcodes = length(cluster_barcodes)
        sample_cluster_matrix[j, i] = n_barcodes
        tasks[[task_id]] = list(
          barcodes = cluster_barcodes,
          barcode_group_name = paste0(sample, "^cl_", cluster)
        )
        task_id = task_id +1
      }

  }

  return(list(tasks=tasks, sample_cluster_matrix = sample_cluster_matrix))

}



get_task_from_annotation = function(annotation_file ="~/Downloads/count/Annotation.csv")

{
  Annotation <- read.csv(annotation_file)
  Annotation[Annotation$Annotation == "", "Annotation"]= "Unassigned"
  cell_groups = unique(Annotation$Annotation)
  n_groups = length(cell_groups)
  print(paste("Total number of groups", n_groups))
  tasks = as.list(integer(n_groups))
  for (i in c(1:n_groups))
  {
    cell_group = cell_groups[i]
    group__barcodes = Annotation[Annotation$Annotation == cell_group, 1]

    tasks[[i]] = list(
      barcodes = group__barcodes,
      barcode_group_name = cell_group
    )
  }
  return(tasks)
}



