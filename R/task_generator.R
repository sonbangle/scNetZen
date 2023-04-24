

#' Title
#'
#' @param aggregate_out_dir
#'
#' @return
#'
#' @examples
get_data = function(aggregate_out_dir)
  {
  ges_list = read_ges(aggregate_out_dir)
  ges = ges_list[[1]]
  feature.names = ges_list[[2]]
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
      barcode_group_name = sample,
      n_replicates = n_replicates
    )
  }
  return(tasks)
}


