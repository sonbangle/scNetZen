

#' Title
#'
#' @param aggregate_out_dir Aggregation outdir, parent folder of folder filtered_feature_bc_matrix
#' @param n_replicates  number of replicates to sample from each sing cell group
#'
#' @return
#' @export
#'
#' @examples
load_ges_slave = function(aggregate_out_dir, n_replicates = 5)
{
  ns = start_mpi()
  mpi.bcast.Robj2slave(aggregate_out_dir)
  mpi.bcast.Robj2slave(n_replicates)
  mpi.bcast.Robj2slave(read_ges)
  mpi.bcast.cmd(ges = read_ges(aggregate_out_dir)[[1]]) # read gene expression profile in each slave node

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
#'
#' @return
#' @export
#'
#' @examples
extract_ges_mpi = function(tasks,
                           min_n_reads_per_cluster = 20000,
                           outdir = "NETZEN_analysis",
                           organism = "human"
)

{
  # This function extracts the percentage and absolute counts of cells in each cluster and sammple, raw number of reads per cluster per sample  and normalized count per million(cpm) per each cluster per each sample
  # Input aggreated_out_dir: "outs" Folder of aggreated sample, parent of filtered_feature_bc_matrix
  # Clusters_file: file with cluster annotation, by default using 10X cluster but, can be used to provided customised clusters such as from Seurat. Table with 2 columns Barcode, Cluster with header.
  #extract_cell_group: TRUE: extract  cluster expression for each sample, cluster. FALSE: extract only sample expression)
  # Happen at the Master Only
  # This is function of master node, rank 0

  ns = start_mpi()
  outs = TaskApply(tasks = tasks, fun = extract_task_ges_global )
  mpi.finalize()
  return(outs)




  ##############################################################Making gene expression profiles for NETZEN

  n_cell_groups = length(tasks)

  cell_group_ges = matrix(data = 0,
                          nrow = nrow(ges) ,
                          ncol = length(n_cell_groups) + 1)
  cell_group_ges = data.frame(cell_group_ges)
  rownames(cell_group_ges) = rownames(ges)
  cell_group_ges[, 1] = rownames(cell_group_ges)
  colnames(cell_group_ges)[1] = "SYMBOL"
  ges_col_index = 2


  cell_group_ges_replicates = matrix(
    data = 0,
    nrow = nrow(ges) ,
    ncol = n_cell_groups * n_replicates + 1
  )
  cell_group_ges_replicates = data.frame(cell_group_ges_replicates)
  rownames(cell_group_ges_replicates) = rownames(ges)
  cell_group_ges_replicates[, 1] = rownames(cell_group_ges_replicates)
  colnames(cell_group_ges_replicates)[1] = "SYMBOL"
  ges_col_index = 2
  ges_col_index_replicates = 2




  # Create data structures to store the results
  cell_group_ges_sum.result = list()
  cell_group_ges_sum_replicates.result = list()

  ###########################
  # send tasks to slaves

  results = as.list(integer(n_cell_groups))
  done_tasks = 0
  current_task = 1
  while (closed_slaves < ns) {
    # Receive a message from a slave
    junk <- 0
    message <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
    message_info <- mpi.get.sourcetag()
    slave_id <- message_info[1]
    tag <- message_info[2]
    print(paste("get mssg from :", slave_id))
    print(paste("tag received:", tag))
    if (tag == 1) {
      # slave ready for a task. Give it the next task, or tell it tasks are done.
      print(paste("n task remained:" , length(tasks)))
      if (current_task <= length(tasks)) {
        # Send a task
        mpi.send.Robj(tasks[[current_task]], slave_id, 1)
        current_task = current_task + 1

      }
      else {
        mpi.send.Robj(junk, slave_id, 2)
      }
    } else if (tag == 2) {
      # Do something with the results. Store in the data structure
      task_id = message[task_id]
      task = tasks[[task_id]]
      sample = task[["sample"]]
      cluster = task[["cluster"]]
      cluster_ges_sum = message[["cluster_ges_sum"]]
      cluster_ges_sum_replicates = message[["cluster_ges_sum_replicates"]]
      sample_index = task[["sample_index"]]
      cluster_index = task[["cluster_index"]]
      print(sample)
      print(cluster)
      print(paste("sample index", sample_index))
      print(paste("cluster index", cluster_index))

      pos =   (sample_index - 1) * n_clusters + cluster_index + 1

      colnames(cell_group_ges)[pos] = paste0(sample, "^", cluster)

      start_pos = (sample_index - 1) * n_clusters * n_replicates + (cluster_index -
                                                                      1) * n_replicates  + 2
      end_pos = start_pos + n_replicates - 1

      colnames(cell_group_ges_replicates)[start_pos:end_pos] = paste0(sample, "^", cluster, "^rep_", c(1:n_replicates))

      if (message[["is_null_result"]] == 0)
        # If result is not null, then add cell_group column to a data frame
      {
        cell_group_ges[, pos] = cluster_ges_sum
        cell_group_ges_replicates[, start_pos:end_pos] = cluster_ges_sum_replicates
      }

    } else if (tag == 3) {
      closed_slaves <- closed_slaves + 1
    }
  }







  # Tell all slaves to close down, and exit the program
  mpi.close.Rslaves(dellog = FALSE)

  print("Done with mpi")

  ###############################################################












  total_n_reads_per_cell_group = colSums(cell_group_ges[2:ncol(cell_group_ges)])
  total_n_reads_per_cell_group = as.data.frame(total_n_reads_per_cell_group)
  cell_group_ges_cpm = cell_group_ges
  total_n_reads_per_cell_group$cell_group = rownames(total_n_reads_per_cell_group)
  total_n_reads_per_cell_group_more_than_threshold = as.data.frame(total_n_reads_per_cell_group[total_n_reads_per_cell_group$total_n_reads_per_cell_group > min_n_reads_per_cluster,])

  for (i in c(2:ncol(cell_group_ges_cpm)))
  {
    cell_group_ges_cpm[, i] =  cell_group_ges[, i] / (total_n_reads_per_cell_group[i -
                                                                                     1, 1] + 1) * 1e6
  }




  cell_group_ges_hugo = cell_group_ges
  cell_group_ges_hugo[, 1] = feature.names[2]
  cell_group_ges_hugo_cpm = cell_group_ges_cpm
  cell_group_ges_hugo_cpm[, 1] = feature.names[2]

  filtered_cell_groups = total_n_reads_per_cell_group_more_than_threshold$cell_group  # Name of sample clusters that have total number of counts more than threshold.
  filtered_cell_group_ges_hugo = cell_group_ges_hugo[, c("SYMBOL", filtered_cell_groups)]
  filtered_cell_group_ges = cell_group_ges[, c("SYMBOL", filtered_cell_groups)]

  total_n_reads_per_cell_group_replicates = colSums(cell_group_ges_replicates[2:ncol(cell_group_ges_replicates)])
  total_n_reads_per_cell_group_replicates = as.data.frame(total_n_reads_per_cell_group_replicates)
  total_n_reads_per_cell_group_replicates$cell_group_replicate = rownames(total_n_reads_per_cell_group_replicates)
  total_n_reads_per_cell_group_replicates_more_than_threshold = as.data.frame(total_n_reads_per_cell_group_replicates[total_n_reads_per_cell_group_replicates$total_n_reads_per_cell_group > min_n_reads_per_cluster,])

  cell_group_ges_replicates = cell_group_ges_replicates

  filtered_cell_groups_replicates = total_n_reads_per_cell_group_replicates_more_than_threshold$cell_group  # Name of sample clusters replicates that have total number of counts more than threshold.
  filtered_cell_group_replicates_ges = cell_group_ges_replicates[, c("SYMBOL", filtered_cell_groups_replicates)]


  if (extract_cell_group)
  {
    prefix ="cell_group"
  }else
  {
    prefix ="sample"
  }


  write.table(
    cell_group_ges,
    file = paste0(outdir, "/", prefix, "_ges.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    cell_group_ges_cpm,
    file = paste0(outdir, "/", prefix, "_ges_cpm.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  write.table(
    cell_group_ges_hugo,
    file = paste0(outdir, "/", prefix, "_ges_translated.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    cell_group_ges_hugo_cpm,
    file = paste0(outdir, "/", prefix, "_ges_cpm_translated.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  write.table(
    filtered_cell_group_ges_hugo,
    file = paste0(outdir, "/", prefix, "_ges_translated.filtered.csv"),
    quote = FALSE,
    row.names = FALSE
  )


  write.table(
    total_n_reads_per_cell_group_more_than_threshold[, c(2, 1)],
    file = paste0(
      outdir,
      "/total_n_reads_per_", prefix, "_more_than_threshold.csv"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )



  write.table(
    total_n_reads_per_cell_group_replicates_more_than_threshold[, c(2, 1)],
    file = paste0(
      outdir,
      "/total_n_reads_per_", prefix, "_replicates_more_than_threshold.csv"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )





  if (extract_cell_group)
  {count_outdir = outdir
  }else
  {
    count_outdir = paste0(outdir, "/SAMPLES_COUNTS")
  }


  normalized_and_export_counts_data(counts = filtered_cell_group_replicates_ges,
                                    outdir = count_outdir,
                                    organism = organism)



  normalized_and_export_counts_data(
    counts =   filtered_cell_group_ges,
    outdir = paste0(count_outdir, "/COUNTS_NO_REPLICATES"),
    organism = organism
  )

  mpi.quit()
  quit()

}









