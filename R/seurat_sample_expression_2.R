
#' Title
#'
#' @param aggregate_out_dir
#' @param min_n_reads_per_cluster
#' @param n_replicates
#' @param outdir
#' @param organism
#' @param clusters_file
#' @param extract_sample_cluster
#'
#' @return
#' @export
#'
#' @examples
extract_ges_mpi = function(aggregate_out_dir  = ".",
                       min_n_reads_per_cluster = 20000,
                       n_replicates = 5,
                       outdir = "NETZEN_analysis",
                       organism = "human",
                       clusters_file = "./analysis/clustering/graphclust/clusters.csv",
                       extract_sample_cluster = FALSE
)

{
  # This function extracts the percentage and absolute counts of cells in each cluster and sammple, raw number of reads per cluster per sample  and normalized count per million(cpm) per each cluster per each sample
  # Input aggreated_out_dir: "outs" Folder of aggreated sample, parent of filtered_feature_bc_matrix
  # Clusters_file: file with cluster annotation, by default using 10X cluster but, can be used to provided customised clusters such as from Seurat. Table with 2 columns Barcode, Cluster with header.
  #extract_sample_cluster: TRUE: extract  cluster expression for each sample, cluster. FALSE: extract only sample expression)
  ns = start_mpi()
  dir.create(outdir)

  # Broadcasting at the beginning so that slaves can read the gene expression file early.

  print("Beginning broadcasting functions and libraries to slave")

  mpi.bcast.Robj2slave(extract_cluster_ges)
  mpi.bcast.Robj2slave(slave_control)
  mpi.bcast.Robj2slave(read_ges)
  mpi.bcast.Robj2slave(aggregate_out_dir)
  mpi.bcast.Robj2slave(n_replicates)


  print("done broadcasting libraries and functions to slaves")

  mpi.bcast.cmd(slave_control()) # Call the function in all the slaves

  ges_list = read_ges(aggregate_out_dir)
  ges = ges_list[[1]]
  feature.names = ges_list[[2]]
  barcode.names = ges_list[[3]]


  barcode.names$samples = substr(barcode.names$V1, 18, nchar(barcode.names$V1))

  if (extract_sample_cluster)
  {
    print(clusters_file)

    clusters <- read.csv(clusters_file, stringsAsFactors = FALSE)
    clusters = as.data.frame(clusters)
    colnames(clusters) = c("Barcode", "Cluster")
    clusters = as.data.frame(clusters)
    clusters[clusters$Cluster == "", 2] = "Unassigned"
    # Remove space in the Cluster names
    clusters$Cluster = gsub(" ", "_", clusters$Cluster)
    uniqclusters = unique(clusters$Cluster)
  }

  aggregation_samples <-
    read.csv(paste0(aggregate_out_dir, "/aggregation.csv"),
             stringsAsFactors = FALSE)
  colnames(aggregation_samples)[1] = "library_id" # Change the first column name to library_id as it can be sample_id initially.
  samples = aggregation_samples$library_id
  if (extract_sample_cluster)
  {
    n_clusters =  length(uniqclusters)
  } else
  {
    n_clusters = 1

  }
  n_samples = length(samples)
  print(paste("Total number of samples", n_samples))
  if (n_samples == 0)
  {
    print(paste("aggregate_outdir:", aggregate_out_dir))
    print(paste(
      "aggreation_samples_file:",
      paste0(aggregate_out_dir, "/aggregation.csv")
    ))
    print(paste("samples:", samples))
    print("Error, there is no sample to process")
    stop()

  }


  if (extract_sample_cluster)
  {
    sample_cluster_matrix = matrix(data = 0,
                                   nrow = n_clusters,
                                   ncol = n_samples)
    sample_cluster_matrix = data.frame(sample_cluster_matrix)
    colnames(sample_cluster_matrix) = samples
    #rownames(sample_cluster_matrix) = paste0("cluster_", c(1:length(uniqclusters)))
    uniqclusters = as.vector(uniqclusters)
    rownames(sample_cluster_matrix) = uniqclusters

  }




  #####################################################################################
  # Create task list
  tasks = list()
  for (i in c(1:n_samples))
  {
    sample = samples[i]
    sample_barcodes = barcode.names[barcode.names$samples == i, 1]

    if (extract_sample_cluster)
    {
      cluster_assignments  = clusters[clusters$Barcode %in% sample_barcodes,]
      for (j  in c(1:n_clusters))
      {
        cluster_barcodes = as.vector(cluster_assignments[cluster_assignments$Cluster == uniqclusters[j], "Barcode"])
        n_barcodes = length(cluster_barcodes)
        sample_cluster_matrix[j, i] = n_barcodes
        tasks[[(i - 1) * n_clusters + j]] = list(
          sample = sample,
          cluster = uniqclusters[j],
          cluster_barcodes = cluster_barcodes,
          i = i,
          j = j,
          n_replicates = n_replicates
        )
      }
    }
    else
    {
      tasks[[i]] = list(
        sample = sample,
        cluster = sample,
        cluster_barcodes = sample_barcodes,
        i = i,
        j = 1,
        n_replicates = n_replicates
      )
    }


  }



  if (extract_sample_cluster)
  {
    # Write proportion table for sample clusters
    total_n_cells = colSums(sample_cluster_matrix)

    sample_cluster_matrix_percentage = sample_cluster_matrix
    for (i in c(1:ncol(sample_cluster_matrix_percentage)))
    {
      sample_cluster_matrix_percentage[, i] =  sample_cluster_matrix[, i] / total_n_cells[i] * 100
    }


    sample_cluster_matrix_assigned_only =  sample_cluster_matrix[rownames(sample_cluster_matrix) != "Unassigned",]
    total_n_cells_assigned_only = colSums(sample_cluster_matrix)
    sample_cluster_matrix_percentage_assigned_only = sample_cluster_matrix_assigned_only
    print("sample_cluster_matrix_percentage_assigned_only ")
    print(sample_cluster_matrix_percentage_assigned_only)
    print(paste(
      "col(sample_cluster_matrix_percentage_assigned_only):",
      ncol(sample_cluster_matrix_percentage_assigned_only)
    ))
    for (i in c(1:ncol(sample_cluster_matrix_percentage_assigned_only)))
    {
      sample_cluster_matrix_percentage_assigned_only[, i] =  sample_cluster_matrix_assigned_only[, i] / total_n_cells_assigned_only[i] * 100
    }


    write.table(
      sample_cluster_matrix_percentage,
      file = paste0(outdir, "/sample_cluster_matrix_percentage.csv"),
      quote = FALSE,
      sep = "\t"
    )
    write.table(
      sample_cluster_matrix,
      file = paste0(outdir, "/sample_cluster_matrix.csv"),
      quote = FALSE,
      sep = "\t"
    )

    write.table(
      sample_cluster_matrix_percentage_assigned_only,
      file = paste0(
        outdir,
        "/sample_cluster_matrix_percentage_assigned_only.csv"
      ),
      quote = FALSE,
      sep = "\t"
    )
    write.table(
      sample_cluster_matrix_assigned_only,
      file = paste0(outdir, "/sample_cluster_matrix_assigned_only.csv"),
      quote = FALSE,
      sep = "\t"
    )

  }


  ##############################################################Making gene expression profiles for NETZEN



  sample_cluster_ges = matrix(data = 0,
                              nrow = nrow(ges) ,
                              ncol = n_clusters * n_samples + 1)
  sample_cluster_ges = data.frame(sample_cluster_ges)
  rownames(sample_cluster_ges) = rownames(ges)
  sample_cluster_ges[, 1] = rownames(sample_cluster_ges)
  colnames(sample_cluster_ges)[1] = "SYMBOL"
  ges_col_index = 2


  sample_cluster_ges_replicates = matrix(
    data = 0,
    nrow = nrow(ges) ,
    ncol = n_clusters * n_samples * n_replicates + 1
  )
  sample_cluster_ges_replicates = data.frame(sample_cluster_ges_replicates)
  rownames(sample_cluster_ges_replicates) = rownames(ges)
  sample_cluster_ges_replicates[, 1] = rownames(sample_cluster_ges_replicates)
  colnames(sample_cluster_ges_replicates)[1] = "SYMBOL"
  ges_col_index = 2
  ges_col_index_replicates = 2




  # Create data structures to store the results
  cluster_ges_sum.result = list()
  cluster_ges_sum_replicates.result = list()

  ###########################
  # send tasks to slaves



  results = list()
  done_tasks = 0
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
      if (length(tasks) > 0) {
        # Send a task, and then remove it from the task list
        mpi.send.Robj(tasks[[1]], slave_id, 1)

        tasks[[1]] <- NULL
      }
      else {
        mpi.send.Robj(junk, slave_id, 2)
      }
    } else if (tag == 2) {
      # Do something with the results. Store in the data structure

      sample = message[["sample"]]
      cluster = message[["cluster"]]
      cluster_ges_sum = message[["cluster_ges_sum"]]
      cluster_ges_sum_replicates = message[["cluster_ges_sum_replicates"]]
      sample_index = message[["sample_index"]]
      cluster_index = message[["cluster_index"]]
      print(sample)
      print(cluster)
      print(paste("sample index", sample_index))
      print(paste("cluster index", cluster_index))

      pos =   (sample_index - 1) * n_clusters + cluster_index + 1

      colnames(sample_cluster_ges)[pos] = paste0(sample, "^", cluster)

      start_pos = (sample_index - 1) * n_clusters * n_replicates + (cluster_index -
                                                                      1) * n_replicates  + 2
      end_pos = start_pos + n_replicates - 1

      colnames(sample_cluster_ges_replicates)[start_pos:end_pos] = paste0(sample, "^", cluster, "^rep_", c(1:n_replicates))

      if (message[["is_null_result"]] == 0)
        # If result is not null, then add sample_cluster column to a data frame
      {
        sample_cluster_ges[, pos] = cluster_ges_sum
        sample_cluster_ges_replicates[, start_pos:end_pos] = cluster_ges_sum_replicates
      }

    } else if (tag == 3) {
      closed_slaves <- closed_slaves + 1
    }
  }



  # Tell all slaves to close down, and exit the program
  mpi.close.Rslaves(dellog = FALSE)

  print("Done with mpi")

  ###############################################################







  total_n_reads_per_sample_cluster = colSums(sample_cluster_ges[2:ncol(sample_cluster_ges)])
  total_n_reads_per_sample_cluster = as.data.frame(total_n_reads_per_sample_cluster)
  sample_cluster_ges_cpm = sample_cluster_ges
  total_n_reads_per_sample_cluster$sample_cluster = rownames(total_n_reads_per_sample_cluster)
  total_n_reads_per_sample_cluster_more_than_threshold = as.data.frame(total_n_reads_per_sample_cluster[total_n_reads_per_sample_cluster$total_n_reads_per_sample_cluster > min_n_reads_per_cluster,])

  for (i in c(2:ncol(sample_cluster_ges_cpm)))
  {
    sample_cluster_ges_cpm[, i] =  sample_cluster_ges[, i] / (total_n_reads_per_sample_cluster[i -
                                                                                                 1, 1] + 1) * 1e6
  }




  sample_cluster_ges_hugo = sample_cluster_ges
  sample_cluster_ges_hugo[, 1] = feature.names[2]
  sample_cluster_ges_hugo_cpm = sample_cluster_ges_cpm
  sample_cluster_ges_hugo_cpm[, 1] = feature.names[2]

  filtered_sample_clusters = total_n_reads_per_sample_cluster_more_than_threshold$sample_cluster  # Name of sample clusters that have total number of counts more than threshold.
  filtered_sample_cluster_ges_hugo = sample_cluster_ges_hugo[, c("SYMBOL", filtered_sample_clusters)]
  filtered_sample_cluster_ges = sample_cluster_ges[, c("SYMBOL", filtered_sample_clusters)]

  total_n_reads_per_sample_cluster_replicates = colSums(sample_cluster_ges_replicates[2:ncol(sample_cluster_ges_replicates)])
  total_n_reads_per_sample_cluster_replicates = as.data.frame(total_n_reads_per_sample_cluster_replicates)
  total_n_reads_per_sample_cluster_replicates$sample_cluster_replicate = rownames(total_n_reads_per_sample_cluster_replicates)
  total_n_reads_per_sample_cluster_replicates_more_than_threshold = as.data.frame(total_n_reads_per_sample_cluster_replicates[total_n_reads_per_sample_cluster_replicates$total_n_reads_per_sample_cluster > min_n_reads_per_cluster,])

  sample_cluster_ges_replicates = sample_cluster_ges_replicates

  filtered_sample_clusters_replicates = total_n_reads_per_sample_cluster_replicates_more_than_threshold$sample_cluster  # Name of sample clusters replicates that have total number of counts more than threshold.
  filtered_sample_cluster_replicates_ges = sample_cluster_ges_replicates[, c("SYMBOL", filtered_sample_clusters_replicates)]


  if (extract_sample_cluster)
  {
    prefix ="sample_cluster"
  }else
  {
    prefix ="sample"
  }


  write.table(
    sample_cluster_ges,
    file = paste0(outdir, "/", prefix, "_ges.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    sample_cluster_ges_cpm,
    file = paste0(outdir, "/", prefix, "_ges_cpm.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  write.table(
    sample_cluster_ges_hugo,
    file = paste0(outdir, "/", prefix, "_ges_translated.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  write.table(
    sample_cluster_ges_hugo_cpm,
    file = paste0(outdir, "/", prefix, "_ges_cpm_translated.csv"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  write.table(
    filtered_sample_cluster_ges_hugo,
    file = paste0(outdir, "/", prefix, "_ges_translated.filtered.csv"),
    quote = FALSE,
    row.names = FALSE
  )


  write.table(
    total_n_reads_per_sample_cluster_more_than_threshold[, c(2, 1)],
    file = paste0(
      outdir,
      "/total_n_reads_per_", prefix, "_more_than_threshold.csv"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )



  write.table(
    total_n_reads_per_sample_cluster_replicates_more_than_threshold[, c(2, 1)],
    file = paste0(
      outdir,
      "/total_n_reads_per_", prefix, "_replicates_more_than_threshold.csv"
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )





  if (extract_sample_cluster)
  {count_outdir = outdir
  }else
  {
    count_outdir = paste0(outdir, "/SAMPLES_COUNTS")
  }


  normalized_and_export_counts_data(counts = filtered_sample_cluster_replicates_ges,
                                    outdir = count_outdir,
                                    organism = organism)



  normalized_and_export_counts_data(
    counts =   filtered_sample_cluster_ges,
    outdir = paste0(count_outdir, "/COUNTS_NO_REPLICATES"),
    organism = organism
  )

  mpi.quit()
  quit()

}











