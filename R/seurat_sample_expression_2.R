
requireNamespace("Matrix", quietly = TRUE)
requireNamespace("optparse", quietly = TRUE)
requireNamespace("Rmpi", quietly = TRUE)
library(Matrix)
library(optparse)
library(Rmpi)

print("some update here")

print("Starting extraction single cell gene expresion profiles")

# run mpi

print("beginning mpi")
closed_slaves <- 0
ns <- mpi.universe.size() - 1
print(paste("Number of slaves:", ns))
if (ns == 0)
  #desktop
{
  ns = 1
}

n_current_slaves =  mpi.comm.size()
print(paste("Current number of unclosed slaves:",  n_current_slaves))
n_slaves_to_spawn = ns - n_current_slaves
if (n_slaves_to_spawn > 0)
{
  mpi.spawn.Rslaves(nslaves = n_slaves_to_spawn)

}

# In case R exits unexpectedly, have it automatically clean up
# resources taken up by Rmpi (slaves, memory, etc...)
.Last <- function() {
  if (is.loaded("mpi_initialize")) {
    if (mpi.comm.size(1) > 0) {
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}
# Tell all slaves to return a message identifying themselves
mpi.bcast.cmd(id <- mpi.comm.rank())
mpi.bcast.cmd(ns <- mpi.comm.size())
mpi.bcast.cmd(host <- mpi.get.processor.name())
mpi.bcast.cmd(library("Matrix"))


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
extract_ges = function(aggregate_out_dir  = ".",
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


}




read_ges = function(aggregate_out_dir)
{
  matrix_dir = paste0(aggregate_out_dir , "/filtered_feature_bc_matrix")

  barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "/features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")
  print("Reading filtered_feature_bc_matrix file, please wait")
  ges <- readMM(file = matrix.path)
  print("Gene expression profile matrix has been  read")
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(ges) = barcode.names$V1
  rownames(ges) = feature.names$V1
  return(list(ges, feature.names, barcode.names))
}


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

slave_control = function()
  #Manage the slave processing
{
  # sent:  1=ready_for_task, 2=done_task, 3=exiting received:  1=task, 2=done_tasks
  ges_list = read_ges(aggregate_out_dir) # Each slave read the data from the hard disk as rmpi cannot pass huge dataset over the network. Error Long Vector not supported
  ges = ges_list[[1]]

  junk <- 0
  done <- 0 # Flag if all tasks are done
  while (done != 1) {
    # Signal being ready to receive a new task
    mpi.send.Robj(junk, 0, 1) # Send to Master: ready for task
    task <-
      mpi.recv.Robj(mpi.any.source(), mpi.any.tag()) # Receive a task
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]
    if (tag == 1) {
      #task exits
      print("processing task")
      results = extract_cluster_ges(task, ges)
      print("done for this task, sending back result")
      # Send a results message back to the master, tag: done_task
      mpi.send.Robj(results, 0, 2)
    }
    else if (tag == 2) {
      done <- 1
    } # We'll just ignore any unknown messages
  }
  mpi.send.Robj(junk, 0, 3)
  # Exiting received

}

normalized_and_export_counts_data = function(counts = NULL,
                                             outdir = ".",
                                             organism = "human")
{
  # From counts table, with columns "SYMBOL", samples where SYMBOL is ENSEMBL gene ID, normalize counts and save normalized rpkm, tpm, cpm, count data into the CONSOLIDATED_COUNTS folder

  normalized_data = get_normalized_data(counts,
                                        organism = organism)
  CONSOLIDATED_COUNTS_dir  = paste0(outdir, "/CONSOLIDATED_COUNTS/")
  print(paste("outdir:", outdir))
  dir.create(outdir)
  dir.create(CONSOLIDATED_COUNTS_dir)
  for (count_type in c("count", "cpm", "rpkm", "tpm"))
  {
    print(count_type)
    dir.create(paste0(CONSOLIDATED_COUNTS_dir, count_type))
    dir.create(paste0(CONSOLIDATED_COUNTS_dir, count_type, "/gene_level"))

    export_datasets = normalized_data[[count_type]]
    export_dataset = export_datasets[["ensemble_id"]]
    export_dataset_translated = export_datasets[["normal_id"]]
    write.table(
      export_dataset,
      file = paste0(
        outdir,
        "/CONSOLIDATED_COUNTS/",
        count_type,
        "/gene_level/consolidated_",
        count_type,
        "_table.csv"
      ),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
    write.table(
      export_dataset_translated,
      file = paste0(
        outdir,
        "/CONSOLIDATED_COUNTS/",
        count_type,
        "/gene_level/consolidated_",
        count_type,
        "_table_translated.csv"
      ),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
  }

}



get_normalized_data = function(counts = NULL,
                               organism = "human")
  # Convert counts data into cpm, rpkm, tpm data
  # Input: counts is dataframe where column 1 contain Ensemble ID , other columns are samples, rows are number of reads
  # Output: list of normalized data, where for each count type in c(count, cpm, rpkm, tpm) there are two types of normalized data: by Ensemble ID and by HUGO ID. Have not tested for mouse yet.
{
  library(biomaRt)
  library(edgeR)
  print(counts[1:5, 1:5])
  ensembl_list <- counts[, 1]
  if (organism == "human")
  {
    mart_dataset = mart_dataset = "hsapiens_gene_ensembl"
  }
  if (organism == "mouse")
  {
    mart_dataset = mart_dataset = "mmusculus_gene_ensembl"
  }
  if (!exists("mart"))
  {
    mart <-
      useMart(
        "ensembl",
        dataset = mart_dataset,
        host = "https://www.ensembl.org"
      )
  }

  if (organism == "human")
  {gene_name_field = "hgnc_symbol"}else
  {gene_name_field = "external_gene_name"}

  gene_coords = getBM(
    attributes = c(
      gene_name_field,
      "ensembl_gene_id",
      "start_position",
      "end_position"
    ),
    filters = "ensembl_gene_id",
    values = ensembl_list,
    mart = mart

  )
  gene_coords$Length = (gene_coords$end_position - gene_coords$start_position) /
    1000 # Length in kilobase

  # The gene coord  names may be not in the order as the input gene list, some genes may missing.

  # Replace empty or repeated hgnc symbol with ensebmle gene id as missing or repeated hgnc symbol will interfere with downstream task that use
  for (i in c(1:nrow(gene_coords)))
  {

    if (gene_coords[i, gene_name_field] == "")
    {
      gene_coords[i, gene_name_field] = gene_coords[i, "ensembl_gene_id"]
    }

  }
  repeated_hgnc_symbols = table(gene_coords[, gene_name_field])
  repeated_hgnc_symbols =  as.data.frame(repeated_hgnc_symbols)
  repeated_hgnc_symbols = repeated_hgnc_symbols[repeated_hgnc_symbols$Freq >
                                                  1, 1]
  repeated_hgnc_symbols  = as.vector(repeated_hgnc_symbols)

  x = counts[gene_coords$ensembl_gene_id, 2:ncol(counts)]
  x_rpkm <- rpkm(x, gene_coords$Length * 1000)
  x_cpm = cpm(x)
  gene_coords[gene_coords$hgnc_symbol %in% repeated_hgnc_symbols, gene_name_field] = gene_coords[gene_coords$hgnc_symbol %in% repeated_hgnc_symbols, "ensembl_gene_id"]


  tpm_reads_per_length = x / gene_coords$Length
  reads_per_length_sum = colSums(tpm_reads_per_length) / 10 ^ 6

  x_tpm = sweep(tpm_reads_per_length, 2, reads_per_length_sum, '/')


  normalized_data = list()
  for (count_type in c("count", "cpm", "rpkm", "tpm"))
  {
    switch (
      count_type,
      "count" = {
        export_dataset = x
      },
      "cpm" = {
        export_dataset = x_cpm
      },
      "rpkm" = {
        export_dataset = x_rpkm
      },
      "tpm" = {
        export_dataset = x_tpm
      }
    )
    export_dataset = as.data.frame(export_dataset)
    export_dataset$SYMBOL = gene_coords$ensembl_gene_id
    ncols = ncol(export_dataset)
    export_dataset = export_dataset[, c(ncols, 1:ncols - 1)]
    export_dataset_translated = export_dataset
    export_dataset_translated$SYMBOL = toupper(gene_coords[, gene_name_field])
    print(paste("done with", count_type))


    normalized_data[[count_type]] = list(ensemble_id = export_dataset, normal_id = export_dataset_translated)

  }

  return(normalized_data)




}




option_list <- list(
  make_option(
    c("-d", "--indir"),
    type = "character",
    default = ".",
    help = "input directory for the script. The outs Folder of aggregated sample, parent of filtered_feature_bc_matrix"
  ),
  make_option(
    c("-m", "--min_n_reads_per_cluster"),
    type = "integer",
    default = 1000,
    help = "Minimal number of reads per cluster to filter out low read clusters"
  ),
  make_option(
    c("-r", "--n_replicates"),
    type = "integer",
    default = 5,
    help = "number of replicates to sample  when creating count sample cluster matrix"
  ),

  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "NETZEN_analysis",
    help = "output directory"
  ),
  make_option(
    c("--organism"),
    type = "character",
    default = "human",
    help = "output directory"
  ),
  make_option(
    c("--clusters_file"),
    type = "character",
    default = "./ClusterAnnotations.csv",
    help = "file with cluster annotation, by default using 10X cluster but, can be used to provided customised clusters such as from Seurat. Table with 2 columns Barcode, Cluster with header."
  )
)


# opt <- parse_args(OptionParser(option_list = option_list))
#
# extract_ges(
#   aggregate_out_dir  = opt$indir,
#   min_n_reads_per_cluster = opt$min_n_reads_per_cluster,
#   n_replicates = opt$n_replicates,
#   outdir = opt$outdir,
#   organism = opt$organism,
#   clusters_file = opt$clusters_file
# )
#
#
# mpi.quit()
# quit()
