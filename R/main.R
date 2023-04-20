
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

start_mpi()


opt <- parse_args(OptionParser(option_list = option_list))

extract_ges(
  aggregate_out_dir  = opt$indir,
  min_n_reads_per_cluster = opt$min_n_reads_per_cluster,
  n_replicates = opt$n_replicates,
  outdir = opt$outdir,
  organism = opt$organism,
  clusters_file = opt$clusters_file
)


mpi.quit()
quit()

