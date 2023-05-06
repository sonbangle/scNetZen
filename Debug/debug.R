


aggregation_samples_file = "~/Downloads/count/aggregation.csv"
aggregate_out_dir = "~/Downloads/count"
barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
clusters_file = "~/Downloads/count/analysis/clustering/graphclust/clusters.csv"
organism = "mouse"
outdir =  "~/Downloads/count/GES"
barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
tasks = get_task_for_samples_clusters(aggregation_samples_file = aggregation_samples_file , clusters_file = clusters_file,  barcode.names = barcode.names)[["tasks"]]
saveRDS(tasks, "./Debug/tasks.RDS" )


prefix = "sample_cluster"
method = "singlecore"
n_replicates = 5
out = get_cell_group_gene_expression_from_tasks(tasks,
                                                aggregate_out_dir,
                                                outdir,
                                                prefix,
                                                organism,
                                                method = method,
                                                close_slaves = close_slaves,
                                                n_replicates = n_replicates,
                                                min_n_reads_per_cell_group = min_n_reads_per_cell_group)


