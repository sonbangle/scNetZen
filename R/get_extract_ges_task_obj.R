
get_extract_ges_tasks_obj = function(outfile="/home/sonle/Downloads/scNetZen/tests/extract_ges_obj.RDS",
                                     aggregation_samples_file = "~/Downloads/count/aggregation.csv",
                                     aggregate_out_dir = "~/Downloads/count",
                                     method = "multicore"
                                     )
{
barcode.names = get_barcode_names(aggregate_out_dir = aggregate_out_dir)
tasks = get_task_from_aggregation_samples(aggregation_samples_file = aggregation_samples_file, barcode.names = barcode.names)
ges<-read_ges(aggregate_out_dir)[[1]]
out = extract_ges_parallel(tasks, method = method)
saveRDS(out, file = outfile)
invisible(out)
return(out)
}
