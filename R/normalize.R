
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
