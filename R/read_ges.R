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
