library(Matrix)

read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # the matrix read has cells in rows
  genes <- readLines(file(paste0(dir, "/", "genes.txt")))
  rownames(m) <- genes
  return(m)
}

matrix <- read_count_output("./quant_unfiltered", name = "matrix.abundance.gene.tpm")
sub<-matrix[,1]
df <- as.data.frame(sub)
write.table(df,file="gene_tpm_matrix.txt",col.names=F,row.names=T,sep='\t')