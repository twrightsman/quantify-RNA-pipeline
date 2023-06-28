library(readr)
library(tximport)

tx2gene_path <- snakemake@input[["tx2gene"]]
tx2gene <- read_table(tx2gene_path, col_names = c("transcript", "gene"), col_types = "cc", progress = FALSE)

quants_tx_path <- snakemake@input[["quants_tx"]]
txi <- tximport(quants_tx_path, type = "salmon", tx2gene = tx2gene)
colnames(txi$abundance) <- "TPM"
colnames(txi$counts) <- "counts"
colnames(txi$length) <- "length"

abundance_and_counts <- data.frame(merge(txi$abundance, txi$counts, by = "row.names"), row.names = 1)
abundance_counts_length <- data.frame(merge(abundance_and_counts, txi$length, by = "row.names"), row.names = 1)
abundance_counts_length$gene <- rownames(abundance_counts_length)

output_path <- snakemake@output[["quants_gene"]]
write.table(abundance_counts_length[,c("gene", "TPM", "counts", "length")], file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
