args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: prepare_cll_rds_export.R INPUT_RDS OUTPUT_DIRECTORY")
}

x <- readRDS(args[[1]])
out_dir <- args[[2]]
required <- c("mRNA", "Methylation", "Drugs", "Mutations")
missing_views <- setdiff(required, names(x))
if (length(missing_views)) {
    stop(paste("Missing RDS views:", paste(missing_views, collapse = ", ")))
}

for (view in required) {
    matrix <- x[[view]]
    frame <- data.frame(feature = rownames(matrix), matrix, check.names = FALSE)
    output <- gzfile(file.path(out_dir, paste0(view, ".tsv.gz")), "wb")
    write.table(
        frame,
        output,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        na = "NA"
    )
    close(output)
}
