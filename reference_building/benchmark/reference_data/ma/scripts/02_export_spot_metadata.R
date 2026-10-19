#!/usr/bin/env Rscript
# Exports per-spot metadata from Ma et al.'s Seurat object list (fetched by
# 01_fetch_al_dalahmah_zenodo_data.py); bayes_sample_key joins to 03_build_benchmark.py's
# RCTD table.

suppressPackageStartupMessages(library(SeuratObject))
source("utils/logging_extra.R")

data_dir <- "reference_building/benchmark/reference_data/ma/data"
infile_name <- file.path(data_dir, "original/processed_st_obj.rds")
outfile_name <- file.path(data_dir, "reformatted/spot_metadata_all_samples.csv.gz")

args <- commandArgs(trailingOnly = TRUE)
rds_path <- if (length(args) >= 1) args[1] else infile_name
output_csv <- if (length(args) >= 2) args[2] else outfile_name

log_run_command(file.path(data_dir, paste0(tools::file_path_sans_ext(basename(this_file())), ".log")))

obj_list <- readRDS(rds_path)

all_meta <- do.call(rbind, lapply(names(obj_list), function(sample_name) {
  md <- obj_list[[sample_name]]@meta.data
  data.frame(
    sample = sample_name,
    barcode = rownames(md),
    region = md$Region,
    condition = md$Condition,
    factor_1 = md$Factor_1,
    factor_5 = md$Factor_5,
    bayes_cluster = md$Bayes_Clusters,
    stringsAsFactors = FALSE
  )
}))
all_meta$bayes_sample_key <- paste0(all_meta$sample, "_", all_meta$bayes_cluster)

cat("Total spots:", nrow(all_meta), "\n")
print(table(all_meta$sample))

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
con <- gzfile(output_csv, "w")
write.csv(all_meta, con, row.names = FALSE)
close(con)
cat("Written", output_csv, "\n")
