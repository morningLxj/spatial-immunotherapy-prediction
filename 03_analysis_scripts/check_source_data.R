#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grepl(file_arg, args)][1])
candidate_roots <- unique(c(
  getwd(),
  if (!is.na(script_path)) file.path(dirname(normalizePath(script_path, mustWork = FALSE)), "..") else NA
))
candidate_roots <- candidate_roots[!is.na(candidate_roots)]
has_data <- file.exists(file.path(candidate_roots, "data", "S1_source_data_and_machine_readable_tables.zip"))
if (!any(has_data)) {
  stop("Could not locate repository root containing data/S1_source_data_and_machine_readable_tables.zip.", call. = FALSE)
}
repo_root <- normalizePath(candidate_roots[which(has_data)[1]], mustWork = TRUE)

archive <- file.path(repo_root, "data", "S1_source_data_and_machine_readable_tables.zip")
if (!file.exists(archive)) {
  stop("Missing source-data archive: ", archive, call. = FALSE)
}

items <- utils::unzip(archive, list = TRUE)
release_dir <- file.path(repo_root, "release")
dir.create(release_dir, showWarnings = FALSE, recursive = TRUE)

out <- file.path(release_dir, "source_data_inventory.csv")
utils::write.csv(items, out, row.names = FALSE)

cat("Source-data archive found:\n")
cat("  ", archive, "\n", sep = "")
cat("Files in archive: ", nrow(items), "\n", sep = "")
cat("Inventory written:\n")
cat("  ", out, "\n", sep = "")
