#' @title Export Ancestral Genome to Standard Formats
#'
#' @description Translates the final LINGUINE Ancestor `.rds` object into universally readable 
#' `.tsv` (tabular gene order) and `.bed` formats for use in external tools like SynVisio.
#'
#' @param root_node character. The ID of the root node.
#' @param node_daughters character vector. The two daughters of the root node.
#' @param config list. The LINGUINE configuration list.
#'
#' @return Invisible NULL. Writes files to the results directory.
#' @export
export_ancestral_genome <- function(root_node, node_daughters, config) {
  results_dir <- config$paths$results
  
  # Identify the correct final rds file
  rds_file <- file.path(results_dir, paste0("ancestral_genome_", node_daughters[1], "_", node_daughters[2], ".rds"))
  if (!file.exists(rds_file)) {
    rds_file <- file.path(results_dir, paste0("ancestral_genome_", node_daughters[2], "_", node_daughters[1], ".rds"))
  }
  
  if (!file.exists(rds_file)) {
    message("Warning: Root ancestral state not found. Skipping BED/TSV export.")
    return(invisible(NULL))
  }
  
  message(sprintf("Exporting %s to TSV and BED formats...", basename(rds_file)))
  ancestral_df <- readRDS(rds_file)
  
  # Ensure necessary structure exists
  if (!all(c("linkage_group_name", "orthogroups_in_block") %in% names(ancestral_df))) {
    message("Warning: Ancestral object does not contain expected columns. Skipping export.")
    return(invisible(NULL))
  }
  
  # 1. Export to detailed TSV
  # Unnest orthogroups to show ordered content
  tsv_df <- ancestral_df |>
    dplyr::select(linkage_group_name, orthogroups_in_block) |>
    tidyr::unnest(orthogroups_in_block) |>
    dplyr::rename(Ancestral_LG = linkage_group_name, Orthogroup = orthogroups_in_block) |>
    dplyr::group_by(Ancestral_LG) |>
    dplyr::mutate(Relative_Order = dplyr::row_number()) |>
    dplyr::ungroup()
  
  tsv_file <- file.path(results_dir, paste0("ancestral_genome_", root_node, "_exported.tsv"))
  write.table(tsv_df, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 2. Export to pseudo-BED
  # Because an ancestor doesn't have literal base pairs, we use Orthogroup index as coordinates
  bed_df <- tsv_df |>
    dplyr::group_by(Ancestral_LG) |>
    dplyr::summarise(
      chromStart = min(Relative_Order) - 1, # 0-based for BED
      chromEnd = max(Relative_Order),
      name = paste0(unique(Ancestral_LG), "_block")
    ) |>
    dplyr::ungroup() |>
    dplyr::select(chrom = Ancestral_LG, chromStart, chromEnd, name)
  
  bed_file <- file.path(results_dir, paste0("ancestral_genome_", root_node, "_exported.bed"))
  write.table(bed_df, file = bed_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message("Export complete. Files saved in results directory.")
  return(invisible(NULL))
}