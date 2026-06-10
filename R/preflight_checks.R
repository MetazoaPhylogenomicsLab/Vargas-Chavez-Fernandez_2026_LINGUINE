#' @title Pre-Flight Checks for LINGUINE Pipeline
#'
#' @description Validates inputs before running the main LINGUINE pipeline.
#' Checks if directories exist, validates GFF formatting (ensuring 'ID=' is present),
#' and verifies that tree tips match OrthoFinder headers.
#'
#' @param config A `linguine_config` object.
#' @param speciesTree An `ape::phylo` object representing the species tree.
#'
#' @return Invisible TRUE if all checks pass. Halts execution with descriptive error if checks fail.
#' @export
run_preflight_checks <- function(config, speciesTree) {
  message("\n=======================================================")
  message("          Running Pre-flight Checks                    ")
  message("=======================================================")

  tips <- speciesTree$tip.label
  
  # 1. Directory and File Checks
  raw_dir <- config$paths$raw_data
  if (!dir.exists(raw_dir)) {
    stop(sprintf("Error: Raw data directory not found at: %s", raw_dir))
  }
  
  orthofinder_file <- file.path(raw_dir, config$orthologs$filename)
  if (!file.exists(orthofinder_file)) {
    stop(sprintf("Error: OrthoFinder file not found at: %s", orthofinder_file))
  }
  
  # 2. GFF Checks for all species
  message("Verifying GFF files for all species in the tree...")
  for (sp in tips) {
    gff_path <- file.path(raw_dir, config$prefixes$gff, paste0(sp, config$suffixes$gff))
    if (!file.exists(gff_path)) {
      stop(sprintf("Error: Missing GFF file for species '%s' at: %s", sp, gff_path))
    }
    
    # Peek at the first few lines to verify 'ID=' usage
    first_lines <- readLines(gff_path, n = 500)
    # Filter out comments
    data_lines <- first_lines[!grepl("^#", first_lines)]
    if (length(data_lines) > 0) {
      has_id <- any(grepl("ID=", data_lines))
      if (!has_id) {
        stop(sprintf("Error: GFF for species '%s' does not appear to use 'ID=' in the attributes column. LINGUINE requires standard 'ID=' tags.", sp))
      }
    }
  }
  
  # 3. OrthoFinder Header Verification
  message("Verifying OrthoFinder headers match tree tips...")
  # Read just the header
  ortho_header <- strsplit(readLines(orthofinder_file, n = 1), "\t")[[1]]
  
  missing_tips <- setdiff(tips, ortho_header)
  if (length(missing_tips) > 0) {
    stop(sprintf(
      "Error: The following species from the phylogenetic tree were NOT found in the OrthoFinder file headers:\n%s\nEnsure tree tip labels exactly match OrthoFinder columns.",
      paste(missing_tips, collapse = ", ")
    ))
  }
  
  message("Pre-flight checks passed successfully.\n")
  return(invisible(TRUE))
}