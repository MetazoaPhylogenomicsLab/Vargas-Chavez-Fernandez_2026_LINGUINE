#' @title Process Extant Genome Information
#'
#' @description Extracts chromosomal coordinates from FASTA files and filters gene
#' models from GFF files for extant tip species. Skips processing if intermediate
#' files are already detected.
#'
#' @param species_name Character. Target species identifier.
#' @param is_inode Logical. True if the target is an internal node (skips GFF processing).
#' @param config A \code{linguine_config} object created by \code{create_linguine_config()}.
#'
#' @return A list containing the file paths to the processed chromosome sizes and gene dataframes.
#' @export
process_extant_genome <- function(species_name, is_inode = FALSE, config) {

  # 1. Validate the configuration object
  if (!inherits(config, "linguine_config")) {
    stop("Error: The 'config' argument must be a valid linguine_config object.")
  }

  # 2. Dynamically construct paths using the config object
  fasta_path <- file.path(config$paths$raw_data, paste0(species_name, ".fna"))
  chr_sizes_path <- file.path(config$paths$processed_data, paste0(species_name, "_chromosome_sizes.rds"))
  ref_genes_df_path <- file.path(config$paths$processed_data, paste0(species_name, "_genes_df.rds"))

  # 3. Coordinate Extraction (FASTA Parsing)
  if (file.exists(chr_sizes_path)) {
    message("Info: Processed chromosome metrics for ", species_name, " detected. Skipping extraction.")
  } else {
    if (file.exists(fasta_path)) {
      extract_chr_sizes(
        ref_fasta_path = fasta_path,
        output_chr_sizes_path = chr_sizes_path,
        min_length = config$min_chromosome_length_bp
      )
    } else {
      stop("Error: Required FASTA file for ", species_name, " not found at ", fasta_path)
    }
  }

  # 4. Gene Model Pre-processing (GFF Filtering)
  if (file.exists(ref_genes_df_path)) {
    message("Info: Processed gene models for ", species_name, " detected. Skipping filtration.")
  } else {
    if (!is_inode) {
      invisible(filter_gff(species_name, config))
    }
  }

  # 5. Return the locations of the generated data
  return(list(
    chr_sizes_path = chr_sizes_path,
    genes_df_path = if (!is_inode) ref_genes_df_path else NULL
  ))
}
