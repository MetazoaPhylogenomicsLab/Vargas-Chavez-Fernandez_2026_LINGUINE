#' @title Auto-Standardize GFF Files
#'
#' @description Reads raw, messy GFF/GFF3 files from an input directory, dynamically extracts
#' the true gene IDs (hunting for 'ID=', 'Name=', 'gene_id=', or 'locus_tag='), and outputs
#' strictly standardized, LINGUINE-compatible GFF3 files where the attributes column is strictly 'ID=...'.
#' This eliminates the need for manual bash/awk formatting by the user.
#'
#' @param input_dir character. Path to the directory containing raw `.gff` or `.gff3` files.
#' @param output_dir character. Path to the directory where standardized GFFs will be saved.
#'
#' @return Invisible NULL. Writes standardized files to the output directory.
#' @export
standardize_linguine_inputs <- function(input_dir, output_dir) {
  if (!dir.exists(input_dir)) stop("Input directory does not exist.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  gff_files <- list.files(input_dir, pattern = "\\.gff3?$", full.names = TRUE)
  if (length(gff_files) == 0) {
    message("No GFF files found in the specified input directory.")
    return(invisible(NULL))
  }

  message(sprintf("Found %d GFF files. Standardizing...", length(gff_files)))

  for (file in gff_files) {
    species_name <- sub("\\.gff3?$", "", basename(file))
    message(sprintf("  Processing %s...", species_name))

    # Read the raw file
    raw_lines <- readLines(file, warn = FALSE)
    
    # Separate headers and data
    is_comment <- grepl("^#", raw_lines)
    comments <- raw_lines[is_comment]
    data_lines <- raw_lines[!is_comment]
    
    if (length(data_lines) == 0) next
    
    # Read as dataframe
    gff_df <- read.delim(text = data_lines, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
    
    # Ensure it has exactly 9 columns
    if (ncol(gff_df) < 9) {
      message(sprintf("    Warning: %s does not have 9 columns. Skipping.", basename(file)))
      next
    }
    
    colnames(gff_df) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
    
    # Filter for genes to reduce noise
    gff_df <- gff_df |> dplyr::filter(type == "gene")
    
    if (nrow(gff_df) == 0) {
      message(sprintf("    Warning: No 'gene' features found in %s.", basename(file)))
      next
    }
    
    # Vectorized regex extraction to find the true ID
    # Priority: ID > gene_id > Name > locus_tag > gene
    extracted_ids <- dplyr::coalesce(
      stringr::str_extract(gff_df$attributes, "(?<=ID=)[^;]+"),
      stringr::str_extract(gff_df$attributes, "(?<=gene_id=)[^;]+"),
      stringr::str_extract(gff_df$attributes, "(?<=Name=)[^;]+"),
      stringr::str_extract(gff_df$attributes, "(?<=locus_tag=)[^;]+"),
      stringr::str_extract(gff_df$attributes, "(?<=gene=)[^;]+")
    )
    
    if (any(is.na(extracted_ids))) {
      missing_count <- sum(is.na(extracted_ids))
      message(sprintf("    Warning: Could not parse an ID for %d genes. They will be skipped.", missing_count))
    }
    
    # Standardize the attributes column
    gff_df$attributes <- ifelse(!is.na(extracted_ids), paste0("ID=", extracted_ids), NA_character_)
    
    # Filter valid lines
    gff_df <- gff_df |> dplyr::filter(!is.na(attributes))
    
    # Write the standardized file
    output_file <- file.path(output_dir, paste0(species_name, ".gff"))
    
    writeLines(c("##gff-version 3", "# Standardized by LINGUINE"), output_file)
    write.table(gff_df, file = output_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  message("\nStandardization complete. Outputs saved to output directory.")
  return(invisible(NULL))
}