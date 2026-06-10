# ==============================================================================
# Script: R/utils.R
# Purpose: Core utility functions for the Ancestral Genome Reconstruction Pipeline.
# ==============================================================================

#' @title Retrieve Node Label
#' @description Extracts the standardized label for a specific tip or internal node
#' within an 'ape' phylogenetic tree structure.
#'
#' @param index integer. The numeric index of the tip or node.
#' @param tree phylo. The phylogenetic tree object.
#'
#' @return character. The assigned tip/node label, or a generated "INode_X" string.
#' @noRd
get_label <- function(index, tree) {
  num_tips <- length(tree$tip.label)

  if (index <= num_tips) {
    return(tree$tip.label[index])
  } else {
    if (is.null(tree$node.label)) {
      return(paste0("INode_", index))
    } else {
      label_index <- index - num_tips
      return(tree$node.label[label_index])
    }
  }
}

#' @title Extract Chromosome Sizes from FASTA
#' @description Parses a reference genome FASTA file, calculates sequence lengths,
#' and filters out short scaffolds based on a global minimum length threshold.
#'
#' @param ref_fasta_path character. Path to the input FASTA file (.fna).
#' @param output_chr_sizes_path character. Destination path for the RDS output.
#' @param min_length numeric. Minimum chromosome length filter in base pairs.
#'
#' @return NULL. Saves a dataframe to disk containing `ref_chromosome` and `chromosome_length_bp`.
#' @noRd
extract_chr_sizes <- function(ref_fasta_path, output_chr_sizes_path, min_length = 0) {
  if (!file.exists(ref_fasta_path)) {
    stop(paste0("Error: Reference FASTA file not found at ", ref_fasta_path))
  }

  message("--- Extracting Reference Genome Chromosome Information ---")
  message(paste0("Loading reference genome FASTA from: ", ref_fasta_path))

  ref_genome <- Biostrings::readDNAStringSet(ref_fasta_path)

  ref_chromosome_sizes_df <- dplyr::tibble(
    ref_chromosome = names(ref_genome),
    chromosome_length_bp = as.numeric(Biostrings::width(ref_genome))
  )

  ref_chromosome_sizes_df <- ref_chromosome_sizes_df |>
    dplyr::arrange(gtools::mixedsort(ref_chromosome)) |>
    dplyr::filter(chromosome_length_bp >= min_length)

  message("\nTotal chromosomes extracted (filtered > ", min_length, " bp): ", nrow(ref_chromosome_sizes_df))

  saveRDS(ref_chromosome_sizes_df, output_chr_sizes_path)
  message("\nChromosome sizes saved to: ", output_chr_sizes_path)
  message("\nChromosome information extraction complete.")
}

#' @title Parse GFF3 Attributes
#' @description Extracts the primary unique identifier from a standard GFF3 attributes string.
#'
#' @param attributes_string character. The unparsed string from the 9th column of a GFF3.
#'
#' @return character. The extracted ID (matching 'ID=', 'gene=', 'Name=', or 'locus_tag='),
#' or NA if unparseable.
#' @noRd
parse_gff_attributes <- function(attributes_string) {
  id_match <- stringr::str_match(attributes_string, "ID=([^;]+)")
  if (!is.na(id_match[2])) return(id_match[2])

  gene_match <- stringr::str_match(attributes_string, "gene=([^;]+)")
  if (!is.na(gene_match[2])) return(gene_match[2])

  name_match <- stringr::str_match(attributes_string, "Name=([^;]+)")
  if (!is.na(name_match[2])) return(name_match[2])

  locus_match <- stringr::str_match(attributes_string, "locus_tag=([^;]+)")
  if (!is.na(locus_match[2])) return(locus_match[2])

  if (!stringr::str_detect(attributes_string, "=") && nchar(attributes_string) > 0) {
    return(attributes_string)
  }
  return(NA_character_)
}

#' @title Extract Gene Features from GFF3
#' @description Parses a generic GFF file, filters for 'gene' feature types,
#' and extracts standardized genomic coordinates.
#'
#' @param filepath character. Path to the target GFF/GFF3 file.
#' @param species_prefix character. Identifier prefix for logging output.
#'
#' @return tibble. A deduplicated dataframe of gene models (`gene_id`, `chromosome`, `start`, `end`).
#' @noRd
read_gff_genes <- function(filepath, species_prefix) {
  message(paste0("Reading GFF file: ", filepath))

  gff_data <- read.delim(filepath, header = FALSE, comment.char = "#", sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

  gene_features <- gff_data |>
    dplyr::filter(type == "gene") |>
    dplyr::mutate(
      gene_id = sapply(attributes, parse_gff_attributes),
      chromosome = seqid,
      start = as.integer(start),
      end = as.integer(end)
    ) |>
    dplyr::select(gene_id, chromosome, start, end) |>
    dplyr::filter(!is.na(gene_id) & gene_id != "") |>
    dplyr::distinct(gene_id, .keep_all = TRUE)

  if (nrow(gene_features) == 0) {
    stop(paste0("No valid 'gene' features found in: ", filepath))
  }

  message(paste0("Found ", nrow(gene_features), " unique gene features from ", filepath, "."))
  return(gene_features)
}

#' @title Filter and Sanitize GFF Models
#' @description Sanitizes identifiers for OrthoFinder compatibility and filters
#' gene models against the pre-computed chromosome size thresholds.
#'
#' @param species_name character. Target species identifier.
#' @param config A `linguine_config` object.
#'
#' @return tibble. The filtered gene dataframe.
#' @noRd
filter_gff <- function(species_name, config) {

  gff_path <- file.path(config$paths$raw_data, paste0(species_name, ".gff"))
  if (!file.exists(gff_path)) {
    gff_path <- file.path(config$paths$raw_data, paste0(species_name, ".gff3"))
  }
  if (!file.exists(gff_path)) stop(paste("GFF file not found:", gff_path))

  genes_raw <- read_gff_genes(gff_path, species_name)

  genes_raw <- genes_raw |>
    dplyr::mutate(
      gene_id = gsub(":", "_", gene_id),
      chromosome = gsub(":", "_", chromosome)
    )

  chr_sizes_path <- file.path(config$paths$processed_data, paste0(species_name, "_chromosome_sizes.rds"))
  if (!file.exists(chr_sizes_path)) {
    stop(paste0("Error: Chromosome sizes missing for ", species_name, ". Execute extraction prior to filtering."))
  }
  chromosome_sizes_df <- readRDS(chr_sizes_path)

  chromosome_sizes_df <- chromosome_sizes_df |>
    dplyr::mutate(ref_chromosome = gsub(":", "_", ref_chromosome))

  message("\n--- Filtering out genes on scaffolds < ", scales::comma(config$min_chromosome_length_bp), " bp ---")

  initial_genes_count <- nrow(genes_raw)
  initial_chromosomes_count <- length(unique(genes_raw$chromosome))

  genes_filtered <- genes_raw |>
    dplyr::inner_join(dplyr::select(chromosome_sizes_df, chromosome = ref_chromosome), by = "chromosome")

  filtered_genes_count <- nrow(genes_filtered)
  filtered_chromosomes_count <- length(unique(genes_filtered$chromosome))

  message(paste0(species_name, " genes filtered: Removed ", initial_genes_count - filtered_genes_count, " genes."))
  message(paste0(species_name, " chromosomes filtered: Reduced from ", initial_chromosomes_count, " to ", filtered_chromosomes_count, " target sequences."))

  if (nrow(genes_filtered) == 0) {
    stop(paste0("Error: Zero genes retained for ", species_name, " post-filtering. Re-evaluate length thresholds."))
  }

  output_path <- file.path(config$paths$processed_data, paste0(species_name, "_genes_df.rds"))
  saveRDS(genes_filtered, file = output_path)
  message(paste0("Processed gene data saved to '", output_path, "'."))

  return(genes_filtered)
}

#' @title Vectorized Orthogroup Mapping
#' @description Maps a list of gene IDs to their corresponding Hierarchical Orthogroups (HOGs).
#'
#' @param gene_list character vector. A list of gene identifiers.
#' @param mapping_df data.frame. A dictionary linking Gene_ID to Target_HOG.
#'
#' @return character vector. A unique list of mapped HOGs.
#' @noRd
map_genes_to_hogs <- function(gene_list, mapping_df) {
  if (length(gene_list) == 0) return(character(0))
  dict <- setNames(mapping_df$Target_HOG, mapping_df$Gene_ID)
  mapped_hogs <- dict[gene_list]
  mapped_hogs <- unique(mapped_hogs[!is.na(mapped_hogs)])
  return(as.character(mapped_hogs))
}

#' @title Filter Chromosomal Observations
#' @description Extracts valid HMM state observations (those starting with S_ or ON_) from a vector.
#'
#' @param obs character vector. Raw observations from the HMM pipeline.
#'
#' @return character vector. Valid on-chromosome observations.
#' @noRd
get_on_chrom_observations <- function(obs) {
  obs[stringr::str_detect(obs, "^S_|^ON_")]
}

#' @title Determine Dominant Chromosome
#' @description Identifies the most frequent physical chromosome target within a rolling window.
#'
#' @param observations character vector. HMM state assignments in a local window.
#'
#' @return character. The name of the dominant target chromosome, or 'NON_SYN_BLOCK' if none.
#' @noRd
get_dominant_chr <- function(observations) {
  on_chrom_obs <- get_on_chrom_observations(observations)
  if (length(on_chrom_obs) == 0) {
    return("NON_SYN_BLOCK")
  }
  comp_chrs <- sub("S_|ON_", "", on_chrom_obs)
  dominant <- names(sort(table(comp_chrs), decreasing = TRUE))[1]
  return(dominant)
}

#' @title Calculate Synteny Block Metrics
#' @description Computes purity and abundance metrics for a localized sequence of orthology matches.
#'
#' @param raw_observations character vector. Hidden Markov Model output states.
#'
#' @return list. Contains `purity` (numeric), `count` (integer), and `dominant` (character).
#' @noRd
calculate_block_metrics <- function(raw_observations) {
  on_chrom_obs <- get_on_chrom_observations(raw_observations)
  if (length(on_chrom_obs) == 0) {
    return(list(purity = 0, count = 0L, dominant = NA_character_))
  }
  comp_chrs <- sub("ON_|S_", "", on_chrom_obs)
  chr_table <- table(comp_chrs)
  sorted_chrs <- names(sort(chr_table, decreasing = TRUE))
  max_freq <- if (length(chr_table) == 0) 0 else max(chr_table)
  total_on_genes <- sum(chr_table)
  metrics <- list(
    purity = max_freq / total_on_genes,
    count = as.integer(length(unique(comp_chrs))),
    dominant = sorted_chrs[1]
  )
  return(metrics)
}

#' @title Delineate Syntenic Blocks via Sliding Window
#' @description Scans a chromosomal dataframe to partition the sequence into contiguous
#' syntenic blocks based on local density and purity thresholds.
#'
#' @param chr_data tibble. Spatial gene coordinates and initial `inferred_state_hmm1`.
#' @param config A `linguine_config` object.
#'
#' @return tibble. Consolidated coordinate blocks with finalized syntenic states.
#' @export
classify_chromosome_blocks <- function(chr_data, config) {

  window_size <- config$thresholds$window_size
  purity_threshold <- config$thresholds$purity_threshold
  min_final_block_size <- config$thresholds$min_final_block_size

  # 1. Establish localized consensus using a rolling window
  chr_data_local_dominance <- chr_data |>
    dplyr::arrange(start) |>
    dplyr::mutate(
      local_dominant_chr = purrr::map_chr(
        dplyr::row_number(),
        ~get_dominant_chr(
          inferred_state_hmm1[max(1, .x - floor(window_size/2)):min(dplyr::n(), .x + floor(window_size/2))]
        )
      )
    )

  # 2. Extract discrete continuous regions via Run Length Encoding (RLE)
  initial_blocks_summary <- chr_data_local_dominance |>
    dplyr::mutate(
      rle_id = with(rle(local_dominant_chr), rep(seq_along(lengths), lengths))
    ) |>
    dplyr::group_by(rle_id) |>
    dplyr::summarise(
      chromosome = dplyr::first(chromosome),
      start_pos = min(start),
      end_pos = max(end),
      block_length = dplyr::n(),
      raw_observations = list(inferred_state_hmm1),
      .groups = "drop"
    )

  # 3. Apply state classifications based on established metrics
  initial_block_classified <- initial_blocks_summary |>
    dplyr::mutate(
      purity_metric = purrr::map_dbl(raw_observations, ~calculate_block_metrics(.x)$purity),
      unique_on_chroms_count = purrr::map_int(raw_observations, ~calculate_block_metrics(.x)$count),
      dominant_chr = purrr::map_chr(raw_observations, ~calculate_block_metrics(.x)$dominant, .default = NA_character_)
    ) |>
    dplyr::select(-raw_observations, -rle_id) |>
    dplyr::mutate(
      inferred_state_final = dplyr::case_when(
        block_length < min_final_block_size ~ "UNCLEAR_NOISE",
        purity_metric >= purity_threshold ~ paste0("SYN_TO_", dominant_chr),
        purity_metric < purity_threshold & unique_on_chroms_count == 2 ~ "TWO_CHR_MIX",
        purity_metric < purity_threshold & unique_on_chroms_count > 2 ~ "COMPLEX_REARRANGED",
        TRUE ~ "UNCLEAR_NOISE"
      )
    )

  # 4. Consolidate sequentially adjacent blocks sharing identical states
  main_blocks_to_consolidate <- initial_block_classified |>
    dplyr::filter(inferred_state_final != "UNCLEAR_NOISE") |>
    dplyr::arrange(start_pos)

  if (nrow(main_blocks_to_consolidate) == 0) {
    return(dplyr::tibble(
      chromosome = dplyr::first(chr_data$chromosome), start_pos = integer(), end_pos = integer(),
      block_length = integer(), inferred_state_final = character(), dominant_chr = character(),
      purity_metric = double(), unique_on_chroms_count = integer()
    ))
  }

  final_consolidated_report <- main_blocks_to_consolidate |>
    dplyr::mutate(
      group_key = cumsum(inferred_state_final != dplyr::lag(inferred_state_final, default = dplyr::first(inferred_state_final)))
    ) |>
    dplyr::group_by(group_key) |>
    dplyr::summarise(
      chromosome = dplyr::first(chromosome),
      start_pos = min(start_pos),
      end_pos = max(end_pos),
      block_length = sum(block_length),
      inferred_state_final = dplyr::first(inferred_state_final),
      .groups = "drop"
    ) |>
    dplyr::select(-group_key)

  return(final_consolidated_report)
}

#' @title Spatially Map Blocks to Genes
#' @description Attaches syntenic block classifications back to the individual gene models using fast interval trees.
#'
#' @param genes_df tibble. The individual gene coordinates.
#' @param blocks_df tibble. The consolidated structural blocks.
#'
#' @return tibble. The gene dataframe with inferred block states appended.
#' @noRd
map_blocks_to_genes <- function(genes_df, blocks_df) {
  if (is.null(blocks_df) || nrow(blocks_df) == 0) {
    result_df <- genes_df |> dplyr::mutate(inferred_state_blocks = NA_character_)
    if ("block_id" %in% names(blocks_df)) result_df <- result_df |> dplyr::mutate(block_id = NA_character_)
    return(result_df)
  }

  genes_dt <- data.table::as.data.table(genes_df)
  blocks_dt <- data.table::as.data.table(blocks_df)

  # foverlaps requires integer coordinates and keys
  blocks_dt[, `:=`(block_start_pos = as.integer(block_start_pos), block_end_pos = as.integer(block_end_pos))]
  genes_dt[, `:=`(start = as.integer(start), end = as.integer(end))]
  
  data.table::setkey(blocks_dt, block_start_pos, block_end_pos)
  
  # perform overlap
  overlaps <- data.table::foverlaps(genes_dt, blocks_dt, by.x = c("start", "end"), by.y = c("block_start_pos", "block_end_pos"), type = "any", mult = "first")
  
  genes_df <- genes_df |>
    dplyr::mutate(inferred_state_blocks = overlaps$inferred_state_blocks)

  if ("block_id" %in% names(blocks_dt)) {
    genes_df <- genes_df |>
      dplyr::mutate(block_id = overlaps$block_id)
  }

  return(genes_df)
}

#' @title Synthesize Broad Ancestral Regions
#' @description Converts fully reconstructed ancestral nodes back into pseudo-chromosomal blocks
#' for downstream cross-node structural comparisons.
#'
#' @param ancestral_df tibble. A fully reconstructed ancestral genome object.
#'
#' @return tibble. A standardized block-mapping dataframe compatible with the paralogy detector.
#' @noRd
create_synthetic_broad_regions <- function(ancestral_df) {
  ancestral_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      broad_ref_gene_ids = list(unique(c(unlist(A_genes), unlist(B_genes)))),
      broad_comp_gene_ids = list(character(0)),
      total_unique_og_count = length(unique(c(unlist(A_orthogroups), unlist(B_orthogroups)))),
      broad_ref_start = 1L,
      broad_ref_end = total_unique_og_count,
      broad_comp_start = 1L,
      broad_comp_end = total_unique_og_count,
      num_strict_blocks_merged = 1L
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      linkage_group_name = ancestral_lg_name,
      chromosome = ancestral_lg_name,
      broad_ref_start, broad_ref_end,
      target_comp_chromosome = ancestral_lg_name,
      broad_comp_start, broad_comp_end,
      num_strict_blocks_merged,
      broad_ref_gene_ids, broad_comp_gene_ids
    )
}

#' @title Enrich Blocks with Orthogroup Membership
#' @description Maps underlying physical genes within a block to their global orthogroups.
#'
#' @param blocks_df tibble. The spatial linkage blocks.
#' @param orthogroups_df tibble. The global orthology dictionary.
#' @param gene_id_column character. The target column containing the gene IDs.
#'
#' @return tibble. The blocks dataframe with a nested list column of `orthogroups_in_block`.
#' @noRd
map_orthogroups_to_blocks <- function(blocks_df, orthogroups_df, gene_id_column) {
  if (!gene_id_column %in% names(blocks_df)) {
    stop(paste0("Error: Column '", gene_id_column, "' not found in blocks_df."))
  }
  gene_id_list <- blocks_df[[gene_id_column]]
  og_list <- purrr::map(gene_id_list, function(g_ids) {
    if (is.null(g_ids) || length(g_ids) == 0) return(character(0))
    genes_vec <- unlist(g_ids)
    orthogroups_df |>
      dplyr::filter(gene_id %in% genes_vec) |>
      dplyr::pull(Orthogroup) |>
      unique()
  })
  blocks_df$orthogroups_in_block <- og_list
  return(blocks_df)
}

#' @title Consolidate Fragmented Linkage Groups
#' @description Merges physically fragmented collinear blocks that map to the same ancestral state.
#'
#' @param broad_regions tibble. Unconsolidated spatial regions.
#'
#' @return tibble. Unified spatial regions.
#' @noRd
merge_fragmented_lgs <- function(broad_regions) {
  if (nrow(broad_regions) == 0) return(broad_regions)
  message(paste0("Merging fragments... Input: ", nrow(broad_regions), " blocks."))
  merged_df <- broad_regions |>
    dplyr::group_by(linkage_group_name) |>
    dplyr::summarise(
      chromosome = dplyr::first(chromosome),
      broad_ref_start = min(broad_ref_start),
      broad_ref_end   = max(broad_ref_end),
      target_comp_chromosome = dplyr::first(target_comp_chromosome),
      broad_comp_start = min(broad_comp_start),
      broad_comp_end   = max(broad_comp_end),
      num_strict_blocks_merged = sum(num_strict_blocks_merged),
      broad_ref_gene_ids = list(unique(unlist(broad_ref_gene_ids))),
      broad_comp_gene_ids = if("broad_comp_gene_ids" %in% names(broad_regions)) list(unique(unlist(broad_comp_gene_ids))) else list(NULL),
      orthogroups_in_block = if("orthogroups_in_block" %in% names(broad_regions)) list(unique(unlist(orthogroups_in_block))) else list(NULL),
      .groups = "drop"
    )
  message(paste0("Merging complete. Output: ", nrow(merged_df), " blocks."))
  return(merged_df)
}

#' @title Execute Paralogy Enrichment Testing
#' @description Performs combinatorial Fisher's Exact tests across all blocks within a lineage
#' to identify statistically significant paralogy or whole genome duplication signatures.
#'
#' @param blocks_with_og tibble. Linkage blocks mapped with their nested orthogroups.
#' @param species_name character. The target species identifier.
#' @param all_orthogroups character vector. The global universe of valid orthogroups.
#' @param config A `linguine_config` object containing p-value and odds thresholds.
#'
#' @return tibble. A list of statistically enriched block pairs indicating paralogy.
#' @noRd
detect_paralogy_signals <- function(blocks_with_og, species_name, all_orthogroups, config) {
  num_blocks <- nrow(blocks_with_og)
  block_pairs <- expand.grid(i = 1:num_blocks, j = 1:num_blocks) |> dplyr::filter(i < j)
  orthogroup_list <- blocks_with_og$orthogroups_in_block
  total_orthogroups <- length(unique(all_orthogroups))

  message(paste0("Scanning for paralogy in ", species_name, " (", num_blocks, " blocks)..."))

  results_df <- purrr::map_dfr(1:nrow(block_pairs), function(k) {
    i <- block_pairs$i[k]
    j <- block_pairs$j[k]
    og1 <- unlist(orthogroup_list[i])
    og2 <- unlist(orthogroup_list[j])

    if (length(og1) == 0 || length(og2) == 0) return(dplyr::tibble(p_value=1, odds_ratio=1))

    shared <- length(intersect(og1, og2))
    uniq1 <- length(og1) - shared
    uniq2 <- length(og2) - shared
    rest <- total_orthogroups - (uniq1 + uniq2 + shared)

    mat <- matrix(c(shared, uniq1, uniq2, rest), nrow=2)
    if (any(mat < 0)) return(dplyr::tibble(p_value=1, odds_ratio=1))

    ft <- fisher.test(mat, alternative = "greater")
    dplyr::tibble(p_value = ft$p.value, odds_ratio = ft$estimate)
  })
  if (nrow(results_df) == 0) return(dplyr::tibble())

  final_results <- dplyr::bind_cols(
    block1_id = blocks_with_og$linkage_group_name[block_pairs$i],
    block2_id = blocks_with_og$linkage_group_name[block_pairs$j],
    results_df
  ) |>
    dplyr::mutate(p_value_adjusted = p.adjust(p_value, method = "BH")) |>
    dplyr::filter(p_value_adjusted < config$thresholds$paralogy_p & odds_ratio > config$thresholds$paralogy_odds) |>
    dplyr::arrange(p_value_adjusted)

  return(final_results)
}

#' @title Compute Jaccard Similarity Matrix
#' @description Calculates fractional orthology overlap between the structural blocks
#' of two divergent lineages to establish homology bridges.
#'
#' @param broad_regions_ref tibble. Reference structural blocks.
#' @param broad_regions_comp tibble. Comparison structural blocks.
#'
#' @return matrix. A bidirectional intersection overlap matrix.
#' @noRd
calculate_overlap_matrix <- function(broad_regions_ref, broad_regions_comp) {
  if (nrow(broad_regions_ref) == 0 || nrow(broad_regions_comp) == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }
  ref_ogs <- broad_regions_ref$orthogroups_in_block
  comp_ogs <- broad_regions_comp$orthogroups_in_block

  ref_ogs <- lapply(ref_ogs, function(x) if (is.null(x)) character(0) else x)
  comp_ogs <- lapply(comp_ogs, function(x) if (is.null(x)) character(0) else x)

  overlap_list <- purrr::map(ref_ogs, function(r_og) {
    purrr::map_dbl(comp_ogs, function(c_og) {
      if (length(r_og) == 0 || length(c_og) == 0) return(0)
      length(intersect(r_og, c_og)) / length(union(r_og, c_og))
    })
  })

  matrix_out <- do.call(rbind, overlap_list)
  colnames(matrix_out) <- broad_regions_comp$linkage_group_name
  rownames(matrix_out) <- broad_regions_ref$linkage_group_name

  return(matrix_out)
}

#' @title Extract Syntenic Partners
#' @description Returns the most highly correlated structural targets across a homology matrix.
#'
#' @param lg_name character. The target reference linkage group.
#' @param overlap_matrix matrix. The pre-computed Jaccard overlap matrix.
#'
#' @return character vector. Names of matching homologous comparison blocks.
#' @noRd
get_syntenic_partners <- function(lg_name, overlap_matrix) {
  if (lg_name %in% rownames(overlap_matrix)) {
    partners <- sort(overlap_matrix[lg_name, ], decreasing = TRUE)
    return(names(partners[partners > 0]))
  }
  return(character(0))
}

#' @title Classify Duplication Timing
#' @description Evaluates intra-lineage paralogy against inter-lineage homology
#' to determine whether duplications occurred prior to speciation (ancestral) or post-speciation.
#'
#' @param paralogy_A tibble. Significant paralog blocks in Lineage A.
#' @param paralogy_B tibble. Significant paralog blocks in Lineage B.
#' @param broad_regions_A tibble. Base structural blocks in Lineage A.
#' @param broad_regions_B tibble. Base structural blocks in Lineage B.
#'
#' @return list. Categorized duplications (`ancestral_dups`, `lineage_A`, `lineage_B`).
#' @noRd
infer_duplication_events <- function(paralogy_A, paralogy_B, broad_regions_A, broad_regions_B) {
  matrix_A_B <- calculate_overlap_matrix(broad_regions_A, broad_regions_B)
  matrix_B_A <- calculate_overlap_matrix(broad_regions_B, broad_regions_A)

  ancestral_pairs <- dplyr::tibble(A_lg1=character(), A_lg2=character(), B_lg1=character(), B_lg2=character())
  lineage_A <- dplyr::tibble(A_lg1=character(), A_lg2=character(), B_lg=character())
  lineage_B <- dplyr::tibble(B_lg1=character(), B_lg2=character(), A_lg=character())

  if (nrow(paralogy_A) > 0 && nrow(paralogy_B) > 0) {
    new_anc <- purrr::map_dfr(1:nrow(paralogy_A), function(i) {
      A1 <- paralogy_A$block1_id[i]; A2 <- paralogy_A$block2_id[i]
      B_partners_1 <- get_syntenic_partners(A1, matrix_A_B)
      B_partners_2 <- get_syntenic_partners(A2, matrix_A_B)

      if (length(B_partners_1) > 0 && length(B_partners_2) > 0) {
        match <- paralogy_B |>
          dplyr::filter((block1_id %in% B_partners_1 & block2_id %in% B_partners_2) |
                          (block1_id %in% B_partners_2 & block2_id %in% B_partners_1))

        if (nrow(match) > 0) {
          return(dplyr::tibble(A_lg1=A1, A_lg2=A2, B_lg1=match$block1_id[1], B_lg2=match$block2_id[1]))
        }
      }
      return(NULL)
    })
    ancestral_pairs <- dplyr::bind_rows(ancestral_pairs, new_anc) |> dplyr::distinct()
  }

  used_A <- unique(c(ancestral_pairs$A_lg1, ancestral_pairs$A_lg2))
  used_B <- unique(c(ancestral_pairs$B_lg1, ancestral_pairs$B_lg2))

  if (nrow(paralogy_A) > 0) {
    new_lin_A <- purrr::map_dfr(1:nrow(paralogy_A), function(i) {
      A1 <- paralogy_A$block1_id[i]; A2 <- paralogy_A$block2_id[i]
      if (A1 %in% used_A || A2 %in% used_A) return(NULL)
      B_partners_1 <- get_syntenic_partners(A1, matrix_A_B)
      B_partners_2 <- get_syntenic_partners(A2, matrix_A_B)
      common <- intersect(B_partners_1, B_partners_2)
      if (length(common) > 0) return(dplyr::tibble(A_lg1=A1, A_lg2=A2, B_lg=common[1]))
      return(NULL)
    })
    lineage_A <- dplyr::bind_rows(lineage_A, new_lin_A) |> dplyr::distinct()
  }

  if (nrow(paralogy_B) > 0) {
    new_lin_B <- purrr::map_dfr(1:nrow(paralogy_B), function(i) {
      B1 <- paralogy_B$block1_id[i]; B2 <- paralogy_B$block2_id[i]
      if (B1 %in% used_B || B2 %in% used_B) return(NULL)
      A_partners_1 <- get_syntenic_partners(B1, matrix_B_A)
      A_partners_2 <- get_syntenic_partners(B2, matrix_B_A)
      common <- intersect(A_partners_1, A_partners_2)
      if (length(common) > 0) return(dplyr::tibble(B_lg1=B1, B_lg2=B2, A_lg=common[1]))
      return(NULL)
    })
    lineage_B <- dplyr::bind_rows(lineage_B, new_lin_B) |> dplyr::distinct()
  }

  return(list(ancestral_dups=ancestral_pairs, lineage_A=lineage_A, lineage_B=lineage_B))
}

#' @title Graph-Based Linkage Group Consolidation
#' @description Merges distinct spatial blocks into unified sub-graph components
#' utilizing inferred paralogy networks via the igraph engine.
#'
#' @param broad_regions tibble. Raw spatial blocks.
#' @param pairs_to_collapse tibble. Significant paralogous pair linkages.
#'
#' @return tibble. Unified linkage groups post-consolidation.
#' @noRd
collapse_lg_pairs <- function(broad_regions, pairs_to_collapse) {
  if (nrow(broad_regions) == 0) {
    return(dplyr::tibble(linkage_group_name = character(), orthogroups_in_block = list(), genes_in_block = list(), original_regions = list()))
  }

  broad_regions <- dplyr::ungroup(broad_regions)

  if ("broad_ref_gene_ids" %in% names(broad_regions)) {
    names(broad_regions)[names(broad_regions) == "broad_ref_gene_ids"] <- "genes_in_block"
  }

  if (nrow(pairs_to_collapse) == 0) {
    return(broad_regions |>
             dplyr::mutate(original_regions = purrr::pmap(list(linkage_group_name, chromosome, broad_ref_start, broad_ref_end),
                                                          function(lg, chr, s, e) dplyr::tibble(linkage_group_name=lg, chromosome=chr, broad_ref_start=s, broad_ref_end=e))) |>
             dplyr::select(linkage_group_name, orthogroups_in_block, genes_in_block, original_regions))
  }

  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package 'igraph' is required.")

  g <- igraph::graph_from_data_frame(pairs_to_collapse, directed = FALSE)
  comps <- igraph::components(g)

  cluster_map <- data.frame(
    old_lg = names(comps$membership),
    cluster_id = as.numeric(comps$membership),
    stringsAsFactors = FALSE
  )

  new_names <- aggregate(old_lg ~ cluster_id, data = cluster_map, FUN = function(x) paste(sort(unique(x)), collapse = "_"))
  names(new_names)[2] <- "new_lg"

  cluster_map <- merge(cluster_map, new_names, by = "cluster_id")

  lgs_in_clusters <- cluster_map$old_lg
  to_collapse <- broad_regions |> dplyr::filter(linkage_group_name %in% lgs_in_clusters)
  untouched <- broad_regions |> dplyr::filter(!linkage_group_name %in% lgs_in_clusters)

  if(nrow(to_collapse) > 0) {
    to_collapse_df <- as.data.frame(to_collapse)
    joined_df <- merge(to_collapse_df, cluster_map, by.x = "linkage_group_name", by.y = "old_lg")

    if(!"new_lg" %in% names(joined_df)) stop("Error: Core mapping lost during structural join.")

    split_data <- split(joined_df, joined_df$new_lg)

    merged_list <- lapply(names(split_data), function(grp_name) {
      sub_df <- split_data[[grp_name]]
      all_ogs <- unique(unlist(sub_df$orthogroups_in_block))
      all_genes <- unique(unlist(sub_df$genes_in_block))

      orig_regs <- dplyr::tibble(
        linkage_group_name = sub_df$linkage_group_name,
        chromosome = sub_df$chromosome,
        broad_ref_start = sub_df$broad_ref_start,
        broad_ref_end = sub_df$broad_ref_end
      )

      dplyr::tibble(
        linkage_group_name = grp_name,
        orthogroups_in_block = list(all_ogs),
        genes_in_block = list(all_genes),
        original_regions = list(orig_regs)
      )
    })

    merged <- dplyr::bind_rows(merged_list)
  } else {
    merged <- dplyr::tibble(linkage_group_name = character(), orthogroups_in_block = list(), genes_in_block = list(), original_regions = list())
  }

  untouched_processed <- untouched |>
    dplyr::mutate(original_regions = purrr::pmap(list(linkage_group_name, chromosome, broad_ref_start, broad_ref_end),
                                                 function(lg, chr, s, e) dplyr::tibble(linkage_group_name=lg, chromosome=chr, broad_ref_start=s, broad_ref_end=e))) |>
    dplyr::select(linkage_group_name, orthogroups_in_block, genes_in_block, original_regions)

  result <- dplyr::bind_rows(merged, untouched_processed)
  return(result)
}

#' @title Clean Linkage Group Partitions
#' @description Merges nested data lists following graph-based consolidation to ensure clean downstream mapping.
#'
#' @param collapsed_lg_df tibble. The semi-processed dataframe from graph mapping.
#'
#' @return tibble. Dataframe with fully unified list arrays.
#' @noRd
merge_lg_fragments <- function(collapsed_lg_df) {
  if (nrow(collapsed_lg_df) == 0) {
    return(collapsed_lg_df)
  }
  merged_df <- collapsed_lg_df |>
    dplyr::group_by(linkage_group_name) |>
    dplyr::summarise(
      orthogroups_in_block = list(unique(unlist(orthogroups_in_block))),
      genes_in_block = list(unique(unlist(genes_in_block))),
      original_regions = list(dplyr::bind_rows(original_regions)),
      .groups = 'drop'
    )
  return(merged_df)
}

#' @title Render Oxford Grid (Dotplot)
#' @description Generates the underlying ggplot and caching for 2D ancestral synteny mapping.
#'
#' @param species_name character. Reference species identifier.
#' @param other_species_name character. Target mapping entity.
#' @param chromosome_sizes_df data.frame. Physical scaffold sizes.
#' @param gene_id_col character. Reference gene string column name.
#' @param chromosome_col character. Physical chromosome column name.
#' @param start_col character. Base pair start coordinate column.
#' @param end_col character. Base pair end coordinate column.
#' @param spA_vs_spB_data tibble. Processed intersection data mapping species coordinates.
#' @param og_master_list tibble. Consolidated orthology boundaries.
#' @param LG_boundaries tibble. Synthesized graph block borders.
#' @param LG_labels tibble. Calculated plotting midpoints.
#' @param min_chromosome_length_bp numeric. Configuration filter threshold.
#' @param lg_colors character vector. Extracted color palette dictionary.
#' @param config A `linguine_config` list containing file paths and system parameters.
#'
#' @return ggplot object. The resulting scatter plot topology.
#' @noRd
generate_synteny_plot <- function(
    species_name, other_species_name, chromosome_sizes_df, gene_id_col,
    chromosome_col, start_col, end_col, spA_vs_spB_data, og_master_list,
    LG_boundaries, LG_labels, min_chromosome_length_bp, lg_colors, config
) {
  chromosome_lengths <- chromosome_sizes_df |>
    dplyr::filter(chromosome_length_bp > min_chromosome_length_bp) |>
    dplyr::select(chromosome = dplyr::all_of("ref_chromosome"), length = chromosome_length_bp) |>
    dplyr::arrange(chromosome) |>
    dplyr::mutate(cumulative_length = dplyr::lag(cumsum(length), default = 0)) |>
    dplyr::select(chromosome, length, cumulative_length)

  chromosome_midpoints <- chromosome_lengths |>
    dplyr::distinct(chromosome, length, cumulative_length) |>
    dplyr::mutate(midpoint = cumulative_length + (length / 2))

  plotting_data <- spA_vs_spB_data |>
    dplyr::left_join(chromosome_lengths, by = setNames("chromosome", chromosome_col)) |>
    dplyr::mutate(cumulative_start = .data[[start_col]] + cumulative_length) |>
    dplyr::left_join(dplyr::select(og_master_list, orthogroup, x_position), by = "orthogroup") |>
    dplyr::filter(!is.na(x_position)) |>
    dplyr::distinct(dplyr::across(dplyr::all_of(gene_id_col)), .keep_all = TRUE) |>
    dplyr::filter(assignment_tag != "multi_LG")

  LGs_df <- plotting_data |>
    dplyr::select(gene_id = dplyr::all_of(gene_id_col),
                  chromosome = dplyr::all_of(chromosome_col),
                  start = dplyr::all_of(start_col),
                  end = dplyr::all_of(end_col),
                  orthogroup, LG = assignment_tag)

  # CORRECTED: Save directly to the dynamically generated results folder!
  write.table(
    LGs_df,
    file.path(config$paths$results, paste0("LGs_", species_name, "_vs_", other_species_name, ".txt")),
    quote = FALSE, row.names = FALSE, sep = "\t"
  )

  dotplot <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = x_position, y = cumulative_start, color = assignment_tag)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.5) +
    ggplot2::labs(
      title = paste(species_name, " Genes vs. Ancestral Syntenic Blocks", sep = ""),
      x = "Linkage groups",
      y = "Chromosomes",
      color = "LG Assignment"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = lg_colors) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(vjust = -2),
      axis.title.y = ggplot2::element_text(vjust = 20),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(1, 1, 3, 5), "lines"),
      strip.background = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    ) +
    ggplot2::geom_vline(
      data = LG_boundaries[LG_boundaries$assignment_tag != "multi_LG", ],
      ggplot2::aes(xintercept = end_position + 0.5),
      color = "grey", linetype = "dotted", linewidth = 0.3
    ) +
    ggplot2::geom_hline(
      data = chromosome_lengths,
      ggplot2::aes(yintercept = cumulative_length + length),
      color = "grey", linetype = "dotted", linewidth = 0.3
    ) +
    ggplot2::annotate(
      "text",
      x = LG_labels$midpoint[LG_labels$assignment_tag != "multi_LG"],
      y = -Inf,
      label = LG_labels$assignment_tag[LG_labels$assignment_tag != "multi_LG"],
      vjust = -0.5, angle = 90, size = 2.5
    ) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = chromosome_midpoints$midpoint,
      label = chromosome_midpoints$chromosome,
      hjust = 1.1, size = 3
    ) +
    ggplot2::coord_cartesian(clip = "off")

  return(dotplot)
}

#' @title Extract Numerical Tree Index
#' @description Maps a string label safely to its underlying Ape Tree internal numeric index.
#'
#' @param tree phylo. The phylogenetic species tree.
#' @param label character. The node or species label to search.
#'
#' @return integer. The parsed index, or NULL if unavailable.
#' @noRd
get_node_index <- function(tree, label) {
  if (label %in% tree$tip.label) {
    return(which(tree$tip.label == label))
  }
  if (!is.null(tree$node.label) && label %in% tree$node.label) {
    return(length(tree$tip.label) + which(tree$node.label == label))
  }
  if (grepl("^INode_", label)) {
    idx <- as.integer(sub("INode_", "", label))
    if (!is.na(idx)) return(idx)
  }
  return(NULL)
}

#' @title Dynamic Outgroup Identification
#' @description Travels backwards up the species tree from the Most Recent Common Ancestor
#' to locate a suitable sister lineage for outgroup polarization mapping.
#'
#' @param tree phylo. The species tree.
#' @param spA character. Reference node/lineage.
#' @param spB character. Comparison node/lineage.
#'
#' @return character. The sister outgroup identifier, or NULL if at the root.
#' @noRd
get_outgroup_species <- function(tree, spA, spB) {
  idxA <- get_node_index(tree, spA)
  idxB <- get_node_index(tree, spB)

  if (is.null(idxA) || is.null(idxB)) return(NULL)

  mrca <- ape::getMRCA(tree, c(idxA, idxB))
  if (is.null(mrca)) return(NULL)

  if (mrca == (length(tree$tip.label) + 1)) return(NULL)

  parent_edge <- which(tree$edge[,2] == mrca)
  if (length(parent_edge) == 0) return(NULL)
  parent <- tree$edge[parent_edge, 1]

  siblings <- tree$edge[which(tree$edge[,1] == parent), 2]
  sister <- setdiff(siblings, mrca)

  if (length(sister) == 0) return(NULL)

  if (sister[1] <= length(tree$tip.label)) {
    return(tree$tip.label[sister[1]])
  } else {
    tips_in_sister <- ape::extract.clade(tree, sister[1])$tip.label
    if (length(tips_in_sister) > 0) return(tips_in_sister[1])
  }
  return(NULL)
}

#' @title Genome Retrieval Agent
#' @description Safely retrieves processed extant genome data if the request involves a tip.
#'
#' @param species_name character. Target string label.
#' @param is_node logical. Flag indicating if target is internal graph node.
#' @param processed_dir character. System path.
#' @param results_dir character. System path.
#'
#' @return list. Contains the type and corresponding dataframe, or NULL.
#' @noRd
load_genome_content <- function(species_name, is_node, processed_dir, results_dir) {
  if (!is_node) {
    fpath <- file.path(processed_dir, paste0(species_name, "_genes_df.rds"))
    if (file.exists(fpath)) {
      return(list(type = "tip", data = readRDS(fpath)))
    }
  } else {
    return(NULL)
  }
  return(NULL)
}

#' @title Isolate Valid Syntenic Spans
#' @description Drops highly fragmented linkage blocks falling below mathematical relevance thresholds.
#'
#' @param df tibble. Raw spatial regions.
#' @param min_ogs numeric. Minimum hard count of valid orthogroups inside block.
#' @param rel_threshold numeric. Proportion threshold for dynamic retention.
#'
#' @return tibble. The filtered broad regions frame.
#' @noRd
clean_broad_regions <- function(df, min_ogs = 20, rel_threshold = 0.01) {
  total_ogs_in_species <- sum(sapply(df$orthogroups_in_block, length))

  df |>
    dplyr::rowwise() |>
    dplyr::mutate(og_count = length(orthogroups_in_block)) |>
    dplyr::ungroup() |>
    dplyr::filter(og_count >= min_ogs | og_count >= (total_ogs_in_species * rel_threshold)) |>
    dplyr::select(-og_count)
}

#' @title Format Visual Output Vectors
#' @description Flattens nested array strings into semicolons for clean graphical plotting tags.
#'
#' @param lg_list_col list. Array of linkage group strings.
#'
#' @return character vector. Unified, semi-colon separated string outputs.
#' @noRd
flatten_lg_list <- function(lg_list_col) {
  purrr::map_chr(lg_list_col, ~ paste(unique(.x), collapse = ";"))
}

#' @title Parse Global Orthology Data
#' @description Reads and transforms standard OrthoFinder outputs (either Hierarchical
#' Orthogroups or standard flat Orthogroups) into a standardized, long-format
#' relational index.
#'
#' @param config A \code{linguine_config} object.
#'
#' @return tibble. A standardized dataframe mapping Gene_ID to Orthogroup/HOG.
#' @noRd
parse_global_orthology <- function(config) {

  if (config$orthology_type == "HOGs") {
    message("Processing Hierarchical Orthogroups (HOGs)...")

    # Assuming config$orthology_filename points to the directory containing N0.tsv, N1.tsv, etc.
    hog_dir <- file.path(config$paths$raw_data, config$orthology_filename)
    ortholog_data_files <- list.files(hog_dir, pattern = "\\.tsv$", full.names = TRUE)

    if (length(ortholog_data_files) == 0) {
      stop("Error: No TSV files found in HOG directory: ", hog_dir)
    }

    # Safely extract node names (e.g., "N0", "N1") from filenames
    node_names <- sapply(basename(ortholog_data_files), function(x) sub("\\.tsv$", "", x))
    names(ortholog_data_files) <- node_names
    ortholog_data_files <- ortholog_data_files[gtools::mixedorder(names(ortholog_data_files))]

    message(paste0("Processing ", length(ortholog_data_files), " HOG files..."))

    list_of_hog_dfs <- lapply(seq_along(ortholog_data_files), function(i) {
      file_path <- ortholog_data_files[i]
      node_id <- names(ortholog_data_files)[i]

      node_data_raw <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
      hog_col_name <- "HOG"
      new_col_name <- paste0("HOG_", node_id)

      species_cols_start <- which(names(node_data_raw) == "Gene Tree Parent Clade") + 1

      node_hogs_long <- node_data_raw |>
        dplyr::select(Orthogroup = OG, rlang::sym(hog_col_name), dplyr::all_of(names(node_data_raw)[species_cols_start:ncol(node_data_raw)])) |>
        tidyr::pivot_longer(
          cols = -c(Orthogroup, rlang::sym(hog_col_name)),
          names_to = "Species",
          values_to = "Gene_IDs_List"
        ) |>
        dplyr::filter(Gene_IDs_List != "" & Gene_IDs_List != "-") |>
        tidyr::separate_rows(Gene_IDs_List, sep = ",\\s*") |>
        dplyr::mutate(Gene_ID = trimws(Gene_IDs_List)) |>
        dplyr::select(Gene_ID, Species, Orthogroup, !!rlang::sym(new_col_name) := rlang::sym(hog_col_name))

      return(node_hogs_long)
    })

    message("Executing full spatial join across all hierarchical nodes...")
    consolidated_hogs_matrix <- purrr::reduce(list_of_hog_dfs, function(x, y) {
      dplyr::full_join(x, y, by = c("Gene_ID", "Species", "Orthogroup"))
    })

    message(paste0("Hierarchical integration complete. Total unique genes mapped: ", nrow(consolidated_hogs_matrix)))
    return(consolidated_hogs_matrix)

  } else {
    message("Processing Standard flat Orthogroups (OGs)...")

    og_file <- file.path(config$paths$raw_data, config$orthology_filename)
    if (!file.exists(og_file)) {
      stop("Error: Orthogroups file not found at: ", og_file)
    }

    orthogroups_raw <- read.delim(og_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

    orthogroups_final <- orthogroups_raw |>
      tidyr::pivot_longer(cols = -Orthogroup, names_to = "Species", values_to = "Gene_IDs_List") |>
      dplyr::filter(Gene_IDs_List != "") |>
      tidyr::separate_rows(Gene_IDs_List, sep = ", ") |>
      dplyr::mutate(Gene_ID = trimws(Gene_IDs_List)) |>
      dplyr::select(Gene_ID, Species, Orthogroup)

    message("Standard OG parsing complete.")
    return(orthogroups_final)
  }
}
