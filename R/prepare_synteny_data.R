#' @title Prepare Genomic and Orthology Data for HMM
#'
#' @description Merges physical gene coordinates with global orthology data to generate
#' the raw sequence of emission states required by the Hidden Markov Model. It gracefully
#' handles both extant-to-extant (tip_vs_tip) and extant-to-ancestor (tip_vs_node) mappings.
#'
#' @param ref_species Character. Reference lineage identifier.
#' @param comp_species Character. Comparison lineage identifier.
#' @param comparison_type Character. Architecture of the comparison ("tip_vs_tip", "tip_vs_node", etc.).
#' @param parent_node Character. Target internal node representing the most recent common ancestor.
#' @param comp_node_daughters Character vector. The two immediate descending lineages of the comparison node.
#' @param config A \code{linguine_config} object created by \code{create_linguine_config()}.
#'
#' @return Character path to the saved HMM observation RDS file, or NULL if skipped.
#' @export
prepare_synteny_data <- function(ref_species, comp_species, comparison_type, parent_node, comp_node_daughters, config) {

  message("\n--- Integrating Genomic and Orthology Data (", ref_species, " vs ", comp_species, ") ---")

  # Short-circuit immediately for INode to prevent loading extant gene files
  if (comparison_type == "INode_vs_INode") {
    message("Skipping preliminary HMM emission mapping for strictly structural Node-vs-Node analysis.")
    return(NULL)
  }

  # ==============================================================================
  # 1. Load Extant Gene Coordinates
  # ==============================================================================
  ref_genes_df_path <- file.path(config$paths$processed_data, paste0(ref_species, "_genes_df.rds"))
  if (!file.exists(ref_genes_df_path)) stop("Error: Gene mapping missing for ", ref_species)
  ref_genes_raw <- readRDS(ref_genes_df_path)

  if (comparison_type == "tip_vs_tip") {
    comp_genes_df_path <- file.path(config$paths$processed_data, paste0(comp_species, "_genes_df.rds"))
    if (!file.exists(comp_genes_df_path)) stop("Error: Gene mapping missing for ", comp_species)
    comp_genes_raw <- readRDS(comp_genes_df_path)
  }

  # ==============================================================================
  # 2. Ingest and Standardize Global Orthology
  # ==============================================================================
  ortholog_data_path <- file.path(config$paths$processed_data, paste0("ortholog_data_", config$orthology_type, ".rds"))

  if (!file.exists(ortholog_data_path)) {
    if (config$orthology_type == "HOGs") {
      ortholog_data_files <- list.files(file.path(config$paths$raw_data, config$orthology_filename), full.names = TRUE)
      names(ortholog_data_files) <- sapply(strsplit(ortholog_data_files, "/"), function(x) sub(".tsv", "", x[grep(".tsv", x)]))
      ortholog_data_files <- ortholog_data_files[gtools::mixedorder(names(ortholog_data_files))]

      message("Processing ", length(ortholog_data_files), " HOG files...")
      list_of_hog_dfs <- lapply(seq_along(ortholog_data_files), function(i) {
        file_path <- ortholog_data_files[i]
        node_id <- names(ortholog_data_files)[i]
        node_data_raw <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

        hog_col_name <- "HOG"
        new_col_name <- paste0("HOG_", node_id)
        species_cols_start <- which(names(node_data_raw) == "Gene Tree Parent Clade") + 1

        node_data_raw |>
          dplyr::select(Orthogroup = OG, !!dplyr::sym(hog_col_name), dplyr::all_of(names(node_data_raw)[species_cols_start:ncol(node_data_raw)])) |>
          tidyr::pivot_longer(cols = -c(Orthogroup, !!dplyr::sym(hog_col_name)), names_to = "Species", values_to = "Gene_IDs_List") |>
          dplyr::filter(Gene_IDs_List != "" & Gene_IDs_List != "-") |>
          tidyr::separate_rows(Gene_IDs_List, sep = ",\\s*") |>
          dplyr::mutate(Gene_ID = trimws(Gene_IDs_List)) |>
          dplyr::select(Gene_ID, Species, Orthogroup, !!dplyr::sym(new_col_name) := !!dplyr::sym(hog_col_name))
      })

      message("Executing full spatial join across all hierarchical nodes...")
      orthogroups_final <- Reduce(function(x, y) dplyr::full_join(x, y, by = c("Gene_ID", "Species", "Orthogroup")), list_of_hog_dfs)
      saveRDS(orthogroups_final, file = ortholog_data_path)
      message("Hierarchical integration complete. Total unique genes mapped: ", nrow(orthogroups_final))

    } else {
      orthogroups_raw <- read.delim(file.path(config$paths$raw_data, config$orthology_filename), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
      orthogroups_final <- orthogroups_raw |>
        tidyr::pivot_longer(cols = -Orthogroup, names_to = "Species", values_to = "Gene_IDs_List") |>
        dplyr::filter(Gene_IDs_List != "") |>
        tidyr::separate_rows(Gene_IDs_List, sep = ", ") |>
        dplyr::mutate(Gene_ID = trimws(Gene_IDs_List)) |>
        dplyr::select(Gene_ID, Species, Orthogroup)
      saveRDS(orthogroups_final, file = ortholog_data_path)
      message("Standard OG parsing complete.")
    }
  } else {
    orthogroups_final <- readRDS(ortholog_data_path)
  }

  # ==============================================================================
  # 3. Execute Pairwise Filtration
  # ==============================================================================
  ortholog_data_filtered_path <- file.path(config$paths$processed_data, paste0(ref_species, "_vs_", comp_species, "_ortholog_data_filtered.rds"))

  if (!file.exists(ortholog_data_filtered_path)) {
    if (config$orthology_type == "OGs") {
      id_col_name <- "Orthogroup"
    } else {
      id_col_name <- paste0("HOG_", parent_node)
      if (!id_col_name %in% colnames(orthogroups_final)) {
        if (parent_node == "N0") id_col_name <- "Orthogroup" else stop("Error: Column ", id_col_name, " is missing.")
      }
    }

    filtered_df <- orthogroups_final |>
      dplyr::select(Gene_ID, Species, OG_ID = !!dplyr::sym(id_col_name)) |>
      dplyr::filter(!is.na(OG_ID)) |>
      dplyr::group_by(OG_ID) |>
      dplyr::filter(dplyr::n_distinct(Species) > 1) |>
      dplyr::ungroup()

    ref_df <- filtered_df |> dplyr::filter(Species == ref_species) |> dplyr::select(Gene_ID_Ref = Gene_ID, OG_ID)

    if (comparison_type == "tip_vs_tip") {
      comp_df <- filtered_df |> dplyr::filter(Species == comp_species) |> dplyr::select(Gene_ID_Comp = Gene_ID, OG_ID)
      ortholog_temp <- ref_df |> dplyr::left_join(comp_df, by = "OG_ID", relationship = "many-to-many") |> dplyr::select(ref_gene_id = Gene_ID_Ref, comp_gene_id = Gene_ID_Comp, OG_ID)
      ortholog_data_filtered <- ortholog_temp |>
        dplyr::inner_join(ref_genes_raw |> dplyr::select(ref_gene_id = gene_id, ref_chromosome = chromosome, ref_start = start, ref_end = end), by = "ref_gene_id") |>
        dplyr::inner_join(comp_genes_raw |> dplyr::select(comp_gene_id = gene_id, comp_chromosome = chromosome, comp_start = start, comp_end = end), by = "comp_gene_id")
    } else {
      ortholog_temp <- ref_df |> dplyr::select(ref_gene_id = Gene_ID_Ref, OG_ID)
      ortholog_data_filtered <- ortholog_temp |>
        dplyr::inner_join(ref_genes_raw |> dplyr::select(ref_gene_id = gene_id, ref_chromosome = chromosome, ref_start = start, ref_end = end), by = "ref_gene_id")
    }

    if (nrow(ortholog_data_filtered) == 0) stop("Error: No viable ortholog pairs remained post-spatial mapping.")
    saveRDS(ortholog_data_filtered, file = ortholog_data_filtered_path)
  } else {
    ortholog_data_filtered <- readRDS(ortholog_data_filtered_path)
  }

  # ==============================================================================
  # 4. HMM Emission Sequence Generation
  # ==============================================================================
  hmm_observation_path <- file.path(config$paths$intermediate_data, paste0("hmm_observations_detailed_", ref_species, "_ref_", comp_species, "_comp.rds"))

  if (file.exists(hmm_observation_path)) {
    return(hmm_observation_path)
  }

  message("Computing Structural HMM Emission Sequences...")

  # Filter out unplaced/small scaffolds before HMM
  min_genes <- if(is.null(config$thresholds$min_genes_per_scaffold)) 10 else config$thresholds$min_genes_per_scaffold
  
  valid_chroms <- ref_genes_raw |>
    dplyr::count(chromosome) |>
    dplyr::filter(n >= min_genes) |>
    dplyr::pull(chromosome)
    
  unplaced_genes <- ref_genes_raw |>
    dplyr::filter(!chromosome %in% valid_chroms)
    
  ref_genes_filtered <- ref_genes_raw |>
    dplyr::filter(chromosome %in% valid_chroms)

  unplaced_path <- file.path(config$paths$intermediate_data, paste0("unplaced_scaffolds_", ref_species, "_ref_", comp_species, "_comp.rds"))
  saveRDS(unplaced_genes, unplaced_path)

  message(sprintf("Filtering small scaffolds (< %d genes): Kept %d valid chromosomes. Saved %d unplaced genes for post-hoc assignment.", min_genes, length(valid_chroms), nrow(unplaced_genes)))

  if (comparison_type == "tip_vs_tip") {
    gene_comp_chroms <- ortholog_data_filtered |>
      dplyr::group_by(ref_gene_id) |>
      dplyr::count(comp_chromosome, name = "count_chr") |>
      dplyr::filter(count_chr == max(count_chr)) |>
      dplyr::summarise(primary_comp_chr = sample(comp_chromosome, 1), .groups = 'drop')

    gene_ortholog_summary <- ortholog_data_filtered |>
      dplyr::group_by(ref_gene_id) |>
      dplyr::summarise(num_ortholog_hits = dplyr::n(), all_comp_chromosomes = list(unique(comp_chromosome)), Orthogroup = unique(OG_ID), .groups = 'drop')

    hmm_observation_df <- ref_genes_filtered |>
      dplyr::left_join(gene_comp_chroms, by = c("gene_id" = "ref_gene_id")) |>
      dplyr::mutate(temp_observation = NA_character_) |>
      dplyr::left_join(gene_ortholog_summary, by = c("gene_id" = "ref_gene_id")) |>
      dplyr::mutate(
        num_ortholog_hits = dplyr::coalesce(num_ortholog_hits, 0L),
        all_comp_chromosomes = purrr::map(all_comp_chromosomes, function(x) {
          if (is.null(x)) return(character(0)) else return(x)
        })
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        temp_observation = dplyr::case_when(
          num_ortholog_hits == 0 ~ "NO_ORTHOLOG",
          length(all_comp_chromosomes) > 1 ~ "MULTIPLE_B_CHRS",
          length(all_comp_chromosomes) == 1 ~ paste0("ON_", primary_comp_chr),
          TRUE ~ "UNKNOWN_CLASSIFICATION"
        )
      ) |> dplyr::ungroup() |> dplyr::select(gene_id, chromosome, start, end, primary_comp_chr, temp_observation, all_comp_chromosomes, Orthogroup)

  } else {
    anc_A_file <- file.path(config$paths$results, paste0("ancestral_genome_", comp_node_daughters[1], "_", comp_node_daughters[2], ".rds"))
    anc_B_file <- file.path(config$paths$results, paste0("ancestral_genome_", comp_node_daughters[2], "_", comp_node_daughters[1], ".rds"))
    ancestral_genome_path <- if(file.exists(anc_A_file)) anc_A_file else anc_B_file
    if (!file.exists(ancestral_genome_path)) stop("Error: Ancestral genome mapping missing for structural evaluation.")

    ancestral_genome_df <- readRDS(ancestral_genome_path)

    if (config$orthology_type == "HOGs") {
      target_hog_col_name <- paste0("HOG_", parent_node)
      hog_mapping_df <- orthogroups_final |> dplyr::select(Gene_ID, Target_HOG = dplyr::sym(target_hog_col_name)) |> dplyr::filter(!is.na(Target_HOG))
      ancestral_genome_df <- ancestral_genome_df |>
        dplyr::mutate(
          Ancestor_Full_OGs = purrr::map(Ancestor_Full_Genes, map_genes_to_hogs, mapping_df = hog_mapping_df),
          Ancestor_og_count = purrr::map_int(Ancestor_Full_OGs, length)
        )
    }

    inode_og_map_raw <- ancestral_genome_df |>
      dplyr::select(ancestral_lg_name, Ancestor_Full_OGs) |>
      tidyr::pivot_longer(cols = Ancestor_Full_OGs, names_to = "source_side", values_to = "og_list") |>
      tidyr::unnest(og_list) |> dplyr::rename(Orthogroup = og_list) |> dplyr::select(ancestral_lg_name, Orthogroup) |> dplyr::distinct()

    inode_og_map_lookup <- inode_og_map_raw |>
      dplyr::group_by(Orthogroup) |> dplyr::summarise(relevant_lgs = list(unique(ancestral_lg_name)), .groups = 'drop')

    tip_gene_og_map <- ortholog_data_filtered |> dplyr::select(gene_id = ref_gene_id, Orthogroup = OG_ID) |> dplyr::distinct()

    hmm_observation_df <- ref_genes_filtered |>
      dplyr::left_join(tip_gene_og_map, by = "gene_id") |>
      dplyr::left_join(inode_og_map_lookup, by = "Orthogroup") |>
      dplyr::mutate(
        all_comp_chromosomes = dplyr::case_when(
          is.na(Orthogroup) ~ list(character(0)),
          purrr::map_lgl(relevant_lgs, is.null) ~ list(character(0)),
          TRUE ~ relevant_lgs
        ),
        lg_count = purrr::map_dbl(all_comp_chromosomes, length),
        on_lg_string = dplyr::if_else(lg_count == 1, purrr::map_chr(all_comp_chromosomes, ~ paste0("ON_", .x[1])), NA_character_),
        temp_observation = dplyr::case_when(
          lg_count == 0 ~ "NO_ORTHOLOG",
          lg_count > 1 ~ "MULTIPLE_B_CHRS",
          lg_count == 1 ~ on_lg_string,
          TRUE ~ "UNKNOWN_CLASSIFICATION"
        ),
        primary_comp_chr = ifelse(grepl("^ON_", temp_observation), sub("ON_", "", temp_observation), NA_character_)
      ) |> dplyr::select(gene_id, chromosome, start, end, primary_comp_chr, temp_observation, all_comp_chromosomes, Orthogroup)
  }

  hmm_observation_df <- hmm_observation_df |> dplyr::mutate(temp_observation = tidyr::replace_na(temp_observation, "UNKNOWN_CLASSIFICATION"))
  saveRDS(hmm_observation_df, file = hmm_observation_path)

  message("HMM sequence assignment successfully completed.")

  return(hmm_observation_path)
}
