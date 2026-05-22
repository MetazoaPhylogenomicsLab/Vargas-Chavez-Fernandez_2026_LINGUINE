#' @title Reconstruct Ancestral Node Topology
#'
#' @description Executes topological graph clustering to merge, resolve, and polarize
#' ancestral linkage groups using bipartite matching and an outgroup.
#'
#' @param ref_species Character. Reference lineage.
#' @param comp_species Character. Comparison lineage.
#' @param comparison_type Character. Architecture of the comparison.
#' @param parent_node Character. Target internal node being reconstructed.
#' @param species_tree phylo object. The full phylogenetic tree for outgroup polarization.
#' @param config A `linguine_config` object.
#'
#' @return Path to the saved ancestral genome RDS file.
#' @export
reconstruct_node <- function(ref_species, comp_species, comparison_type, parent_node, species_tree, config) {

  message(sprintf("\n--- Initiating Topological Graph Modeling & Reconstruction (%s) ---", parent_node))

  get_anc_file <- function(node_name) {
    tips <- species_tree$tip.label
    idx <- which(species_tree$node.label == node_name) + length(tips)
    children <- species_tree$edge[species_tree$edge[, 1] == idx, 2]
    d_labels <- sapply(children, function(x) if (x <= length(tips)) tips[x] else species_tree$node.label[x - length(tips)])

    p1 <- file.path(config$paths$results, paste0("ancestral_genome_", d_labels[1], "_", d_labels[2], ".rds"))
    p2 <- file.path(config$paths$results, paste0("ancestral_genome_", d_labels[2], "_", d_labels[1], ".rds"))

    if (file.exists(p1)) return(p1)
    if (file.exists(p2)) return(p2)
    stop("CRITICAL ERROR: Could not find ancestral genome file for node: ", node_name)
  }

  # ----------------------------------------------------------------------------
  # 1. Load Structural Foundations
  # ----------------------------------------------------------------------------
  if (comparison_type != "INode_vs_INode") {
    file_A_vs_B <- file.path(config$paths$results, paste0("broad_syntenic_regions_", ref_species, "_ref_", comp_species, "_comp.rds"))
    broad_regions_A_vs_B <- readRDS(file_A_vs_B)

    if (comparison_type == "tip_vs_tip") {
      file_B_vs_A <- file.path(config$paths$results, paste0("broad_syntenic_regions_", comp_species, "_ref_", ref_species, "_comp.rds"))
      broad_regions_B_vs_A <- readRDS(file_B_vs_A)
    } else {
      anc_path <- get_anc_file(comp_species)
      broad_regions_B_vs_A <- create_synthetic_broad_regions(readRDS(anc_path))
    }
  } else {
    anc_A <- readRDS(get_anc_file(ref_species))
    anc_B <- readRDS(get_anc_file(comp_species))

    if(config$orthology_type == "HOGs") {
      orthogroups_full <- readRDS(file.path(config$paths$processed_data, "ortholog_data_HOGs.rds"))
      hog_col_name <- if(parent_node == "N0") "Orthogroup" else paste0("HOG_", parent_node)

      if (!hog_col_name %in% colnames(orthogroups_full)) {
        if (parent_node == "N0") {
          hog_col_name <- "Orthogroup"
        } else {
          stop(paste("CRITICAL ERROR: Column", hog_col_name, "is missing from the OrthoFinder data."))
        }
      }

      hog_mapping_df <- orthogroups_full |> dplyr::select(Gene_ID, Target_HOG = dplyr::sym(hog_col_name)) |> dplyr::filter(!is.na(Target_HOG))
      anc_A <- anc_A |> dplyr::mutate(Ancestor_Full_OGs = lapply(Ancestor_Full_Genes, function(x) map_genes_to_hogs(unlist(x), hog_mapping_df)))
      anc_B <- anc_B |> dplyr::mutate(Ancestor_Full_OGs = lapply(Ancestor_Full_Genes, function(x) map_genes_to_hogs(unlist(x), hog_mapping_df)))
    }
    broad_regions_A_vs_B <- create_synthetic_broad_regions(anc_A)
    broad_regions_B_vs_A <- create_synthetic_broad_regions(anc_B)
  }

  # ----------------------------------------------------------------------------
  # 2. Execute Ubiquity / Noise Base Filter
  # ----------------------------------------------------------------------------
  ortho_path <- file.path(config$paths$processed_data, paste0("ortholog_data_", config$orthology_type, ".rds"))
  orthogroups_final <- readRDS(ortho_path)

  id_col_name <- if(config$orthology_type == "HOGs") paste0("HOG_", parent_node) else "Orthogroup"

  if (config$orthology_type == "HOGs" && !id_col_name %in% colnames(orthogroups_final)) {
    if (parent_node == "N0") {
      id_col_name <- "Orthogroup"
    } else {
      stop(paste("CRITICAL ERROR: Column", id_col_name, "is missing from the OrthoFinder data."))
    }
  }

  orthogroups_df <- orthogroups_final |>
    dplyr::select(Orthogroup = dplyr::sym(id_col_name), species = Species, gene_id = Gene_ID) |>
    dplyr::filter(!is.na(Orthogroup))

  temp_A <- map_orthogroups_to_blocks(broad_regions_A_vs_B, orthogroups_df, "broad_ref_gene_ids")
  temp_B <- map_orthogroups_to_blocks(broad_regions_B_vs_A, orthogroups_df, "broad_ref_gene_ids")

  og_lg_counts <- dplyr::bind_rows(
    temp_A |> tidyr::unnest(orthogroups_in_block) |> dplyr::group_by(orthogroups_in_block) |> dplyr::summarize(n_lgs = dplyr::n_distinct(linkage_group_name)),
    temp_B |> tidyr::unnest(orthogroups_in_block) |> dplyr::group_by(orthogroups_in_block) |> dplyr::summarize(n_lgs = dplyr::n_distinct(linkage_group_name))
  ) |> dplyr::group_by(orthogroups_in_block) |> dplyr::summarize(max_lg_span = max(n_lgs))

  abundance_threshold <- stats::quantile(og_lg_counts$max_lg_span, probs = config$thresholds$og_abundance/100, na.rm = TRUE)
  filtered_og_list <- og_lg_counts |> dplyr::filter(max_lg_span <= abundance_threshold) |> dplyr::pull(orthogroups_in_block)

  orthogroups_df <- orthogroups_df |> dplyr::filter(Orthogroup %in% filtered_og_list)
  all_orthogroups <- unique(orthogroups_df$Orthogroup)

  # ----------------------------------------------------------------------------
  # 3. Paralogy Consolidation
  # ----------------------------------------------------------------------------
  broad_A_og <- clean_broad_regions(merge_fragmented_lgs(map_orthogroups_to_blocks(broad_regions_A_vs_B, orthogroups_df, "broad_ref_gene_ids")))
  broad_B_og <- clean_broad_regions(merge_fragmented_lgs(map_orthogroups_to_blocks(broad_regions_B_vs_A, orthogroups_df, "broad_ref_gene_ids")))

  paralogy_A <- detect_paralogy_signals(broad_A_og, ref_species, all_orthogroups, config)
  paralogy_B <- detect_paralogy_signals(broad_B_og, comp_species, all_orthogroups, config)

  dup_events <- infer_duplication_events(paralogy_A, paralogy_B, broad_A_og, broad_B_og)

  pairs_A <- dplyr::bind_rows(dup_events$ancestral_dups |> dplyr::select(block1_id=A_lg1, block2_id=A_lg2), dup_events$lineage_A |> dplyr::select(block1_id=A_lg1, block2_id=A_lg2))
  pairs_B <- dplyr::bind_rows(dup_events$ancestral_dups |> dplyr::select(block1_id=B_lg1, block2_id=B_lg2), dup_events$lineage_B |> dplyr::select(block1_id=B_lg1, block2_id=B_lg2))

  collapsed_lg_A <- merge_lg_fragments(collapse_lg_pairs(broad_A_og, pairs_A))
  collapsed_lg_B <- merge_lg_fragments(collapse_lg_pairs(broad_B_og, pairs_B))

  # --- PARALOGY DIAGNOSTICS & WGD ALERTS ---
  blocks_in_dups_A <- length(unique(c(pairs_A$block1_id, pairs_A$block2_id)))
  if (nrow(broad_A_og) > 0) {
    perc_A <- (blocks_in_dups_A / nrow(broad_A_og)) * 100
    message(sprintf("\n[Paralogy Analysis] %s:", ref_species))
    message(sprintf("  -> Initial genomic blocks: %d", nrow(broad_A_og)))
    message(sprintf("  -> Blocks involved in duplication: %d (%.1f%%)", blocks_in_dups_A, perc_A))
    message(sprintf("  -> Total blocks after paralog collapse: %d", nrow(collapsed_lg_A)))
    if (perc_A >= 60) {
      message(sprintf("  *** HIGH PARALOGY DETECTED: The massive duplication signal (%.1f%% of blocks) strongly suggests the %s lineage underwent a Whole Genome Duplication (WGD). ***", perc_A, ref_species))
    }
  }

  blocks_in_dups_B <- length(unique(c(pairs_B$block1_id, pairs_B$block2_id)))
  if (nrow(broad_B_og) > 0) {
    perc_B <- (blocks_in_dups_B / nrow(broad_B_og)) * 100
    message(sprintf("\n[Paralogy Analysis] %s:", comp_species))
    message(sprintf("  -> Initial genomic blocks: %d", nrow(broad_B_og)))
    message(sprintf("  -> Blocks involved in duplication: %d (%.1f%%)", blocks_in_dups_B, perc_B))
    message(sprintf("  -> Total blocks after paralog collapse: %d", nrow(collapsed_lg_B)))
    if (perc_B >= 60) {
      message(sprintf("  *** HIGH PARALOGY DETECTED: The massive duplication signal (%.1f%% of blocks) strongly suggests the %s lineage underwent a Whole Genome Duplication (WGD). ***", perc_B, comp_species))
    }
  }
  # -----------------------------------------

  gene_to_chr_A <- broad_regions_A_vs_B |> dplyr::select(chromosome, broad_ref_gene_ids) |> tidyr::unnest(broad_ref_gene_ids) |> dplyr::rename(gene_id = broad_ref_gene_ids) |> dplyr::distinct()
  gene_to_chr_B <- broad_regions_B_vs_A |> dplyr::select(chromosome, broad_ref_gene_ids) |> tidyr::unnest(broad_ref_gene_ids) |> dplyr::rename(gene_id = broad_ref_gene_ids) |> dplyr::distinct()

  collapsed_lg_A <- collapsed_lg_A |>
    dplyr::mutate(physical_chr = purrr::map_chr(genes_in_block, function(g) {
      chrs <- gene_to_chr_A |> dplyr::filter(gene_id %in% g) |> dplyr::pull(chromosome)
      if(length(chrs) == 0) return("Unknown")
      names(sort(table(chrs), decreasing = TRUE))[1]
    })) |>
    dplyr::filter(physical_chr != "Unknown") |>
    dplyr::mutate(linkage_group_name = paste0(physical_chr, "::", linkage_group_name))

  collapsed_lg_B <- collapsed_lg_B |>
    dplyr::mutate(physical_chr = purrr::map_chr(genes_in_block, function(g) {
      chrs <- gene_to_chr_B |> dplyr::filter(gene_id %in% g) |> dplyr::pull(chromosome)
      if(length(chrs) == 0) return("Unknown")
      names(sort(table(chrs), decreasing = TRUE))[1]
    })) |>
    dplyr::filter(physical_chr != "Unknown") |>
    dplyr::mutate(linkage_group_name = paste0(physical_chr, "::", linkage_group_name))

  # ----------------------------------------------------------------------------
  # 4. Outgroup Prep & Bipartite Graph
  # ----------------------------------------------------------------------------
  species_C_name <- get_outgroup_species(species_tree, ref_species, comp_species)
  common_ancestor_hog_col <- NULL

  if (!is.null(species_C_name)) {
    idx_ref <- get_node_index(species_tree, ref_species)
    idx_C <- get_node_index(species_tree, species_C_name)
    if (!is.null(idx_ref) && !is.null(idx_C)) {
      common_node_idx <- ape::getMRCA(species_tree, c(idx_ref, idx_C))
      if (!is.null(species_tree$node.label)) {
        common_node_label <- species_tree$node.label[common_node_idx - length(species_tree$tip.label)]
      } else {
        common_node_label <- paste0("N", common_node_idx)
      }
      common_ancestor_hog_col <- paste0("HOG_", common_node_label)

      if (!common_ancestor_hog_col %in% colnames(orthogroups_final)) {
        common_ancestor_hog_col <- "Orthogroup"
      }
    }
  }

  C_lg_map <- NULL
  if (!is.null(species_C_name)) {
    is_C_node <- !(species_C_name %in% species_tree$tip.label)
    if (is_C_node) {
      anc_outgroup_path <- get_anc_file(species_C_name)
      outgroup_data <- readRDS(anc_outgroup_path)
      C_genes_list <- unlist(outgroup_data$Ancestor_Full_Genes)
      node_gene_map <- outgroup_data |> dplyr::select(LG_name = ancestral_lg_name, Ancestor_Full_Genes) |> tidyr::unnest(Ancestor_Full_Genes) |> dplyr::rename(gene_id = Ancestor_Full_Genes)
      C_lg_map <- orthogroups_final |> dplyr::filter(Gene_ID %in% C_genes_list) |> dplyr::select(Target_HOG = !!dplyr::sym(common_ancestor_hog_col), gene_id = Gene_ID) |> dplyr::inner_join(node_gene_map, by = "gene_id") |> dplyr::select(Target_HOG, LG_name) |> dplyr::filter(!is.na(Target_HOG)) |> dplyr::distinct()
    } else {
      C_content <- load_genome_content(species_C_name, is_C_node, config$paths$processed_data, config$paths$results)
      if (!is.null(C_content)) {
        C_genes <- C_content$data |> dplyr::select(gene_id, chromosome)
        C_lg_map <- orthogroups_final |> dplyr::filter(Species == species_C_name) |> dplyr::select(Target_HOG = !!dplyr::sym(common_ancestor_hog_col), gene_id = Gene_ID) |> dplyr::inner_join(C_genes, by = "gene_id") |> dplyr::select(Target_HOG, LG_name = chromosome) |> dplyr::filter(!is.na(Target_HOG)) |> dplyr::distinct()
      }
    }
  }

  hog_to_parent_map <- NULL
  if (!is.null(common_ancestor_hog_col)) {
    hog_to_parent_map <- orthogroups_final |> dplyr::filter(!is.na(!!dplyr::sym(id_col_name))) |> dplyr::select(Local_HOG = !!dplyr::sym(id_col_name), Parent_HOG = !!dplyr::sym(common_ancestor_hog_col)) |> dplyr::filter(!is.na(Parent_HOG)) |> dplyr::distinct()
  }

  num_A <- nrow(collapsed_lg_A)
  num_B <- nrow(collapsed_lg_B)
  pairs_grid <- expand.grid(i = 1:num_A, j = 1:num_B)
  total_orthos <- length(unique(all_orthogroups))

  results_list <- purrr::map_dfr(1:nrow(pairs_grid), function(k) {
    i <- pairs_grid$i[k]; j <- pairs_grid$j[k]
    og_A <- unlist(collapsed_lg_A$orthogroups_in_block[i])
    og_B <- unlist(collapsed_lg_B$orthogroups_in_block[j])
    shared <- length(intersect(og_A, og_B))
    if (shared == 0) return(NULL)
    mat <- matrix(c(shared, length(og_A)-shared, length(og_B)-shared, total_orthos - (length(og_A) + length(og_B) - shared)), nrow=2)
    if (min(mat) < 0) return(NULL)
    ft <- stats::fisher.test(mat, alternative = "greater")
    return(dplyr::tibble(A_lg = collapsed_lg_A$linkage_group_name[i], B_lg = collapsed_lg_B$linkage_group_name[j], p_val = ft$p.value, odds = ft$estimate))
  })

  dyn_odds <- max(1.5, min(10, ((num_A + num_B)/2) / 2))
  edges <- results_list |> dplyr::mutate(p_adj = stats::p.adjust(p_val, method = "BH")) |>
    dplyr::filter(p_adj < config$thresholds$fisher_p_lgs, odds > dyn_odds) |>
    dplyr::mutate(node_A = paste0("A|", A_lg), node_B = paste0("B|", B_lg)) |> dplyr::select(node_A, node_B, odds)

  if(nrow(edges) == 0) {
    edges <- results_list |> dplyr::mutate(p_adj = stats::p.adjust(p_val, method="BH")) |> dplyr::filter(p_adj < 1e-02, odds > 1) |> dplyr::mutate(node_A = paste0("A|", A_lg), node_B = paste0("B|", B_lg)) |> dplyr::select(node_A, node_B, odds)
  }

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  missing_nodes <- setdiff(c(paste0("A|", collapsed_lg_A$linkage_group_name), paste0("B|", collapsed_lg_B$linkage_group_name)), igraph::V(g)$name)
  if(length(missing_nodes) > 0) g <- igraph::add_vertices(g, nv = length(missing_nodes), name = missing_nodes)

  # ----------------------------------------------------------------------------
  # 5. Component Polarization & Reassembly
  # ----------------------------------------------------------------------------
  components <- igraph::components(g)
  mem_list <- split(names(components$membership), components$membership)
  final_ancestral_lgs_list <- list()

  for (comp_id in names(mem_list)) {
    nodes <- mem_list[[comp_id]]
    lgs_A_in_comp <- sub("^A\\|", "", nodes[grep("^A\\|", nodes)])
    lgs_B_in_comp <- sub("^B\\|", "", nodes[grep("^B\\|", nodes)])

    genes_A_pool <- collapsed_lg_A |> dplyr::filter(linkage_group_name %in% lgs_A_in_comp) |> dplyr::pull(genes_in_block) |> unlist() |> unique()
    genes_B_pool <- collapsed_lg_B |> dplyr::filter(linkage_group_name %in% lgs_B_in_comp) |> dplyr::pull(genes_in_block) |> unlist() |> unique()
    ogs_A_pool <- collapsed_lg_A |> dplyr::filter(linkage_group_name %in% lgs_A_in_comp) |> dplyr::pull(orthogroups_in_block) |> unlist() |> unique()
    ogs_B_pool <- collapsed_lg_B |> dplyr::filter(linkage_group_name %in% lgs_B_in_comp) |> dplyr::pull(orthogroups_in_block) |> unlist() |> unique()
    all_comp_ogs <- unique(c(ogs_A_pool, ogs_B_pool))

    decision <- "FUSE"
    partitions <- list(); partition_map <- NULL
    is_1to1_consensus <- (length(lgs_A_in_comp) == 1 && length(lgs_B_in_comp) == 1)

    if (!is_1to1_consensus && !is.null(C_lg_map) && length(all_comp_ogs) > 0) {
      if (!is.null(hog_to_parent_map)) {
        relevant_parent_hogs <- hog_to_parent_map |> dplyr::filter(Local_HOG %in% all_comp_ogs) |> dplyr::pull(Parent_HOG) |> unique()
        hits <- C_lg_map |> dplyr::filter(Target_HOG %in% relevant_parent_hogs)

        if (nrow(hits) > 0) {
          dominant_LGs <- hits |> dplyr::count(LG_name) |> dplyr::mutate(freq = n/sum(n)) |> dplyr::filter(freq > config$thresholds$outgroup_dominance)

          if (nrow(dominant_LGs) > 1) {
            decision <- "SPLIT"
            partitions <- dominant_LGs$LG_name
            partition_map <- hits |> dplyr::filter(LG_name %in% partitions) |> dplyr::group_by(LG_name) |> dplyr::summarize(parent_hogs = list(unique(Target_HOG)), .groups = "drop")
          }
        }
      }
    }

    if (decision == "FUSE") {
      final_ancestral_lgs_list[[length(final_ancestral_lgs_list)+1]] <- dplyr::tibble(
        A_lg = list(lgs_A_in_comp), B_lg = list(lgs_B_in_comp),
        A_orthogroups = list(unique(ogs_A_pool)), B_orthogroups = list(unique(ogs_B_pool)),
        A_genes = list(unique(genes_A_pool)), B_genes = list(unique(genes_B_pool)),
        Ancestor_Full_OGs = list(unique(all_comp_ogs)), Ancestor_Full_Genes = list(unique(c(genes_A_pool, genes_B_pool))),
        A_og_count = length(unique(ogs_A_pool)), B_og_count = length(unique(ogs_B_pool)),
        Decision = "FUSE"
      )
    } else {
      candidates <- list()
      for (i in 1:nrow(partition_map)) {
        part_lg_name <- partition_map$LG_name[i]
        part_parent_hogs <- unlist(partition_map$parent_hogs[i])
        target_local_hogs <- hog_to_parent_map |> dplyr::filter(Parent_HOG %in% part_parent_hogs, Local_HOG %in% all_comp_ogs) |> dplyr::pull(Local_HOG)

        get_best_match <- function(lgs_in_comp, collapsed_df, target_ogs) {
          overlap <- collapsed_df |> dplyr::filter(linkage_group_name %in% lgs_in_comp) |> dplyr::rowwise() |> dplyr::mutate(hits = length(intersect(unlist(orthogroups_in_block), target_ogs))) |> dplyr::ungroup() |> dplyr::filter(hits > 0) |> dplyr::mutate(prop = hits / sum(hits)) |> dplyr::arrange(dplyr::desc(hits))
          if(nrow(overlap) > 0 && overlap$prop[1] > config$thresholds$parent_assignment) return(overlap$linkage_group_name[1]) else return("Ambiguous")
        }
        candidates[[i]] <- dplyr::tibble(scav_part = part_lg_name, part_ogs = list(target_local_hogs), parentA = get_best_match(lgs_A_in_comp, collapsed_lg_A, target_local_hogs), parentB = get_best_match(lgs_B_in_comp, collapsed_lg_B, target_local_hogs))
      }
      cand_df <- dplyr::bind_rows(candidates)

      supported_A_LGs <- unique(setdiff(cand_df$parentA, "Ambiguous"))
      supported_B_LGs <- unique(setdiff(cand_df$parentB, "Ambiguous"))
      winner <- "FUSE"
      if (length(supported_A_LGs) > 1 && length(supported_B_LGs) <= 1) winner <- "A"
      else if (length(supported_B_LGs) > 1 && length(supported_A_LGs) <= 1) winner <- "B"
      else if (length(supported_A_LGs) > 1 && length(supported_B_LGs) > 1) winner <- if(length(supported_A_LGs) >= length(supported_B_LGs)) "A" else "B"

      if (winner == "FUSE") {

        final_ancestral_lgs_list[[length(final_ancestral_lgs_list)+1]] <- dplyr::tibble(
          A_lg = list(lgs_A_in_comp), B_lg = list(lgs_B_in_comp),
          A_orthogroups = list(unique(ogs_A_pool)), B_orthogroups = list(unique(ogs_B_pool)),
          A_genes = list(unique(genes_A_pool)), B_genes = list(unique(genes_B_pool)),
          Ancestor_Full_OGs = list(unique(all_comp_ogs)), Ancestor_Full_Genes = list(unique(c(genes_A_pool, genes_B_pool))),
          A_og_count = length(unique(ogs_A_pool)), B_og_count = length(unique(ogs_B_pool)),
          Decision = "FUSE_NO_SUPPORT"
        )
      } else {

        primary_LGs <- if(winner == "A") supported_A_LGs else supported_B_LGs
        all_LGs_in_cluster <- if(winner == "A") lgs_A_in_comp else lgs_B_in_comp
        orphan_LGs <- setdiff(all_LGs_in_cluster, primary_LGs)

        merge_map <- list()
        for(p in primary_LGs) merge_map[[p]] <- c(p)
        source_df <- if(winner == "A") collapsed_lg_A else collapsed_lg_B

        for(orphan in orphan_LGs) {
          orphan_ogs <- source_df |> dplyr::filter(linkage_group_name == orphan) |> dplyr::pull(orthogroups_in_block) |> unlist()
          best_match_primary <- NULL; max_overlap <- 0
          for(i in 1:nrow(cand_df)) {
            owner <- if(winner == "A") cand_df$parentA[i] else cand_df$parentB[i]
            if(owner == "Ambiguous" || !(owner %in% primary_LGs)) next
            hits <- length(intersect(orphan_ogs, unlist(cand_df$part_ogs[i])))
            if(hits > max_overlap) { max_overlap <- hits; best_match_primary <- owner }
          }
          if(!is.null(best_match_primary)) merge_map[[best_match_primary]] <- c(merge_map[[best_match_primary]], orphan) else merge_map[[orphan]] <- c(orphan)
        }

        for (primary_key in names(merge_map)) {
          lg_group <- merge_map[[primary_key]]
          lg_data <- source_df |> dplyr::filter(linkage_group_name %in% lg_group)
          target_ogs <- unique(unlist(lg_data$orthogroups_in_block))
          valid_genes_A <- orthogroups_df |> dplyr::filter(gene_id %in% genes_A_pool, Orthogroup %in% target_ogs) |> dplyr::pull(gene_id)
          valid_genes_B <- orthogroups_df |> dplyr::filter(gene_id %in% genes_B_pool, Orthogroup %in% target_ogs) |> dplyr::pull(gene_id)

          final_ancestral_lgs_list[[length(final_ancestral_lgs_list)+1]] <- dplyr::tibble(
            A_lg = list(if(winner == "A") lg_group else lgs_A_in_comp), B_lg = list(if(winner == "B") lg_group else lgs_B_in_comp),
            A_orthogroups = list(unique(intersect(ogs_A_pool, target_ogs))), B_orthogroups = list(unique(intersect(ogs_B_pool, target_ogs))),
            A_genes = list(unique(valid_genes_A)), B_genes = list(unique(valid_genes_B)),
            Ancestor_Full_OGs = list(unique(target_ogs)), Ancestor_Full_Genes = list(unique(c(valid_genes_A, valid_genes_B))),
            A_og_count = length(unique(valid_genes_A)), B_og_count = length(unique(valid_genes_B)),
            Decision = paste0("SPLIT_", winner, "_", primary_key, (if(length(lg_group)>1) "_RESCUED" else ""))
          )
        }
      }
    }
  }

  # Initial compilation of the ancestral genome
  ancestral_genome <- dplyr::bind_rows(final_ancestral_lgs_list) |>
    dplyr::mutate(Ancestor_og_count = purrr::map_int(Ancestor_Full_OGs, length)) |>
    dplyr::arrange(dplyr::desc(Ancestor_og_count)) |>
    dplyr::mutate(
      ancestral_lg_name = paste0("LG_", stringr::str_pad(dplyr::row_number(), width = 3, pad = "0")),
      Shared_og_count = Ancestor_og_count,
      is_strict_1to1 = (purrr::map_int(A_lg, length) == 1 & purrr::map_int(B_lg, length) == 1 & !grepl("SPLIT", Decision)),
      A_original_regions = list(dplyr::tibble()),
      B_original_regions = list(dplyr::tibble())
    )

  # ==============================================================================
  # 6. Apply User-Defined Noise & Paralogy Filters
  # ==============================================================================

  # A. The Proportional Size Filter (Microchromosome Sweeper)
  total_unique_ogs <- length(unique(unlist(ancestral_genome$Ancestor_Full_OGs)))
  min_required_ogs <- ceiling(total_unique_ogs * config$thresholds$min_lg_fraction)

  ancestral_genome <- ancestral_genome |>
    dplyr::filter(Ancestor_og_count >= min_required_ogs)

  # B. The Multi-mapped Resolution Filter
  resolution_strategy <- config$thresholds$resolve_multimapped

  if (resolution_strategy == "drop") {
    message("Applying STRICT multimapping resolution: Dropping duplicated orthogroups.")

    og_occurrences <- table(unlist(ancestral_genome$Ancestor_Full_OGs))
    strict_ogs <- names(og_occurrences[og_occurrences == 1])

    ancestral_genome <- ancestral_genome |>
      dplyr::mutate(Ancestor_Full_OGs = purrr::map(Ancestor_Full_OGs, ~intersect(.x, strict_ogs)))

  } else if (resolution_strategy == "random") {
    message("Applying RANDOM multimapping resolution: Randomly assigning duplicated orthogroups to a single LG.")

    # Use a fixed seed for reproducible pipeline runs
    set.seed(42)

    # Unnest, randomly pick one LG per OG, and renest
    resolved_ogs <- ancestral_genome |>
      dplyr::select(ancestral_lg_name, Ancestor_Full_OGs) |>
      tidyr::unnest(Ancestor_Full_OGs) |>
      dplyr::group_by(Ancestor_Full_OGs) |>
      dplyr::slice_sample(n = 1) |>
      dplyr::ungroup() |>
      dplyr::group_by(ancestral_lg_name) |>
      dplyr::summarise(Resolved_OGs = list(Ancestor_Full_OGs), .groups = "drop")

    ancestral_genome <- ancestral_genome |>
      dplyr::left_join(resolved_ogs, by = "ancestral_lg_name") |>
      dplyr::mutate(
        Ancestor_Full_OGs = ifelse(purrr::map_lgl(Resolved_OGs, is.null), list(character(0)), Resolved_OGs)
      ) |>
      dplyr::select(-Resolved_OGs)

  } else if (resolution_strategy == "keep") {
    message("Applying KEEP multimapping resolution: Preserving WGD/paralogy signals across multiple LGs.")
    # We do nothing to the lists! They remain exactly as the graph algorithm built them.
  } else if (resolution_strategy == "consensus") {
    message("Applying CONSENSUS multimapping resolution: Assigning duplicated OGs to their majority LG (Ties broken randomly).")

    # 1. Unnest ONLY the genes safely
    expanded_lgs <- ancestral_genome |>
      dplyr::select(ancestral_lg_name, Ancestor_Full_Genes) |>
      tidyr::unnest(Ancestor_Full_Genes) |>
      dplyr::rename(gene_id = Ancestor_Full_Genes)

    # 2. Map genes to OGs, count the weight, and resolve ties randomly
    resolved_ogs <- expanded_lgs |>
      dplyr::inner_join(orthogroups_df |> dplyr::select(gene_id, Orthogroup), by = "gene_id") |>
      dplyr::group_by(Orthogroup, ancestral_lg_name) |>
      dplyr::summarise(weight = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(Orthogroup) |>
      # Grab the max weight, and explicitly KEEP any ties
      dplyr::slice_max(order_by = weight, n = 1, with_ties = TRUE) |>
      # Randomly pick exactly 1 from the remaining row(s)
      dplyr::slice_sample(n = 1) |>
      dplyr::ungroup() |>
      dplyr::group_by(ancestral_lg_name) |>
      dplyr::summarise(Resolved_OGs = list(Orthogroup), .groups = "drop")

    # 3. Apply the resolved list back to the ancestral genome
    ancestral_genome <- ancestral_genome |>
      dplyr::left_join(resolved_ogs, by = "ancestral_lg_name") |>
      dplyr::mutate(
        Ancestor_Full_OGs = ifelse(purrr::map_lgl(Resolved_OGs, is.null), list(character(0)), Resolved_OGs)
      ) |>
      dplyr::select(-Resolved_OGs)
    } else {
    stop("CRITICAL ERROR: Unknown resolve_multimapped strategy in config. Use 'drop', 'random', 'keep' or 'consensus'.")
  }

  # C. Final Recalculation & Cleanup
  ancestral_genome <- ancestral_genome |>
    dplyr::mutate(
      Ancestor_og_count = purrr::map_int(Ancestor_Full_OGs, length),
      Shared_og_count = Ancestor_og_count
    ) |>
    # Drop any LGs that became completely empty after filtering
    dplyr::filter(Ancestor_og_count > 0) |>
    dplyr::arrange(dplyr::desc(Ancestor_og_count)) |>
    dplyr::mutate(
      # Rename them cleanly from 1 to N based on their final curated sizes
      ancestral_lg_name = paste0("LG_", stringr::str_pad(dplyr::row_number(), width = 3, pad = "0"))
    ) |>
    dplyr::select(ancestral_lg_name, dplyr::everything())

  output_path <- file.path(config$paths$results, paste0("ancestral_genome_", ref_species, "_", comp_species, ".rds"))
  saveRDS(ancestral_genome, output_path)

  message(sprintf("Graph modeling successfully deployed. (Retained %d Curated Linkage Groups)", nrow(ancestral_genome)))

  return(output_path)
}
