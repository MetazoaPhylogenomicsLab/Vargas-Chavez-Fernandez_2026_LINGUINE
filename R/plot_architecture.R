#' @title Plot Ancestral Synteny Visualizations
#'
#' @description Generates Oxford Grids and high-fidelity linear/ribbon plots mapping
#' algorithmic ancestral states to extant sequences. Dynamically adapts the coordinate
#' projection based on whether the comparison involves tips or internal ancestral nodes.
#'
#' @param ref_species Character. Reference lineage identifier.
#' @param comp_species Character. Comparison lineage identifier.
#' @param comparison_type Character. Architecture of the comparison ("tip_vs_tip", "tip_vs_node", "INode_vs_INode").
#' @param parent_node Character. Target internal node being reconstructed.
#' @param config A `linguine_config` object containing path and threshold configurations.
#'
#' @return Invisible NULL. Saves plots and underlying mapping dataframes to the configured plots and results directories.
#' @export
plot_architecture <- function(ref_species, comp_species, comparison_type, parent_node, config) {

  message("\n--- Generating Matrix Topologies & Displays ---")

  anc_file <- file.path(config$paths$results, paste0("ancestral_genome_", ref_species, "_", comp_species, ".rds"))
  if (!file.exists(anc_file)) stop("Error: Ancestral model missing: ", anc_file)
  ancestral_genome <- readRDS(anc_file)

  # Centralized Orthology Dictionary
  ortho_path <- file.path(config$paths$processed_data, paste0("ortholog_data_", config$orthology_type, ".rds"))
  id_col_name <- if (config$orthology_type == "HOGs") paste0("HOG_", parent_node) else "Orthogroup"
  orthogroups_final <- readRDS(ortho_path)
  if (!id_col_name %in% colnames(orthogroups_final)) id_col_name <- "Orthogroup"

  hog_dict <- orthogroups_final |>
    dplyr::select(gene_id = Gene_ID, HOG = !!dplyr::sym(id_col_name)) |>
    dplyr::filter(!is.na(HOG)) |> dplyr::distinct()

  # Helper function to dynamically find the TRUE raw input files for nodes
  get_node_raw_file <- function(node_name) {
    tree_path <- file.path(config$paths$raw_data, config$tree_filename)
    if (!file.exists(tree_path)) tree_path <- file.path(dirname(config$paths$raw_data), config$tree_filename)
    if (file.exists(tree_path)) {
      species_tree <- ape::read.tree(tree_path)
      tips <- species_tree$tip.label
      idx <- which(species_tree$node.label == node_name) + length(tips)
      if(length(idx) > 0) {
        children <- species_tree$edge[species_tree$edge[, 1] == idx, 2]
        d_labels <- sapply(children, function(x) if (x <= length(tips)) tips[x] else species_tree$node.label[x - length(tips)])
        p1 <- file.path(config$paths$results, paste0("ancestral_genome_", d_labels[1], "_", d_labels[2], ".rds"))
        p2 <- file.path(config$paths$results, paste0("ancestral_genome_", d_labels[2], "_", d_labels[1], ".rds"))
        if (file.exists(p1)) return(p1)
        if (file.exists(p2)) return(p2)
      }
    }
    fb <- file.path(config$paths$results, paste0("ancestral_genome_", node_name, ".rds"))
    if(file.exists(fb)) return(fb)
    stop("Error: Could not locate raw ancestral file for node: ", node_name)
  }

  # ==============================================================================
  # 1. INode vs INode (Square Cloud Oxford Grid)
  # ==============================================================================
  if (comparison_type == "INode_vs_INode") {
    anc_A_file <- get_node_raw_file(ref_species)
    anc_B_file <- get_node_raw_file(comp_species)

    anc_A <- readRDS(anc_A_file)
    anc_B <- readRDS(anc_B_file)

    # Extract raw genes and bridge them into the PARENT node's HOG space
    df_A <- anc_A |> dplyr::select(A_lg = ancestral_lg_name, genes = Ancestor_Full_Genes) |>
      tidyr::unnest(genes) |> dplyr::inner_join(hog_dict, by = c("genes" = "gene_id")) |>
      dplyr::select(A_lg, HOG) |> dplyr::distinct()

    df_B <- anc_B |> dplyr::select(B_lg = ancestral_lg_name, genes = Ancestor_Full_Genes) |>
      tidyr::unnest(genes) |> dplyr::inner_join(hog_dict, by = c("genes" = "gene_id")) |>
      dplyr::select(B_lg, HOG) |> dplyr::distinct()

    plot_data <- dplyr::inner_join(df_A, df_B, by = "HOG")

    # Get colors from the generated ancestral object
    parent_map <- ancestral_genome |> dplyr::select(Parent_LG = ancestral_lg_name, Ancestor_Full_OGs) |>
      tidyr::unnest(Ancestor_Full_OGs) |> dplyr::rename(HOG = Ancestor_Full_OGs) |> dplyr::distinct()

    plot_data <- plot_data |> dplyr::left_join(parent_map, by = "HOG") |>
      dplyr::mutate(Parent_LG = tidyr::replace_na(Parent_LG, "Unassigned"))

    if(nrow(plot_data) > 0) {
      a_levels <- sort(unique(df_A$A_lg))
      b_levels <- sort(unique(df_B$B_lg))

      # Sort so Unassigned is drawn first (in the back)
      plot_data <- plot_data |>
        dplyr::arrange(Parent_LG != "Unassigned") |>
        dplyr::mutate(x_center = match(A_lg, a_levels), y_center = match(B_lg, b_levels)) |>
        dplyr::mutate(x_coord = x_center + stats::runif(dplyr::n(), min = -0.4, max = 0.4),
                      y_coord = y_center + stats::runif(dplyr::n(), min = -0.4, max = 0.4))

      lg_tags <- sort(unique(plot_data$Parent_LG))
      lg_tags <- lg_tags[lg_tags != "Unassigned"]
      custom_colors <- setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, length(lg_tags))), "Set1"))(length(lg_tags)), lg_tags)
      custom_colors["Unassigned"] <- "#D1CFCF"

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_coord, y = y_coord, color = Parent_LG)) +
        ggplot2::geom_point(alpha = 0.6, size = 1.5, stroke = 0) +
        ggplot2::scale_color_manual(values = custom_colors) +
        ggplot2::scale_x_continuous(breaks = 1:length(a_levels), labels = a_levels, expand = c(0.02, 0)) +
        ggplot2::scale_y_continuous(breaks = 1:length(b_levels), labels = b_levels, expand = c(0.02, 0)) +
        ggplot2::geom_vline(xintercept = (1:length(a_levels)) + 0.5, color = "grey90", linewidth = 0.3) +
        ggplot2::geom_hline(yintercept = (1:length(b_levels)) + 0.5, color = "grey90", linewidth = 0.3) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y = ggplot2::element_text(size = 8), panel.grid = ggplot2::element_blank(),
          legend.position = "right", legend.title = ggplot2::element_text(size = 10, face = "bold"),
          legend.text = ggplot2::element_text(size = 6)
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1))) +
        ggplot2::labs(title = paste("Ancestral Oxford Grid:", ref_species, "vs", comp_species),
                      x = paste(ref_species, "Linkage Groups"), y = paste(comp_species, "Linkage Groups"), color = paste(parent_node, "LG"))

      ggplot2::ggsave(file.path(config$paths$plots, paste0("ancestral_vs_ancestral_HOG_dotplot_", ref_species, "_", comp_species, ".png")), plot = p, width=15, height=15)
      message("Ancestral positional visualization deployed.")
    }
    return(invisible(NULL))
  }

  # ==============================================================================
  # 2. Tip vs Node (Cloud Columns Oxford Grid)
  # ==============================================================================
  if (comparison_type != "tip_vs_tip" && comparison_type != "INode_vs_INode") {

    tip_genes_raw <- readRDS(file.path(config$paths$processed_data, paste0(ref_species, "_genes_df.rds")))
    tip_chr_sizes <- readRDS(file.path(config$paths$processed_data, paste0(ref_species, "_chromosome_sizes.rds"))) |>
      dplyr::filter(chromosome_length_bp > config$min_chromosome_length_bp) |> dplyr::arrange(ref_chromosome) |>
      dplyr::mutate(cum_y_offset = dplyr::lag(cumsum(as.numeric(chromosome_length_bp)), default = 0), y_mid = cum_y_offset + chromosome_length_bp / 2)

    tip_genes <- tip_genes_raw |> dplyr::mutate(chromosome = as.character(chromosome)) |>
      dplyr::inner_join(tip_chr_sizes |> dplyr::select(chromosome = ref_chromosome, cum_y_offset), by = "chromosome") |>
      dplyr::mutate(y_coord = start + cum_y_offset)

    anc_B_file <- get_node_raw_file(comp_species)
    anc_B <- readRDS(anc_B_file)

    node_hogs <- anc_B |> dplyr::select(Node_LG = ancestral_lg_name, A_genes, B_genes) |>
      dplyr::rowwise() |> dplyr::mutate(gene_id = list(unique(c(unlist(A_genes), unlist(B_genes))))) |>
      dplyr::ungroup() |> dplyr::select(Node_LG, gene_id) |> tidyr::unnest(gene_id) |>
      dplyr::inner_join(hog_dict, by = "gene_id") |> dplyr::select(Node_LG, HOG) |> dplyr::distinct()

    parent_colors <- ancestral_genome |> dplyr::select(Parent_LG = ancestral_lg_name, Ancestor_Full_OGs) |>
      tidyr::unnest(Ancestor_Full_OGs) |> dplyr::rename(HOG = Ancestor_Full_OGs) |> dplyr::distinct()

    plot_data <- tip_genes |> dplyr::inner_join(hog_dict, by = "gene_id") |>
      dplyr::inner_join(node_hogs, by = "HOG", relationship = "many-to-many") |>
      dplyr::left_join(parent_colors, by = "HOG", relationship = "many-to-many") |>
      dplyr::mutate(Parent_LG = tidyr::replace_na(Parent_LG, "Unassigned"))

    if(nrow(plot_data) > 0) {
      node_levels <- sort(unique(plot_data$Node_LG))

      plot_data <- plot_data |>
        dplyr::arrange(Parent_LG != "Unassigned") |>
        dplyr::mutate(node_index = match(Node_LG, node_levels)) |>
        dplyr::mutate(x_coord = node_index + stats::runif(dplyr::n(), min = -0.4, max = 0.4))

      x_breaks <- data.frame(Node_LG = node_levels, x_mid = 1:length(node_levels), x_max = (1:length(node_levels)) + 0.5)

      lg_tags <- sort(unique(plot_data$Parent_LG))
      lg_tags <- lg_tags[lg_tags != "Unassigned"]
      custom_colors <- setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, length(lg_tags))), "Set1"))(length(lg_tags)), lg_tags)
      custom_colors["Unassigned"] <- "#D1CFCF"

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_coord, y = y_coord, color = Parent_LG)) +
        ggplot2::geom_point(alpha = 0.6, size = 1.5, stroke = 0) +
        ggplot2::scale_color_manual(values = custom_colors) +
        ggplot2::scale_x_continuous(breaks = x_breaks$x_mid, labels = x_breaks$Node_LG, expand = c(0.01, 0)) +
        ggplot2::scale_y_continuous(breaks = tip_chr_sizes$y_mid, labels = tip_chr_sizes$ref_chromosome, expand = c(0.01, 0), trans = "reverse") +
        ggplot2::geom_vline(xintercept = x_breaks$x_max, color = "grey80", linetype = "dotted") +
        ggplot2::geom_hline(yintercept = tip_chr_sizes$cum_y_offset, color = "grey80", linetype = "dotted") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y = ggplot2::element_text(size = 6), panel.grid = ggplot2::element_blank(),
          legend.position = "right", legend.title = ggplot2::element_text(size = 10, face = "bold"),
          legend.text = ggplot2::element_text(size = 6)
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 3))) +
        ggplot2::labs(title = paste("Diagnostic Oxford Grid:", ref_species, "vs", comp_species),
                      x = paste("True", comp_species, "Linkage Groups"), y = paste(ref_species, "Physical Chromosomes"), color = paste("Generated", parent_node, "LG"))

      ggplot2::ggsave(file.path(config$paths$plots, paste0(ref_species, "_vs_", comp_species, "_synteny_plot.png")), plot = p, width = 20, height = 15, dpi = 300)
    }
    return(invisible(NULL))
  }

  # ==============================================================================
  # 3. Tip vs Tip (Physical Extant Data Prep)
  # ==============================================================================
  if (comparison_type == "tip_vs_tip") {
    ref_genes_raw <- readRDS(file.path(config$paths$processed_data, paste0(ref_species, "_genes_df.rds"))) |> dplyr::mutate(chromosome = as.character(chromosome))
    ref_chr_sizes <- readRDS(file.path(config$paths$processed_data, paste0(ref_species, "_chromosome_sizes.rds"))) |> dplyr::mutate(ref_chromosome = as.character(ref_chromosome))
    comp_genes_raw <- readRDS(file.path(config$paths$processed_data, paste0(comp_species, "_genes_df.rds"))) |> dplyr::mutate(chromosome = as.character(chromosome))
    comp_chr_sizes <- readRDS(file.path(config$paths$processed_data, paste0(comp_species, "_chromosome_sizes.rds"))) |> dplyr::mutate(ref_chromosome = as.character(ref_chromosome))

    orthogroups_df <- orthogroups_final |> dplyr::filter(Species %in% c(ref_species, comp_species)) |>
      dplyr::select(Orthogroup = !!dplyr::sym(id_col_name), species = Species, gene_id = Gene_ID)

    og_mapping_data <- ancestral_genome |> dplyr::select(ancestral_lg_name, A_lg, B_lg, A_orthogroups, B_orthogroups) |>
      tidyr::unnest_longer(A_orthogroups, indices_include = FALSE) |> tidyr::unnest_longer(B_orthogroups, indices_include = FALSE) |>
      dplyr::filter(A_orthogroups == B_orthogroups) |> dplyr::mutate(A_lg = flatten_lg_list(A_lg), B_lg = flatten_lg_list(B_lg)) |>
      dplyr::select(ancestral_lg_name, orthogroup = A_orthogroups, ancestral_lg_A = A_lg, ancestral_lg_B = B_lg) |>
      dplyr::mutate(found_in_both_species = TRUE)

    unassigned_ogs_df <- dplyr::tibble(orthogroup = setdiff(unique(orthogroups_df$Orthogroup), unique(og_mapping_data$orthogroup)),
                                       ancestral_lg_A = NA_character_, ancestral_lg_B = NA_character_, found_in_both_species = FALSE)
    og_mapping_data_complete <- dplyr::bind_rows(og_mapping_data, unassigned_ogs_df)

    og_assignments_with_final_tags <- og_mapping_data_complete |> dplyr::distinct(orthogroup, ancestral_lg_name, found_in_both_species) |>
      dplyr::group_by(orthogroup) |> dplyr::summarise(n_total_lg = dplyr::n_distinct(ancestral_lg_name, na.rm = TRUE), n_matched_in_both = sum(found_in_both_species), assigned_lg_name = list(unique(ancestral_lg_name[found_in_both_species])), .groups = 'drop') |>
      dplyr::mutate(assignment_tag = dplyr::case_when(n_matched_in_both == 1 ~ purrr::map_chr(assigned_lg_name, 1, .default = NA_character_), n_matched_in_both > 1 ~ "multi_LG", n_matched_in_both == 0 & n_total_lg > 0 ~ "complex", TRUE ~ "no_matches")) |>
      dplyr::select(orthogroup, assignment_tag) |> dplyr::filter(!assignment_tag %in% c("multi_LG", "no_matches"))

    og_master_list <- og_assignments_with_final_tags |> dplyr::filter(!assignment_tag %in% c("no_matches", "complex")) |> dplyr::arrange(assignment_tag, orthogroup) |> dplyr::mutate(x_position = dplyr::row_number())
    LG_boundaries <- og_master_list |> dplyr::group_by(assignment_tag) |> dplyr::summarise(end_position = max(x_position), .groups = 'drop')
    LG_labels <- og_master_list |> dplyr::group_by(assignment_tag) |> dplyr::summarise(midpoint = mean(x_position), .groups = 'drop')

    gene_to_og_mapping <- orthogroups_df |> dplyr::select(gene_id = gene_id, orthogroup = Orthogroup)

    ref_genes_with_tags <- ref_genes_raw |> dplyr::left_join(gene_to_og_mapping, by = "gene_id") |> dplyr::left_join(og_assignments_with_final_tags, by = "orthogroup") |>
      dplyr::mutate(assignment_tag = dplyr::if_else(is.na(assignment_tag), "no_info", assignment_tag)) |> dplyr::select(chromosome_A = chromosome, start_A = start, dplyr::everything())

    spA_for_join <- ref_genes_with_tags |>
      dplyr::filter(!is.na(orthogroup)) |>
      dplyr::select(gene_id_A = gene_id, orthogroup, chromosome_A, start_A, end_A = end, assignment_tag) |>
      dplyr::mutate(chromosome_A = as.character(chromosome_A))

    comp_genes_for_join <- comp_genes_raw |> dplyr::left_join(gene_to_og_mapping, by = "gene_id") |>
      dplyr::filter(!is.na(orthogroup)) |>
      dplyr::select(gene_id_B = gene_id, orthogroup, chromosome_B = chromosome, start_B = start, end_B = end) |>
      dplyr::mutate(chromosome_B = as.character(chromosome_B))

    spA_vs_spB_data <- spA_for_join |>
      dplyr::inner_join(comp_genes_for_join, by = "orthogroup", relationship = "many-to-many")

    # Unified Unassigned category for valid genes that were too complex
    spA_vs_spB_data <- spA_vs_spB_data |>
      dplyr::mutate(assignment_tag = dplyr::if_else(assignment_tag %in% c("no_info", "multi_LG", "no_matches", "complex"), "Unassigned", assignment_tag))

    # Dynamic Coloring Map
    syntenic_lg_tags <- unique(spA_vs_spB_data$assignment_tag)
    syntenic_lg_tags <- syntenic_lg_tags[syntenic_lg_tags != "Unassigned"]
    lg_colors <- c("Unassigned"="#D1CFCF")
    lg_colors_ideogram <- c("Unassigned"="d1cfcf")
    if(length(syntenic_lg_tags) > 0){
      syntenic_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, length(syntenic_lg_tags))), "Set1"))(length(syntenic_lg_tags))
      lg_colors <- c(setNames(syntenic_colors, sort(syntenic_lg_tags)), lg_colors)
      lg_colors_ideogram <- tolower(gsub("#", "", lg_colors))
    }

    # Rendering
    if (nrow(spA_vs_spB_data) > 0) {
      plot_data_A <- spA_vs_spB_data |> dplyr::filter(!is.na(start_A))
      dotplotA <- generate_synteny_plot(ref_species, comp_species, ref_chr_sizes, "gene_id_A", "chromosome_A", "start_A", "end_A", plot_data_A, og_master_list, LG_boundaries, LG_labels, config$min_chromosome_length_bp, lg_colors, config)
      ggplot2::ggsave(file.path(config$paths$plots, paste0(ref_species, "_vs_", comp_species, "_dotplot.png")), plot = dotplotA, width = 45, height = 30, units = "in", dpi = 300)

      plot_data_B <- spA_vs_spB_data |> dplyr::filter(!is.na(start_B))
      dotplotB <- generate_synteny_plot(comp_species, ref_species, comp_chr_sizes, "gene_id_B", "chromosome_B", "start_B", "end_B", plot_data_B, og_master_list, LG_boundaries, LG_labels, config$min_chromosome_length_bp, lg_colors, config)
      ggplot2::ggsave(file.path(config$paths$plots, paste0(comp_species, "_vs_", ref_species, "_dotplot.png")), plot = dotplotB, width = 45, height = 30, units = "in", dpi = 300)

      chr_lens_A <- ref_chr_sizes |> dplyr::filter(chromosome_length_bp > config$min_chromosome_length_bp) |> dplyr::mutate(cumulative_length = dplyr::lag(cumsum(chromosome_length_bp), default = 0), midpoint = cumulative_length + chromosome_length_bp/2) |> dplyr::filter(!is.na(cumulative_length))
      chr_lens_B <- comp_chr_sizes |> dplyr::filter(chromosome_length_bp > config$min_chromosome_length_bp) |> dplyr::mutate(cumulative_length = dplyr::lag(cumsum(chromosome_length_bp), default = 0), midpoint = cumulative_length + chromosome_length_bp/2) |> dplyr::filter(!is.na(cumulative_length))

      plot_data_AvB <- spA_vs_spB_data |>
        dplyr::inner_join(chr_lens_A |> dplyr::select(chromosome_A = ref_chromosome, cum_A = cumulative_length), by = "chromosome_A") |>
        dplyr::inner_join(chr_lens_B |> dplyr::select(chromosome_B = ref_chromosome, cum_B = cumulative_length), by = "chromosome_B") |>
        dplyr::mutate(cum_start_A = start_A + cum_A, cum_start_B = start_B + cum_B) |> dplyr::filter(!is.na(cum_start_A) & !is.na(cum_start_B))

      if(nrow(plot_data_AvB) > 0) {

        plot_data_AvB <- plot_data_AvB |>
          dplyr::arrange(assignment_tag != "Unassigned") |>
          dplyr::mutate(is_noise = ifelse(assignment_tag == "Unassigned", "noise", "signal"))

        dotplot_AvB <- ggplot2::ggplot(plot_data_AvB, ggplot2::aes(x = cum_start_A, y = cum_start_B, color = assignment_tag, size = is_noise, alpha = is_noise)) +
          ggplot2::geom_point(stroke = 0) +
          ggplot2::scale_color_manual(values = lg_colors) +
          ggplot2::scale_size_manual(values = c("noise" = 0.4, "signal" = 1.5), guide = "none") +
          ggplot2::scale_alpha_manual(values = c("noise" = 0.2, "signal" = 0.8), guide = "none") +
          ggplot2::theme_minimal() +
          ggplot2::geom_vline(data = chr_lens_A, ggplot2::aes(xintercept = cumulative_length + chromosome_length_bp), color = "grey", linetype = "dotted") +
          ggplot2::geom_hline(data = chr_lens_B, ggplot2::aes(yintercept = cumulative_length + chromosome_length_bp), color = "grey", linetype = "dotted") +
          ggplot2::scale_x_continuous(breaks = chr_lens_A$midpoint, labels = chr_lens_A$ref_chromosome, expand = c(0,0)) +
          ggplot2::scale_y_continuous(breaks = chr_lens_B$midpoint, labels = chr_lens_B$ref_chromosome, expand = c(0,0)) +
          ggplot2::labs(x = ref_species, y = comp_species, color = "Ancestral LG") +
          ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")

        ggplot2::ggsave(file.path(config$paths$plots, paste0(ref_species, "_vs_", comp_species, "_synteny_plot.png")), plot = dotplot_AvB, width = 20, height = 20, units = "in", dpi = 300)
      }

      # RIdeogram
      karyotype_A <- ref_chr_sizes |> dplyr::filter(chromosome_length_bp > config$min_chromosome_length_bp) |> dplyr::select(Chr = ref_chromosome, End = chromosome_length_bp) |> dplyr::mutate(Start=1, species=ref_species, fill="969696", size=12, color="252525") |> dplyr::select(Chr, Start, End, fill, species, size, color)
      karyotype_B <- comp_chr_sizes |> dplyr::filter(chromosome_length_bp > config$min_chromosome_length_bp) |> dplyr::select(Chr = ref_chromosome, End = chromosome_length_bp) |> dplyr::mutate(Start=1, species=comp_species, fill="969696", size=12, color="252525") |> dplyr::select(Chr, Start, End, fill, species, size, color)
      karyo_df <- dplyr::bind_rows(karyotype_A, karyotype_B)
      mapping_A <- karyotype_A |> dplyr::mutate(index_A = as.integer(dplyr::row_number())) |> dplyr::select(Chr, index_A)
      mapping_B <- karyotype_B |> dplyr::mutate(index_B = as.integer(dplyr::row_number())) |> dplyr::select(Chr, index_B)

      filtered_ogs <- plot_data_AvB |> dplyr::count(orthogroup) |> dplyr::filter(n <= config$thresholds$ribbon_max_links) |> dplyr::pull(orthogroup)

      # Explicitly filter out Unassigned from Ribbon mapping to avoid chaos
      link_data <- plot_data_AvB |>
        dplyr::filter(orthogroup %in% filtered_ogs, assignment_tag != "Unassigned") |>
        dplyr::left_join(mapping_A, by=c("chromosome_A"="Chr")) |> dplyr::left_join(mapping_B, by=c("chromosome_B"="Chr")) |>
        dplyr::mutate(Species_1 = index_A, Start_1 = start_A, End_1 = end_A, Species_2 = index_B, Start_2 = start_B, End_2 = end_B, fill = lg_colors_ideogram[assignment_tag]) |>
        dplyr::select(Species_1, Start_1, End_1, Species_2, Start_2, End_2, fill) |> dplyr::filter(!is.na(Species_1) & !is.na(Species_2) & !is.na(fill))

      if(nrow(link_data) > 0) {
        svg_path <- file.path(config$paths$plots, paste0(ref_species, "_vs_", comp_species, "_ribbon.svg"))
        png_path <- file.path(config$paths$plots, paste0(ref_species, "_vs_", comp_species, "_ribbon.png"))
        RIdeogram::ideogram(karyotype = as.data.frame(karyo_df), synteny = as.data.frame(link_data), width = 1000, output = svg_path)
        if(requireNamespace("rsvg", quietly = TRUE)) rsvg::rsvg_png(svg_path, png_path, width = 2000)
      }
    }
  }
  message("Visual Topologies Successfully Rendered.")
}
