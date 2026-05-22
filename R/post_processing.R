library(dplyr)
library(tidyr)
library(purrr)

# Add chrom name toggle, either no names, names as 1, 2, 3 with an output table of the equivalences or the names as they are
# Automatically find all files, only require run directory
# Automatically find the column that will be plotted from the HOGs or allow the user to select which column to plot by using the tree


plot_rideogram_exact <- function(species_order, processed_dir, ancestral_rds_path, ortholog_rds_path, id_col_name = "HOG_N3", output_filename = "Ideogram.svg", min_chr_size = 5000000, resolve_plot_multimapping = "random") {

  message("\n--- Building plot ---")

  # --- 1. Load Ancestral LGs ---
  anc_data <- readRDS(ancestral_rds_path)
  lg_mapping <- anc_data |>
    dplyr::select(ancestral_lg_name, Ancestor_Full_Genes) |>
    tidyr::unnest(Ancestor_Full_Genes) |>
    dplyr::rename(gene_id = Ancestor_Full_Genes, LG = ancestral_lg_name) |>
    dplyr::filter(!is.na(LG), LG != "Unassigned") |>
    dplyr::distinct(gene_id, LG)

  # --- 2. Load Orthogroups ---
  ortho_df <- readRDS(ortholog_rds_path)
  hog_dict <- ortho_df |>
    dplyr::select(gene_id = Gene_ID, orthogroup = !!dplyr::sym(id_col_name)) |>
    dplyr::filter(!is.na(orthogroup)) |>
    dplyr::distinct()

  chr_layouts <- list()
  genes_list <- list()

  # --- 3. Pixel-Perfect Coordinate Mapping ---
  svg_width <- 1800
  margin_left <- 150
  margin_right <- 50
  gap_px <- 12

  for (i in seq_along(species_order)) {
    sp <- species_order[i]

    chr_file <- file.path(processed_dir, paste0(sp, "_chromosome_sizes.rds"))
    chrs <- readRDS(chr_file) |>
      dplyr::rename(seq_id = ref_chromosome, length = chromosome_length_bp) |>
      dplyr::filter(length >= min_chr_size)

    gene_file <- file.path(processed_dir, paste0(sp, "_genes_df.rds"))
    genes <- readRDS(gene_file) |>
      dplyr::rename(seq_id = chromosome) |>
      dplyr::inner_join(lg_mapping, by = "gene_id") |>
      dplyr::inner_join(hog_dict, by = "gene_id") |>
      dplyr::filter(seq_id %in% chrs$seq_id)

    dom_lg <- genes |>
      dplyr::count(seq_id, LG) |>
      dplyr::group_by(seq_id) |>
      dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) |>
      dplyr::arrange(LG)

    chrs <- chrs |> dplyr::mutate(seq_id = factor(seq_id, levels = dom_lg$seq_id)) |> dplyr::arrange(seq_id)

    total_bp <- sum(chrs$length)
    available_px <- (svg_width - margin_left - margin_right) - ((nrow(chrs) - 1) * gap_px)
    bp_to_px <- available_px / total_bp

    chrs <- chrs |>
      dplyr::mutate(
        width_px = length * bp_to_px,
        svg_x_start = margin_left + cumsum(dplyr::lag(width_px + gap_px, default = 0)),
        svg_x_end = svg_x_start + width_px,
        species = sp
      )
    chr_layouts[[sp]] <- chrs

    genes <- genes |>
      dplyr::inner_join(chrs |> dplyr::select(seq_id, svg_x_start), by = "seq_id") |>
      dplyr::mutate(
        gene_px_start = svg_x_start + (start * bp_to_px),
        gene_px_end = svg_x_start + (end * bp_to_px),
        gene_px_end = ifelse(gene_px_end - gene_px_start < 1.0, gene_px_start + 1.0, gene_px_end),
        species = sp
      )
    genes_list[[sp]] <- genes
  }

  chrs_df <- dplyr::bind_rows(chr_layouts)
  all_genes_df <- dplyr::bind_rows(genes_list)

  # =======================================================================
  # --- 3.5 MULTI-MAPPING PLOT RESOLUTION (USER CONTROLLED) ---
  if (resolve_plot_multimapping == "random") {
    message("Plot Resolution: Randomly assigning multi-mapped OGs to a single LG for visual consistency...")
    set.seed(42)

    resolved_lgs <- all_genes_df |>
      dplyr::select(orthogroup, LG) |>
      dplyr::distinct() |>
      dplyr::group_by(orthogroup) |>
      dplyr::slice_sample(n = 1) |>
      dplyr::ungroup() |>
      dplyr::rename(Resolved_LG = LG)

    all_genes_df <- all_genes_df |>
      dplyr::left_join(resolved_lgs, by = "orthogroup") |>
      dplyr::select(-LG) |>
      dplyr::rename(LG = Resolved_LG)

  } else if (resolve_plot_multimapping == "drop") {
    message("Plot Resolution: Dropping multi-mapped OGs to strictly enforce 1-to-1 visual orthology...")

    # Identify OGs that map to exactly 1 unique LG globally across the plot
    strict_ogs <- all_genes_df |>
      dplyr::group_by(orthogroup) |>
      dplyr::summarise(n_lgs = dplyr::n_distinct(LG), .groups = "drop") |>
      dplyr::filter(n_lgs == 1) |>
      dplyr::pull(orthogroup)

    # Filter the plot data to only include those pure families
    all_genes_df <- all_genes_df |>
      dplyr::filter(orthogroup %in% strict_ogs)

  } else if (resolve_plot_multimapping == "keep") {
    message("Plot Resolution: Keeping all multi-mappings (Note: conflicting translocations will not draw ribbons).")
    # Do nothing, leave the data exactly as it is.
  } else {
    stop("CRITICAL ERROR: Unknown resolve_plot_multimapping strategy. Use 'random', 'drop', or 'keep'.")
  }
  # =======================================================================

  # --- 4. Canvas & Styling ---
  unique_lgs <- sort(unique(all_genes_df$LG))
  palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(max(3, min(9, length(unique_lgs))), "Set1"))(length(unique_lgs))
  lg_colors <- setNames(palette, unique_lgs)

  svg_height <- 200 + (length(species_order) * 200) + 150
  chr_height <- 12
  y_spacing <- (svg_height - 300) / max(1, length(species_order) - 1)
  if(length(species_order) == 2) y_spacing <- 300

  svg_lines <- c(
    '<?xml version="1.0" standalone="no"?>',
    sprintf('<svg width="%.0f" height="%.0f" xmlns="http://www.w3.org/2000/svg">', svg_width, svg_height),
    '<rect width="100%" height="100%" fill="#FFFFFF"/>',
    '<style>',
    '  .chr { fill: #E8E8E8; stroke: #888888; stroke-width: 1.0; }',
    '  .label { font-family: Arial, sans-serif; font-size: 16px; font-weight: bold; fill: #333333; }',
    '  .legend-text { font-family: Arial, sans-serif; font-size: 14px; fill: #333333; }',
    '</style>'
  )

  # --- 5. Draw Individual Gene Ribbons ---
  message("Weaving gene-by-gene ribbons...")

  for (i in 1:(length(species_order) - 1)) {
    sp1 <- species_order[i]
    sp2 <- species_order[i + 1]

    sp1_genes <- all_genes_df |> dplyr::filter(species == sp1) |> dplyr::select(orthogroup, LG, start1 = gene_px_start, end1 = gene_px_end)
    sp2_genes <- all_genes_df |> dplyr::filter(species == sp2) |> dplyr::select(orthogroup, LG, start2 = gene_px_start, end2 = gene_px_end)

    sp_links <- dplyr::inner_join(sp1_genes, sp2_genes, by = c("orthogroup", "LG"), relationship = "many-to-many") |>
      dplyr::group_by(orthogroup) |>
      dplyr::slice_head(n = 2) |>
      dplyr::ungroup()

    y_top_chr <- 80 + ((i - 1) * y_spacing)
    y_start <- y_top_chr + chr_height
    y_end <- 80 + (i * y_spacing)

    # Staggered control points for smooth, non-clipping Bezier curves
    y_ctrl_1 <- y_start + (y_end - y_start) * 0.35
    y_ctrl_2 <- y_start + (y_end - y_start) * 0.65

    for (j in 1:nrow(sp_links)) {
      row <- sp_links[j, ]
      color <- unname(lg_colors[as.character(row$LG)])

      path <- sprintf(
        '<path d="M %.2f %.2f L %.2f %.2f C %.2f %.2f, %.2f %.2f, %.2f %.2f L %.2f %.2f C %.2f %.2f, %.2f %.2f, %.2f %.2f Z" fill="%s" fill-opacity="0.5" />',
        row$start1, y_start,
        row$end1, y_start,
        row$end1, y_ctrl_1, row$end2, y_ctrl_2, row$end2, y_end,
        row$start2, y_end,
        row$start2, y_ctrl_2, row$start1, y_ctrl_1, row$start1, y_start,
        color
      )
      svg_lines <- c(svg_lines, path)
    }
  }

  # --- 6. Draw Chromosomes on top of ribbons ---
  message("Rendering chromosomes...")
  for (i in seq_along(species_order)) {
    sp <- species_order[i]
    y_chr <- 80 + ((i - 1) * y_spacing)

    text_line <- sprintf('<text x="%.2f" y="%.2f" class="label" text-anchor="end">%s</text>',
                         margin_left - 20, y_chr + 11, sp)
    svg_lines <- c(svg_lines, text_line)

    sp_chrs <- chrs_df |> dplyr::filter(species == sp)
    for (j in 1:nrow(sp_chrs)) {
      row <- sp_chrs[j, ]
      rect_line <- sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" class="chr" rx="2" ry="2" />',
                           row$svg_x_start, y_chr, row$width_px, chr_height)
      svg_lines <- c(svg_lines, rect_line)
    }
  }

  # --- 7. Draw The Legend ---
  message("Building color legend...")
  legend_start_y <- 80 + ((length(species_order) - 1) * y_spacing) + 100
  legend_item_width <- 120
  items_per_row <- floor((svg_width - margin_left - margin_right) / legend_item_width)

  for (k in seq_along(unique_lgs)) {
    lg_name <- unique_lgs[k]
    color <- unname(lg_colors[as.character(lg_name)])

    row_idx <- floor((k - 1) / items_per_row)
    col_idx <- (k - 1) %% items_per_row

    x_pos <- margin_left + (col_idx * legend_item_width)
    y_pos <- legend_start_y + (row_idx * 30)

    square <- sprintf('<rect x="%.2f" y="%.2f" width="14" height="14" fill="%s" rx="2" ry="2" />',
                      x_pos, y_pos, color)
    label <- sprintf('<text x="%.2f" y="%.2f" class="legend-text" text-anchor="start">%s</text>',
                     x_pos + 20, y_pos + 12, lg_name)

    svg_lines <- c(svg_lines, square, label)
  }

  svg_lines <- c(svg_lines, '</svg>')

  # --- 8. Save ---
  writeLines(svg_lines, output_filename)
  message(paste0("Success! Plot clone saved as: ", output_filename))

  if (requireNamespace("rsvg", quietly = TRUE)) {
    png_file <- sub("\\.svg$", ".png", output_filename)
    rsvg::rsvg_png(output_filename, png_file, width = svg_width * 2)
  }
}


# 1. Define your 3 species in order
species_to_plot <- c("SCAV", "HMAN", "PCOS", "BLOB", "DELA", "NNAJ", "EAND", "MVUL")

# 2. Define the exact paths to your data
processed_dir <- "~/Documents/GitHub/getLinkageGroups/runs/Annelida_Rpackage_min_chr_size_5000000bp/processed_data/"
rds_file <- "~/Documents/GitHub/getLinkageGroups/runs/Annelida_Rpackage_min_chr_size_5000000bp/results/ancestral_genome_N1_N2.rds"
ortho_file <- "~/Documents/GitHub/getLinkageGroups/runs/Annelida_Rpackage_min_chr_size_5000000bp/processed_data/ortholog_data_HOGs.rds"

# 3. Launch the function!
plot_rideogram_exact(
  species_order = species_to_plot,
  processed_dir = processed_dir,
  ancestral_rds_path = rds_file,
  ortholog_rds_path = ortho_file,
  id_col_name = "Orthogroup",
  output_filename = "Perfect_Ternary_RIdeogram.svg"
)

# 1. Define your 3 species in order
species_to_plot <- c("HCON", "OTIP", "CBRI", "CNIO", "CREM", "CELE", "CINO", "PEXS", "PRPA", "BOKI", "SRAT")

# 2. Define the exact paths to your data
processed_dir <- "~/Documents/GitHub/getLinkageGroups/runs/Nematoda_Rpackage_min_chr_size_4500000bp/processed_data/"
rds_file <- "~/Documents/GitHub/getLinkageGroups/runs/Nematoda_Rpackage_min_chr_size_4500000bp/results/ancestral_genome_N3_N4.rds"
ortho_file <- "~/Documents/GitHub/getLinkageGroups/runs/Nematoda_Rpackage_min_chr_size_4500000bp/processed_data/ortholog_data_HOGs.rds"

# 3. Launch the function!
plot_rideogram_exact(
  species_order = species_to_plot,
  processed_dir = processed_dir,
  ancestral_rds_path = rds_file,
  ortholog_rds_path = ortho_file,
  id_col_name = "HOG_N1",
  output_filename = "Nematoda_RIdeogram.svg"
)
