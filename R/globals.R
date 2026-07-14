# ==============================================================================
# Script: R/globals.R
# Purpose: Global variable declarations to satisfy R CMD check for tidyverse NSE.
# ==============================================================================

# Silence R CMD check notes for tidyverse non-standard evaluation
utils::globalVariables(c(
  ".data", ":=", "A_genes", "A_lg", "A_lg1", "A_lg2", "A_orthogroups", "Ancestor_Full_Genes",
  "Ancestor_Full_OGs", "Ancestor_og_count", "Ancestral_LG", "B_genes", "B_lg", "B_lg1", "B_lg2",
  "B_orthogroups", "Chr", "Decision", "End", "End_1", "End_2", "Gene_ID", "Gene_IDs_List",
  "Gene_ID_Comp", "Gene_ID_Ref", "HOG", "LG", "LG_name", "Local_HOG", "Node_LG", "OG", "OG_ID",
  "Orthogroup", "Parent_HOG", "Parent_LG", "Relative_Order", "Resolved_LG", "Resolved_OGs", "Species", "Species_1", "Species_2",
  "Start", "Start_1", "Start_2", "Target_HOG", "aggregate", "all_comp_chromosomes",
  "ancestral_lg_name", "ancestral_lg_name_A", "assignment_tag", "block1_id", "block2_id",
  "block_end", "block_end_pos", "block_group_key", "block_id", "block_length", "block_start",
  "block_start_pos", "bound", "broad_block_key", "broad_comp_end", "broad_comp_end_initial",
  "broad_comp_gene_ids", "broad_comp_start", "broad_comp_start_initial", "broad_ref_end",
  "broad_ref_gene_ids", "broad_ref_start", "chromEnd", "chromStart", "chromosome", "chromosome_A", "chromosome_B",
  "chromosome_length_bp", "color", "colorRampPalette", "comp_block_end", "comp_block_start",
  "comp_chromosome", "comp_end", "comp_gene_id", "comp_gene_ids_in_block", "comp_start", "count_chr", "cum_A", "cum_B",
  "cum_start_A", "cum_start_B", "cum_y_offset", "cumulative_length", "cumulative_start",
  "end", "end_A", "end_B", "end_pos", "end_position", "fill", "final_synteny_classification",
  "fisher.test", "fisher_key", "found_in_both_species", "freq", "fx", "fy", "gene_id",
  "gene_px_end", "gene_px_start", "genes", "genes_in_block", "group_key", "has_comp_orthologs", "head", "i", "index_A",
  "index_B", "inferred_comp_chr", "inferred_state_blocks", "inferred_state_final",
  "inferred_state_hmm1", "is_noise", "j", "lg_count", "linkage_group_name", "local_dominant_chr",
  "max_lg_span", "max_r", "n", "n_genes", "n_lgs", "na.omit", "name", "node_A", "node_B", "node_index",
  "num_ortholog_hits", "num_ref_OGs_in_block", "num_strict_blocks_merged", "odds", "odds_ratio",
  "og_count", "og_list", "original_regions", "orthogroup", "orthogroups_in_block", "overlapping_block",
  "p.adjust", "p_adj", "p_val", "p_value", "p_value_adjusted", "paralogy_odds_ratio_threshold",
  "paralogy_p_threshold", "percentage", "physical_chr", "primary_comp_chr", "quantile",
  "raw_observations", "read.delim", "ref_block_end", "ref_block_start", "ref_chromosome",
  "ref_end", "ref_gene_id", "ref_gene_ids_in_block", "ref_gene_ids_in_strict_blocks",
  "ref_start", "relevant_orthologs", "rle_id", "sA", "sB", "seqid", "seq_id", "setNames", "size", "species",
  "start", "start_A", "start_B", "start_pos", "svg_x_start", "target_comp_chromosome", "temp_observation",
  "total_genes", "total_unique_og_count", "type", "weight", "width_px", "write.table", "x_center", "x_coord",
  "x_position", "y_center", "y_coord"
))
