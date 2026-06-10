# LINGUINE Architectural & Algorithmic Analysis Summary

**Date:** June 8, 2026

This document summarizes the codebase analysis, algorithmic review, and proposed alternative strategies for the LINGUINE (LINkage GroUps INfErence for Ancestral Genomes) R package.

## 1. Initial Codebase Review & Recommendations
The initial review of the `R/` directory identified several areas for improvement, focusing on performance, robustness, and R package best practices.

*   **Performance (Spatial Overlaps):** The pipeline currently uses `dplyr::rowwise()` and `purrr::pmap()` for spatial overlaps, resulting in $O(N \times M)$ complexity. Refactoring to use interval trees (e.g., `IRanges::findOverlaps()` or `data.table::foverlaps()`) could yield massive speedups.
*   **Robustness (Dangerous Loops):** Usage of `1:length(x)` or `1:nrow(df)` is unsafe and can cause out-of-bounds errors if objects are empty. Replace with `seq_along(x)` and `seq_len(nrow(df))`.
*   **R Package Best Practices:** 
    *   Avoid altering global state with `set.seed(42)` directly; scope it using `withr::with_seed()`.
    *   Explicitly declare `Depends: R (>= 4.1.0)` in `DESCRIPTION` for native pipe `|>` usage.
    *   Remove redundant runtime namespace checks for packages already listed in `Imports` (e.g., `igraph`).

## 2. Logical and Algorithmic Review
The core biological logic—using orthogroups, HMMs, and post-order traversal to infer macro-synteny and resolve WGDs—is sound and achieves its purpose well. However, improvements can be made:

*   **Empty Columns in Ancestral Objects:** It was noted that columns like `comp_gene_ids_in_block` become `character(0)` for ancestral nodes. This is a correct biological abstraction since ancestral nodes lack physical coordinates and are represented by orthogroups. However, storing explicit metadata (e.g., `status_flag`) instead of relying on `NA` or `"Unassigned"` would improve transparency.
*   **Methodological Improvements:**
    *   **Dynamic HMM Transitions:** Instead of uniform static transition matrices, allow transitions to scale with physical distance or recombination rates.
    *   **Weighted Bipartite Matching:** Enhance WGD collapse logic by using weighted bipartite graphs (e.g., Jaccard similarity of shared orthogroups) rather than strict percentage thresholds.
*   **Data Structure Upgrades:** Store intermediate HMM posterior probabilities, outgroup signal strengths, and alternative graph topologies to allow for confidence assessment and sensitivity analysis.

## 3. Blueprint for "LINGUINE 2.0" (Object-Oriented Redesign)
A proposed architectural redesign to modernize the pipeline without altering the core biological goals:

1.  **State Management:** Implement S3/R6 classes (`LINGUINE_Genome`, `LINGUINE_Ancestor`) to handle state polymorphically, avoiding empty columns and fragile dataframe passing.
2.  **High-Performance Engines:**
    *   `engine_interval_math.R`: Uses interval trees for instant spatial mapping.
    *   `engine_hmm.R`: Vectorized HMM observations with dynamic transitions.
    *   `engine_graph_matching.R`: Mathematical resolution of WGDs via Maximum-Weight Bipartite Matching.
3.  **Strict Typing & Caching:** Pure functions, rigorous input validation, and parameter-aware hash-based caching.

## 4. An Alternative Paradigm: Graph-Theory & Adjacencies
A theoretical exploration of a fundamentally different algorithmic strategy, shifting from "Set-Theory/Block-Matching" to "Graph-Theory/Breakpoint Analysis".

*   **Noise Reduction:** Replace sequence-based HMMs with continuous spatial density clustering like **DBSCAN** to identify syntenic diagonals and flag translocations natively.
*   **Formulating Evolution:** Track genomic **Adjacencies** (connections between genes) rather than discrete "Linkage Groups" (buckets of genes).
*   **Ancestral Reconstruction:** Use **Maximum Parsimony / Maximum Likelihood** on adjacencies to minimize the total evolutionary cost (fusions, fissions, breakages) across the entire tree simultaneously.
*   **Handling Paralogy (WGDs):** Resolving multi-chromosome orthogroups becomes an NP-hard problem. It requires algorithms like Maximum Weight Matching and Guided Genome Halving to trace parallel tracks of adjacencies and find the pre-WGD ancestor that minimizes structural breakages across the tree.

**Summary of Paradigms:**
*   **LINGUINE (Current):** Compares "bags of genes" (Macro-blocks) using majority-rule thresholds. Efficient, linear time complexity, heuristic.
*   **Adjacency Strategy (Proposed):** Analyzes "micro-context" (neighbors). Uses parsimony to pair up copies requiring the fewest structural breakages. Computationally heavier (NP-hard approximations) but mathematically rigorous.
