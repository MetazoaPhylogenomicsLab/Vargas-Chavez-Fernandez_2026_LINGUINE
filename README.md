# LINGUINE <img src="man/figures/logo.png" align="left" height="130" />
**LINkage GroUps INfErence**

[![R CMD Check](https://github.com/cvargas88/LINGUINE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cvargas88/LINGUINE/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Reconstructing the chromosomal architecture of ancestral genomes is a central challenge in comparative genomics. Existing methods struggle with rapid genome evolution, extensive chromosomal rearrangements, and whole-genome duplications (WGDs) followed by gene loss. 

**LINGUINE** is a phylogeny-aware R pipeline designed to solve these limitations. Instead of relying purely on single-copy orthologs, LINGUINE leverages full orthogroups to maximize marker density. It employs a Hidden Markov Model (HMM) framework to seamlessly bridge syntenic blocks across gene loss, micro-inversions, and assembly fragmentation.

Crucially, LINGUINE is highly flexible regarding its orthology inputs. It natively ingests outputs from **OrthoFinder**, allowing users to run the pipeline using either standard flat Orthogroups (`Orthogroups.tsv`) or high-resolution Hierarchical Orthogroups (`HOGs`). 

## ⚙️ Pipeline Overview
LINGUINE reconstructs ancestral states via an iterative post-order traversal of your species tree. For every internal node, the pipeline executes the following sequence:

1. **Genomic Pre-Processing:** Parses FASTA and GFF files, dynamically filtering out fragmented micro-scaffolds based on user-defined length thresholds.
2. **Orthology Integration:** Maps global orthogroups (OGs or HOGs) to the physical chromosomal coordinates of the descendant lineages.
3. **HMM Synteny Delineation:** Uses the Viterbi algorithm to decode hidden ancestral states, delineating contiguous syntenic blocks while ignoring local noise.
4. **Bipartite Graph Reconstruction:** Projects syntenic blocks into a bipartite graph, utilizing an outgroup to polarize "FUSE" vs "SPLIT" chromosomal fusion events.
5. **Paralogy Resolution:** Explicitly tests for statistically significant paralogy enrichments, collapsing duplicated linkage groups arising from WGDs prior to ancestral finalization to prevent karyotype inflation.

## ⚠️ Critical Data Requirement
**Your GFF feature IDs must perfectly match your OrthoFinder sequence IDs!**
If your genome annotation pipeline appends suffixes (e.g., `-PA`, `_t1`) to transcripts, but your OrthoFinder dictionary lacks them, LINGUINE will fail to bridge the physical coordinates with the evolutionary orthology. Ensure your GFF `ID=` attributes and OrthoFinder `.tsv` gene names are strictly harmonized before running the pipeline.

## 🛠 Installation

You can install the development version of LINGUINE from GitHub using `devtools`:

```r
# install.packages("pak")
pak::pak("MetazoaPhylogenomicsLab/Vargas-Chavez-Fernandez_2026_LINGUINE")
```

## 🚀 Quick Start & Configuration

LINGUINE is operated entirely through a centralized configuration object. This ensures your runs are highly reproducible.

```r
library(LINGUINE)

# 1. Generate the configuration object
config <- create_linguine_config(
  dataset = "Nematoda_Analysis",
  base_dir = "/path/to/your/working/directory",
  tree_filename = "SpeciesTree_rooted_node_labels.txt",
  orthology_type = "HOGs",
  orthology_filename = "Phylogenetic_Hierarchical_Orthogroups",
  min_chromosome_length_bp = 4500000, 
  resolve_multimapped = "drop" 
)

# 2. Execute the entire pipeline
run_linguine(config)
```

### Key Configuration Parameters
The `create_linguine_config()` function accepts several arguments to tune the pipeline to your specific clade's evolutionary rate:

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `dataset` | *String* | The unique name of your run (used for generating output folders). |
| `base_dir` | *String* | Absolute path to the working directory containing your `raw_data` folder. |
| `tree_filename` | *String* | Name of your Newick species tree file (must be inside `raw_data`). |
| `orthology_type` | *String* | Either `"OGs"` for standard Orthogroups or `"HOGs"` for Hierarchical. |
| `orthology_filename` | *String* | The OrthoFinder output file (e.g., `"Orthogroups.tsv"` or the directory name containing HOG `.tsv` files). |
| `min_chromosome_length_bp`| *Numeric* | Drops physical scaffolds smaller than this threshold (in base pairs) to remove assembly noise. |
| `resolve_multimapped` | *String* | Paralogy strategy: `"drop"` (enforces strict 1-to-1 synteny), `"keep"` (preserves paralogy for WGD hunting), or `"random"` (preserves density without inflation). |
| `min_lg_fraction` | *Numeric* | The minimum fraction of total mapped orthogroups a reconstructed Linkage Group must contain to be retained (e.g., `0.01` for 1%). |
| `purity_threshold` | *Numeric* | Synteny HMM parameter. The minimum fractional purity (0.0 to 1.0) required to confidently classify a contiguous block. Default is `0.80`. |
| `paralogy_p` | *Numeric* | Adjusted p-value threshold for the Fisher's Exact Test when collapsing duplicated blocks. Default is `1e-20`. |

*(Note: Advanced HMM transition and emission probabilities can also be customized. Run `?create_linguine_config` in R for the full parameter list).*

## 📊 Pipeline Outputs
All outputs are automatically sorted into the directory defined in your `base_dir`:
* `processed_data/`: Sanitized gene coordinates, filtered chromosome sizes, and centralized orthology dictionaries.
* `intermediate_data/`: HMM emission matrices and Viterbi state classifications.
* `results/`: The finalized `.rds` objects for the Ancestral Genomes at every internal node.
* `plots/`: High-resolution PNGs and SVGs of phylogenetic trees, Oxford Grids, and physical karyotypes.

## 📖 Citation
If you use LINGUINE in your research, please cite:
> Vargas-Chávez, C., & Fernández, R. (202X). *[Paper Title]*. [Journal Name]. [DOI]
