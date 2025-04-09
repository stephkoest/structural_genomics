# Annotating gene phylogenies with iTOL visualization

This repository contains a script to generate iTOL-compatible annotation datasets for gene phylogenies used in the manuscript:

**Koestlbacher et al., 2025. _Inference of eukaryotic complexity in Asgard archaea using structural modeling_**

The script generates datasets to utilize the Interactive Tree of Life (iTOL) tool for visualization. For more information on iTOL, see:
Letunic I, Bork P. Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool. _Nucleic Acids Res._ 2024;52(W1):W78-W82. doi: 10.1093/nar/gkae268.

## Overview

The script `generate_input_data_iTol.py` takes as input:
- a gene tree (Newick format),
- optional rooting leaves (one or more),
- and metadata tables for Archaea, Bacteria, Eukaryotes, and Asgard archaea.

It outputs:
- a reformatted Newick tree,
- iTOL annotation datasets for domain-based branch coloring, renaming, and paralog labeling.

The script utilizes the Environment for Tree Exploration (ETE 3) toolkit for phylogenetic analysis. For more information on ETE 3, see:
Huerta-Cepas J, Serra F, Bork P. ETE 3: Reconstruction, analysis, and visualization of phylogenomic data. _Mol Biol Evol._ 2016;33(6):1635-1638. doi: 10.1093/molbev/msw046.

## File structure

```
.
├── generate_input_data_iTol.py
└── phylogeny_annotation_files_ITOL/
    ├── ar53_metadata_r207.qscore.family_representative.csv
    ├── bac120_metadata_r207.qscore.family_representative.csv
    ├── Euk5FinalSet.adjust.busco.euk5_tree_abbrev.csv
    └── Asgard_DB_230420.tsv
```

## Requirements

- Python ≥ 3.7
- `ete3`
- `pandas`

## Example usage

The following commands were used to prepare visualizations in our manuscript:

- **COMMD**  
  ```bash
  python generate_input_data_iTol.py \
      -t all_homs1.COMMD.linsi.gappyout.iqtree_modelfinder.treefile \
      -r CAPOWC009071 CAPOWC009884
  ```

- **Umf1**  
  ```bash
  python generate_input_data_iTol.py \
      -t all_homs2_rg.Ufm1.linsi.bmge.iqtree_modelfinder.treefile \
      -r HeGe_bin89_00901 HeNj__B7_G17__00516
  ```

- **CINP**  
  ```bash
  python generate_input_data_iTol.py \
      -t all_homs2.CINPL.PROMALS3D.bmge.iqtree_modelfinder.treefile \
      -r Thor__Bin_478_01378 Thor_M8_58_Bin_368__01268
  ```

- **MVP shoulder**  
  ```bash
  python generate_input_data_iTol.py \
      -t all_homs1.MVP_shoulder.linsi.bmge.iqtree_modelfinder.treefile \
      -r Loki__LW55_83__00113 Bact_Cyanobacteria_Cyanobacteriia_RS_GCF_016446395.1_NZ_CP051168.1_3151
  ```

## Output

Each run creates:
- `*.reformatted`: a cleaned and optionally rerooted Newick tree
- `*.reformatted.iTOL_domain.dataset.txt`: branch colors by domain
- `*.reformatted.iTOL_labels.dataset.txt`: updated taxon names
- `*.reformatted.iTOL_paralogshapes.dataset.txt`: markers for paralogs

## Notes

- Rooting was used for inspection purposes; **unrooted phylogenies** were used in the manuscript.
- Metadata tables are expected in predefined formats (see the `phylogeny_annotation_files_ITOL/` folder for examples).
- iTOL links in the paper refer to visualizations based on these output files.

---
