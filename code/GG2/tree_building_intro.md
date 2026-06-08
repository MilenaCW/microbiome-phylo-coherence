# Phylogenetic Tree Building Introduction
**Date:** Updated Jan 2026

---

## Overview
This markdown explains the conceptual differences between **reference‑tree** approaches (using backbones from GG2 or GTDB, as is used in this project) and **naïve evolutionary model** approaches (e.g., FastTree on your sequences).

---

## Two Strategies for Building Trees

### 1) Reference‑tree approach (e.g., **Greengenes2** or **GTDB**)
**How mapping works**
- Your input sequences (ASVs/OTUs) are **placed onto a prebuilt backbone** comprised of genomes and 16S sequences.
- Multiple ASVs/OTUs can **map to the same backbone tip** (feature). In GG2, this is often intentional: closely related/identical sequences collapse to a single representative (“feature”).
- Output is a **consistent phylogeny across datasets**, because every dataset is ultimately mapped to the same backbone.

**Pros**
- Computationally efficient and scalable; mapping is far faster than de‑novo tree inference on thousands of sequences.
- Cross‑study comparability: the backbone provides a stable coordinate system for phylogenetic distances.
- Backbone tips are tied to genomes + full‑length 16S, which can improve stability of deep branches.

**Cons**
- **Granularity loss:** distinct OTUs/ASVs can collapse to a smaller set of backbone features (many‑to‑one mapping).
- The mapping step may leave some sequences **unmapped** (identity threshold, fragments, quality).

### 2) Naïve evolutionary model approach (e.g., **FastTree** on your MSA)
**How it works**
- Build a multiple sequence alignment (MSA) for your sequences and infer a tree using a substitution model (e.g., FastTree, IQ‑TREE).

**Pros**
- **Maximum granularity:** every unique sequence appears as its own tip.
- No dependency on an external backbone—full control over model & options.

**Cons**
- **Computationally heavy** for large, variable‑length 16S sets; alignment quality and trimming become critical.
- **Study‑specific phylogeny:** distances are not directly comparable across datasets built independently.
- Requires careful parameter choice (model, filtering, masking), and results may vary more across runs/setups.

---

## What GG2/GTDB Do With Multiple ASVs/OTUs Per Tip
- Mapping/placement intentionally **groups near‑identical fragments** to a single representative feature.
- Expect a distribution of **occupancy** (number of original OTUs per GG2 feature). The diagnostic script used in this project visualizes this distribution and quantifies mapping success vs. unmapped.