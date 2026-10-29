# ThreePillar

**Functional dependency, not network topology, carries anticancer drug-synergy signal.**

A reproducible MATLAB benchmark that asks, systematically, which mechanism-aware
features carry genuine drug-synergy signal once drug promiscuity is controlled.
Three independent feature families are tested on public data under strict grouped
cross-validation and a confound-matched permutation test:

1. **Static network motifs** (STRING proximity, SIGNOR feed-forward loops)
2. **Context-specific motifs** (the same motifs on expression-gated, cell-line-specific networks)
3. **Functional genetic dependency** (DepMap CRISPR Chronos scores)

The finding: static and context-specific topology do not beat a trivial
target-count baseline, and their apparent association with synergy is an artefact
of promiscuity. Functional dependency is the one family that carries a real,
confound-independent signal.

**Author:** Dr. Soutrick Das
School of Computer Science and Artificial Intelligence, SR University, Warangal, Telangana 506371, India
Email: soutrickd@gmail.com

---

## Requirements

- MATLAB R2019b or later
- Statistics and Machine Learning Toolbox (for `fitcensemble` / RUSBoost)

No other dependencies. The scripts locate their own folder, so no paths need editing.

---

## Repository layout

```
ThreePillar/
  run_all_threepillar.m        one entry point, runs the whole study
  step1_features.m             family 1: drug-target graph and static motifs
  step2_omics.m                molecular features (drug and cell-line PCA)
  step3_models.m               grouped CV: network ablation (Table 1), predictor (Table 4)
  step4_figures.m              main figures (Fig. 2 to Fig. 6)
  step5_baseline.m             representation and model comparison (Table 5)
  step6_external.m             external validation on the O'Neil screen (Fig. 6)
  step7_deep.m                 DeepSynergy-style MLP baseline
  step9_context_audit.m        family 2: context-specific motif audit (Table 2)
  step10_dependency_audit.m    family 3: functional dependency audit (Table 3)
  plot_pillars.m               three-family comparison figure (Fig. 4)
  auroc.m  auprc.m  read_named.m  read_cellminer.m  normcell.m
  normalize_drug.m  canonical_drug.m  ffl_pair_features.m      (helpers)
  data/                        place all input files here (see below)
  results/                     outputs are written here (figures in results/figures/)
```

---

## Data

All inputs are from public resources. Download them and place them in `data/`
under the exact file names below.

| File (in `data/`) | Source | Description |
|---|---|---|
| `almanac_synergy.csv` | [DrugComb](https://drugcomb.org) | NCI-ALMANAC pairs, ZIP synergy scores |
| `9606.protein.links.v12.0.txt` | [STRING v12](https://string-db.org) | human protein-protein links |
| `9606.protein.info.v12.0.txt` | [STRING v12](https://string-db.org) | protein to gene-symbol map |
| `interactions.tsv` | [DGIdb](https://dgidb.org) | drug-target interactions |
| `all_data_15_09_26.tsv` | [SIGNOR](https://signor.uniroma2.it) | directed causal signalling network |
| `nci60_expr.txt` | [CellMiner](https://discover.nci.nih.gov/cellminer) | NCI-60 composite expression |
| `oneil_synergy.csv` | [O'Neil 2016](https://doi.org/10.1158/1535-7163.MCT-15-0843) / DrugComb | external validation screen |
| `CRISPRGeneEffect.csv` | [DepMap Public 26Q3](https://depmap.org/portal/download/) | Chronos gene-effect matrix (~497 MB) |
| `Model.csv` | [DepMap Public 26Q3](https://depmap.org/portal/download/) | cell-line model metadata |

Data files are not tracked in this repository (see `.gitignore`); download them
from the sources above.

---

## How to run

Open `run_all_threepillar.m` and press **Run**. It executes the full pipeline:

```matlab
>> run_all_threepillar
```

To run only the two audit pillars after the predictor has produced
`results/omics_features.mat` once, call them directly with no arguments:

```matlab
>> step9_context_audit
>> step10_dependency_audit
>> plot_pillars
```

Reading the 497 MB DepMap matrix in `step10` is the slow part (a few minutes).

---

## What each part produces

| Manuscript item | Produced by | Result file |
|---|---|---|
| Table 1 (static ablation) | `step3_models` | console + `model_results.mat` |
| Table 2 (context audit) | `step9_context_audit` | `context_audit.mat` |
| Table 3 (dependency audit) | `step10_dependency_audit` | `dependency_audit.mat` |
| Table 4 (per-cell-line model) | `step3_models` | `model_results.mat` |
| Table 5 (baselines, MLP) | `step5_baseline`, `step7_deep` | `baseline_results.mat`, `deep_results.mat` |
| Fig. 2 (data) | `step4_figures` | `results/figures/Fig2_data.pdf` |
| Fig. 3 (confound) | `step4_figures` | `results/figures/Fig4_mechanism.pdf` |
| Fig. 4 (three families) | `plot_pillars` | `results/figures/Fig7_pillars.pdf` |
| Fig. 5 (prediction) | `step4_figures` | `results/figures/Fig3_results.pdf` |
| Fig. 6 (external) | `step6_external` | `results/figures/Fig6_external.pdf` |

Figure file names follow generation order and differ from the final manuscript
numbering; the mapping above is the reference.

---

## Expected results (console)

| Family | Outcome |
|---|---|
| 1, static motifs | do not beat the target-count baseline (leave-drug-out AUROC ~0.60); no feature survives the confound-matched test |
| 2, context motifs | within-cell AUROC 0.50 to 0.53; every feature returns q = 1.00 |
| 3, dependency | `dep_sum_min` within-cell absolute AUROC ~0.57, q = 0.001, cell+count-matched p = 0.001 |

Predictive increment of the dependency block over chemistry plus expression:
about +0.006 AUROC for unseen cell lines (significant) and about -0.003 for unseen
drugs (not significant). The signal is genuine but predictively marginal, because
cell-line expression already encodes much of the dependency context.

---

## Citation

If you use this code, please cite the accompanying paper (details to be added on
acceptance). A `CITATION.cff` file is included in this repository.

---

## License

Released under the MIT License. See `LICENSE`.
