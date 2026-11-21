# Compartment-resolved immunity in irreversible pulpitis

Analysis code and numerical result tables for:

> **Compartment-resolved immunity in irreversible pulpitis: compositional bulk signals and a venular interface**
>
> Zhongzhen Ren, School of Stomatology, Shandong Medical and Pharmaceutical University, Yantai 264003, China

This repository accompanies the Scientific Reports submission. All analyses are secondary analyses of public data; no new sequencing or qPCR data were generated. The manuscript has **eight main display items (Figs. 1–8) and no main tables**. The four tables formerly planned as main tables are **Supplementary Tables 14–17**.

## 1. Release archive

The complete, hash-manifested release is [`pulpitis_compartment_immunity_code_20261120.zip`](pulpitis_compartment_immunity_code_20261120.zip). It contains 42 R scripts, 2 Python figure/manuscript builders, 75 CSV result tables, `MANIFEST_SHA256.csv`, `CODE_VALIDATION.md` and `package_versions.csv`.

SHA-256: `c4ece488efec9d83f403042f54284691073b8668096d3115ccef12d0189dd58b`

The `R/` and `results/` directories provide a browsable copy of the release layout. The ZIP is authoritative for the final versions, including the scripts added during the final Figure 7 and patient-level revisions.

## 2. Public datasets

| Dataset | Platform / design | Role |
|---|---|---|
| GSE77459 | Affymetrix Human Gene 2.1 ST, GPL17692; 6 healthy and 6 irreversible-pulpitis bulk samples; 27,279 genes | Discovery bulk transcriptome |
| GSE92681 | Agilent-045997 Arraystar human lncRNA microarray V3, GPL16956; 5 normal and 7 pulpitis samples; 66 marker genes annotated by 60-mer matching | Independent bulk M1/M2 sensitivity analysis |
| GSE185222 | Illumina HiSeq 4000; 2 sound/enamel-caries and 2 deep-caries teeth; 6,582 cells | Single-cell discovery and compartment localisation |
| GSE193448 | RT-PCR panel, 808 nm photobiomodulation; control, LPS and LPS + laser; 90 wells (84 targets plus 6 controls) | Descriptive supplementary PBM panel only |
| E-MTAB-14896 | Parse Evercode WT v2, NovaSeq X; 4 healthy, 4 early and 4 advanced patients; 53,407 cells and 60,721 genes | Three-stage validation and patient-level analyses |

Processed public matrices must be obtained from GEO or ArrayExpress/BioStudies. Raw expression matrices and large RDS objects are not redistributed. The scripts document how the included result tables were generated.

## 3. Directory layout

- `R/bulk/`: GSE77459 and GSE92681 differential expression, enrichment, immune scores, WGCNA, machine-learning, PBM, drug-signature and bulk M1/M2 analyses.
- `R/singlecell/`: GSE185222 preparation, annotation, composition, communication, trajectory, deconvolution and myeloid analyses.
- `R/emtab/`: E-MTAB-14896 three-stage composition, cell-type analyses, TF/PROGENy activity, communication, trajectory, robustness and patient-level effect sizes.
- `R/manuscript/`: Python figure assembly and manuscript builders. `60_build_sr_figures.py` builds main Figs. 1–6 and 8 plus Supplementary Figs. 1–28; `61_build_current_figures.py` builds Fig. 7 and Supplementary Figs. 29–33.
- `results/bulk/`: numerical bulk outputs.
- `results/singlecell/`: numerical outputs from both single-cell cohorts and single-cell-informed bulk deconvolution.

The file manifest in the release ZIP is the authoritative inventory.

## 4. Display-item mapping

### Main figures (8 figures, 0 main tables)

| Item | Content | Principal scripts / inputs |
|---|---|---|
| Fig. 1 | Bulk pulpitis transcriptome: volcano plot and Hallmark GSEA | `06_DEG_enrich_GSE77459.R`, `29_GSEA.R` |
| Fig. 2 | Immune infiltration: ssGSEA heatmap across GSE77459 samples | `11_immune_ssgsea.R` |
| Fig. 3 | GSE185222 single-cell atlas: UMAP and cell-type composition | `20_prep_GSE185222.R`, `21_annotate_GSE185222.R`, `22_composition_genes_GSE185222.R` |
| Fig. 4 | Three-stage compartment shifts and endothelial-cluster composition | `36_EMTAB_validate.R`, `54_composition_by_sample.R`, `40_EMTAB_endo_immune.R`, `67_compartment_effect_sizes_patient_figure.R` |
| Fig. 5 | Cell-type inflammatory programmes and the ACKR1+ venular interface | `39_EMTAB_celltype_DE_GSEA.R`, `40_EMTAB_endo_immune.R`, `27_keygene_localization.R` |
| Fig. 6 | Convergent regulatory inference: TF activity and PROGENy pathways | `44_TF_activity.R`, `47_PROGENy.R` |
| Fig. 7 | Bulk M1 elevation as a compositional signal, not a global M1/M2 switch | `58_myeloid_refine_M1M2.R`, `65_GSE92681_probe_annotation_M1M2.R`, `61_build_current_figures.py` |
| Fig. 8 | Descriptive communication and pseudotime analyses in E-MTAB-14896 | `50_EMTAB_comm_focus.R`, `51_pseudotime.R`, `57_LR_patient_perm.R`, `60_build_sr_figures.py` |

Figure assembly is at 300 dpi. The drug-signature panel is a supplementary direction-aware positive control and is not assigned to Fig. 7.

### Supplementary Tables 14–17 (the four former main tables)

| Item | Content | Principal scripts / inputs |
|---|---|---|
| Supplementary Table 14 | Public datasets, platforms, sample sizes and roles | Compiled from GEO/ArrayExpress metadata; PBM details from `12_PBM_crosscheck.R` and `15_PBM_deepdive.R` |
| Supplementary Table 15 | Cell-type composition of GSE185222 | `22_composition_genes_GSE185222.R` |
| Supplementary Table 16 | E-MTAB-14896 compartments and signature genes | `54_composition_by_sample.R`, `40_EMTAB_endo_immune.R` |
| Supplementary Table 17 | Selected transcription-factor activities (DoRothEA A+B, decoupleR ULM) | `44_TF_activity.R`, `68_recheck_TF_stage_stats.R` |

## 5. Statistical and reporting notes

- Patient-level values are used for three-stage composition, TF/pathway activity and communication inference wherever applicable. Cell-level P values are descriptive because cells are nested within patients.
- The total-endothelial advanced-minus-healthy effect is +27.235 percentage points (95% CI −5.120 to 59.591; Welch P = 0.08234). Four patients per stage limit power.
- In the primary M1/M2 marker set, M1 and M2 increase together in both bulk cohorts, but M1-minus-M2 does not change significantly (GSE77459 P = 0.699; GSE92681 P = 0.755). With the alternative marker set, the corresponding P values are 0.041 and 0.106. The result is therefore marker-set dependent and is not described as a consistent global polarisation switch.
- E-MTAB myeloid M1-minus-M2 is not significant with either marker set (P = 0.122 and 0.885). GSE185222 also shows marker-set dependence; activation, interferon/antigen-presentation and M2-core programmes are reported alongside the alternative-set contrast.
- All 118 E-MTAB transcription factors were rechecked with stage-wise Kruskal–Wallis tests and Benjamini–Hochberg correction. STAT1 and MYC have nominal P = 0.04978707; the minimum BH q is 0.3894635 and no TF passes q < 0.05. The complete recalculation is `results/singlecell/TF_stage_KW_BH_recalculated.csv`.
- ROC estimates are apparent, in-sample estimates from the same 12 discovery samples used for candidate selection. Among genes displayed in Supplementary Fig. 6, TNF has AUC = 0.8888889 and the other displayed genes have AUC = 1.00. These values are not diagnostic performance estimates.
- The PBM analysis averages two Ct replicates per group and uses 72 non-missing paired LPS/PBM mean DeltaCt targets. It applies no replicate-concordance exclusion and no replicate-level differential-expression inference. Thirty-two of 72 targets reverse direction (exact binomial P = 0.41); this panel is descriptive only.
- Drug-signature enrichment is a pathway-recovery positive control, not evidence of drug efficacy.

## 6. Software and execution

The analyses were run locally with R 4.5.3 and Bioconductor 3.22. Figure assembly uses Python with Pillow; manuscript checks use python-docx and PyMuPDF. Principal installed package versions are listed in `package_versions.csv` inside the release ZIP.

Scripts retain the original Windows absolute paths for provenance. Before rerunning, edit the `BASE`, `RES`, `SC`, input and output constants to point to a local copy of the public data and this release. Do not densify the full E-MTAB matrix at once; process it patient by patient. Some scripts download annotation or MSigDB resources and therefore require network access.

`CODE_VALIDATION.md` distinguishes parse-only and figure-assembly checks from a complete end-to-end rerun. Raw large datasets are intentionally absent, so the release archive itself is not a one-click containerized rerun.

## 7. Release and citation

Repository: https://github.com/renzhongzhen/pulpitis-compartment-immunity  
Complete release archive: `pulpitis_compartment_immunity_code_20261120.zip`

Code is provided under the MIT License (`LICENSE`). Public datasets remain subject to their original submitters' terms. If using this release, cite the manuscript and the specific release archive.

```
Ren Z. Compartment-resolved immunity in irreversible pulpitis: analysis code and result tables.
GitHub repository. https://github.com/renzhongzhen/pulpitis-compartment-immunity
```
