# Sample inventory — *A diverse and distinct microbiome inside living trees*

This README explains the sample naming conventions and how to find any sample across the project (NCBI SRA, GitHub repo, Figshare archive, and the published analysis).

**If you just want to know what samples exist and how to retrieve them, start with `master_sample_inventory.csv`. Everything else is supporting detail.**

Related resources:
- **NCBI BioProject:** [PRJNA1124946](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1124946) (SRA Study SRP514484)
- **GitHub repo:** [jgewirtzman/tree-microbiome](https://github.com/jgewirtzman/tree-microbiome)
- **Figshare:** [doi.org/...](https://figshare.com/s/b781d3c5c24f1089ed5d)
- **Paper:** Arnold & Gewirtzman et al., *Nature*, 2025 — [doi:10.1038/s41586-025-09316-0](https://doi.org/10.1038/s41586-025-09316-0)

---

## The master sample inventory

`master_sample_inventory.csv` has one row per (sample × marker) pair, joining the SRA submission, the DADA2 input manifests, and the analysis-ready metadata tables. Columns:

| Column | Description |
|---|---|
| `sample_uid` | Canonical identifier, e.g. `WP7Inner_16S_S207`. The same string used in FASTQ filenames. |
| `prefix` | Sample stem without the marker/well suffix, e.g. `WP7Inner`. |
| `marker` | `16S` (bacteria/archaea) or `ITS` (fungi). |
| `tree_id` | The biological individual the sample comes from, e.g. `WP7`. |
| `compartment` | What part of the tree/soil this is, e.g. `Inner`, `Outer`, `Mineral`, `Organic`, `Bark`, `Heart1`, etc. |
| `material` | Higher-level grouping: `Wood`, `Soil`, `Bark`, `Litter`, `Root`, `Wetland/Marginal`, `Mock`, etc. |
| `species_code` | 4-letter species code (see table below). |
| `species_latin` | Full Latin binomial. |
| `common_name` | English common name. |
| `archive_source` | `Forest` for the main forest plot study; `BlackOak` for the QUVE deep-sampling study. |
| `plate_well_S` | The Illumina sample-plate well number. |
| `sra_run_SRR` | NCBI SRA Run accession (the actual data file). |
| `sra_biosample_SAMN` | NCBI BioSample. |
| `sra_sample_SRS` | NCBI SRA Sample. |
| `sra_experiment_SRX` | NCBI SRA Experiment. |
| `sra_spots_reads` | Read count (paired-end spots). |
| `sra_release_date` | When the record was released. |
| `sra_layout` | `PAIRED` or (rarely) `SINGLE`. |
| `sra_status` | `on_SRA` or `missing_from_SRA`. |
| `analysis_status` | `used_in_analysis`, `excluded`, or `not_in_analysis_tables`. |

### Status fields

- **`sra_status`**
  - `on_SRA` — the sample has SRA accessions and can be downloaded from NCBI.
  - `missing_from_SRA` — the sample is in the authors' pipeline but was not deposited to SRA.

- **`analysis_status`**
  - `used_in_analysis` — appears in the published analysis tables with a valid species/tree assignment.
  - `excluded` — present in the metadata but not used in the published analysis (wetland-margin samples, ambiguous-tree samples, subsampling-test replicates, etc.).
  - `not_in_analysis_tables` — appears nowhere in the analysis tables (mocks and MISC controls).

---

## Species code → tree prefix mapping

| Code | Latin name | Common name | Tree-ID prefix(es) | # Trees | # Samples (16S+ITS) | # On SRA |
|---|---|---|---|---|---|---|
| `ACRU` | *Acer rubrum* | Red maple | `RM, RMS` | 17 | 126 | 112 |
| `ACSA` | *Acer saccharum* | Sugar maple | `SM` | 15 | 114 | 114 |
| `BEAL` | *Betula alleghaniensis* | Yellow birch | `CUL, YB, YBYB` | 14 | 108 | 98 |
| `BELE` | *Betula lenta* | Black/sweet birch | `BB, BBGT, BBX` | 19 | 150 | 150 |
| `BEPA` | *Betula papyrifera* | Paper birch | `PB` | 6 | 48 | 48 |
| `CAOV` | *Carya ovata* | Shagbark hickory | `SH` | 1 | 8 | 8 |
| `FAGR` | *Fagus grandifolia* | American beech | `AB, ABT` | 14 | 108 | 108 |
| `FRAM` | *Fraxinus americana* | White ash | `WA, WAX` | 12 | 82 | 82 |
| `KALA` | *Kalmia latifolia* | Mountain laurel | `ML, MLGP` | 2 | 16 | 16 |
| `PIST` | *Pinus strobus* | Eastern white pine | `WP, WPBERENICE, WPX` | 20 | 136 | 122 |
| `PRSE` | *Prunus serotina* | Black cherry | `BC` | 3 | 22 | 22 |
| `QUAL` | *Quercus alba* | White oak | `WO` | 1 | 8 | 8 |
| `QURU` | *Quercus rubra* | Northern red oak | `RO` | 16 | 108 | 96 |
| `QUVE` | *Quercus velutina* | Black oak | `BO, QUVE` | 16 | 128 | 128 |
| `SAAL` | *Sassafras albidum* | Sassafras | `SA` | 1 | 8 | 8 |
| `TSCA` | *Tsuga canadensis* | Eastern hemlock | `H, HX, HZ` | 16 | 118 | 104 |

---

## Sample-name anatomy

A typical sample identifier looks like:

```
WP7Inner_16S_S207
└┬┘ └─┬─┘ └┬┘ └┬─┘
 │    │    │   └── plate well number
 │    │    └────── marker gene (16S = bacteria, ITS = fungi)
 │    └─────────── compartment / material location
 └──────────────── tree identifier
```

The original FASTQ filenames append `_R1_001.fastq.gz` for the forward read.

### Compartments (what part of the tree/soil)

| Compartment | Material | Notes |
|---|---|---|
| `Inner` | Wood | Inner 5 cm of core (operationally defined heartwood) |
| `Outer` | Wood | Outer 5 cm of core (operationally defined sapwood) |
| `Heart`, `Heart1`, `Heart2` | Wood | Whole heartwood |
| `Sap` | Wood | Whole sapwood |
| `Mineral` | Soil | Mineral soil horizon |
| `Organic` | Soil | Organic soil horizon |
| `Bark`, `Bark1`, `Bark2` | Bark | (mostly QUVE Black Oak study) |
| `Litter1`–`Litter4` | Litter | Leaf litter |
| `Branch1`–`Branch4` | Branch | Branch wood |
| `Foliage1`–`Foliage4` | Foliage | Leaves |
| `Coarse*`, `Fine*` | Root | Coarse / fine roots |
| `Muck`, `TopMuck`, `BottomMuck`, `Aquatic`, `MicrobialGoo` | Wetland/Marginal | Wetland-edge and standing-water samples (excluded from analysis) |
| `InsideRot` | Wood | Rotted-interior wood |

---

## File structure (where data lives)

| Source | What's there | Best for |
|---|---|---|
| **NCBI SRA — PRJNA1124946** | Raw FASTQ files (1,233 runs) | Re-running analyses from scratch |
| **GitHub: `Other-Metadata/`** | Plate maps, ddPCR, tree core sectioning data | Lab-bench-level metadata |
| **GitHub: `Other-Metadata/annotated_metadata/`** | Final analysis-ready meta tables | What was actually used in the paper |
| **Figshare: `Forest_Microbiome.zip`** | DADA2 outputs, OTU tables, taxonomy, trees, mapping files for the main forest study | Skipping the FASTQ→OTU step |
| **Figshare: `Black_Oak.zip`** | Same, but for the Black Oak (QUVE) deep-sampling study | Black Oak analyses only |
| **Figshare: GC, ddPCR, metabolisms data** | Gas chromatography, qPCR, FAPROTAX outputs | Functional analyses |

---

## Notes on the data

A few things are worth knowing before you work with the deposit. None of them affect the underlying read data; they're metadata details.

**Some samples on SRA were not used in the published analysis.** 29 sample-prefixes (50 runs) are present on SRA but were excluded from the analysis tables — wetland-margin samples (`*Muck`, `*Aquatic`, `*MicrobialGoo`), subsampling-test replicates, and a few others. They are flagged `analysis_status = excluded` in the master sheet. Two MISC controls and two mock community samples (`Zymo_mock`, `Fungal_mock`) are similarly flagged `not_in_analysis_tables`.

**One sample appears twice on SRA.** `WA5Inner_ITS_S1001` was uploaded in both the 2024 and 2026 submissions (`SRR29444616` and `SRR37072217`), with conflicting layout metadata in the two records. Same underlying data — pick one and ignore the other.

**The SRA `LibraryName` field is unreliable.** 1,232 of the 1,233 records have `LibraryName` containing literal Excel formula text (`= LEFT(A5, ...)` etc.) instead of a library name. Use `SampleName` (which is correct) or join via `Run` / `BioSample` accessions.

---

## How to cite

If you use these data, please cite the paper:

> Arnold, A.E., Gewirtzman, J., et al. (2025). A diverse and distinct microbiome inside living trees. *Nature*. https://doi.org/10.1038/s41586-025-09316-0

And the SRA deposit: NCBI BioProject PRJNA1124946.
