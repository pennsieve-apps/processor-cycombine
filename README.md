# processor-cycombine

Batch correction of flow / mass cytometry data with **cyCombine** (Pedersen et al., Nat
Commun 2022), one of the two correction branches of the *Multi-batch Harmonization & QC*
workflow:

```
fcs-prep -> peacoqc -> { batch-assess(raw), cytonorm -> batch-assess(cytonorm), cycombine -> batch-assess(cycombine) } -> cytometry-report
```

**Why.** Samples stained and acquired on different days (batches) carry technical shifts
in marker intensity that can look like biology. cyCombine removes them without needing a
control sample in every batch (compare `processor-cytonorm`), by assuming that each cell
type should look the same in every batch apart from the biology you tell it to keep.

**How.** cyCombine pools cells from all samples and

1. normalises each batch on its own (`NORM_METHOD`, default z-score per batch) - only so
   that the next step groups cells by cell type instead of by batch;
2. clusters the cells with a self-organising map (kohonen, `XDIM` x `YDIM` = 64 nodes);
3. in every SOM node runs **ComBat** (Johnson et al. 2007, the empirical-Bayes method
   from gene-expression arrays) on the original marker values: for each marker it removes
   each batch's location and scale offset, borrowing strength across markers. Values are
   capped to the node's input range.

**Why the covariate matters.** ComBat explains each marker by batch plus the variables in
its model. If the condition of interest (here stimulated vs unstimulated) is not in the
model and conditions are spread unevenly over batches - one batch with more stimulated
samples - ComBat would read the stimulation effect inside that batch as a batch effect
and remove part of the biology. With `COVARIATE_COLUMN=condition` the condition is in the
model, so only the batch part is removed. cyCombine drops the covariate in a SOM node
where it is confounded with batch or where almost all cells share one condition. Set
`COVARIATE_COLUMN=none` only when conditions are perfectly balanced or absent.

**Subsampled output.** ComBat estimates and applies its correction in one call on the
cells it is given, and cyCombine has no train/apply split (its `batch_correct` /
`correct_data` return corrected values only for the input cells; sva's `ComBat` returns no
reusable model). So the correction cannot be learned on a subsample and applied file by
file. The app therefore corrects a seeded subsample of at most `MAX_CELLS_PER_SAMPLE`
events per sample (default 50,000, file order kept) and **writes only those events**: the
output files are subsampled, `$TOT` is updated and the `PENNSIEVE_NORM` keyword records
it, e.g. `cycombine;covariate=condition;som=8x8;norm=scale;subsample=50000/559685;seed=42`.
`MAX_CELLS_PER_SAMPLE=0` corrects every event, at a memory cost (see Runtime).

Inputs are already arcsinh-transformed by `fcs-prep`; cyCombine's own transform is not
used. Scatter, `Time` and other non-marker channels and all keywords are kept.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `fcs` | input (subfolder) | `<sample_id>.fcs` | Arcsinh-transformed FCS from `fcs-prep` / `peacoqc`. Looked for in `INPUT_DIR/fcs/`, else anywhere under `INPUT_DIR`. |
| `metadata` | input (subfolder) | `samples.csv`, `panel.csv` | `samples.csv`: `sample_id`, `BATCH_COLUMN`, `COVARIATE_COLUMN`. `panel.csv`: `channel`, `marker`, `type`. |
| `fcs` | output | `fcs/<sample_id>.fcs`, `fcs/provenance.json` | The corrected events (subsample, see above), all input channels and keywords, FCS 3.1 float32; marker channels corrected; `PENNSIEVE_NORM` added. |
| `qc` | output | `qc/cycombine/...` | cyCombine's evaluation, density plots, `provenance.json`, all under a `cycombine/` folder: in workflow B the report links the qc outputs of cytonorm and cycombine into one input port, where entries with the same name (`provenance.json`) would collide. |

### `qc/cycombine/` files

cyCombine's evaluation, as in its vignette: the corrected cells are clustered again with a
SOM (same grid and seed); both the uncorrected and corrected values are then compared
within those clusters. To bound memory the evaluation uses at most 20,000 of the
corrected cells per sample (seeded).

| File | Columns / content |
|---|---|
| `evaluation.csv` | `metric, value`: `emd_reduction` (cyCombine `evaluate_emd`: EMDs below 2 both before and after are dropped; reduction = sum(before - after) / sum(before); 1 = all batch distance removed), `mad_score` (`evaluate_mad`: median over cluster x batch x marker of the absolute change in median absolute deviation; small = within-batch spread, i.e. biology, preserved), counts. |
| `emd.csv` | `cluster, marker, batch1, batch2, emd_uncorrected, emd_corrected`: earth mover's distance between two batches' distributions of a marker within a cluster (`compute_emd`, bin size 0.1). |
| `emd_by_marker.csv` | `marker, emd_uncorrected, emd_corrected`: mean of the above per marker. |
| `mad.csv` | `cluster, batch, marker, mad_uncorrected, mad_corrected` (`compute_mad`). |
| `som_node_counts.csv` | `sample_id, batch, som_node, n`: cells per correction SOM node and sample. |
| `density_NN.png` | `plot_density`: per marker and batch, density before (uncorrected) and after (corrected), 9 markers per page, on at most 5,000 cells per sample. |
| `provenance.json` | Versions, parameters, inputs, whether the covariate was used, subsampling, evaluation, timings, per-sample events in/out. |

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `BATCH_COLUMN` | string | `batch` | `samples.csv` column with each sample's batch. |
| `COVARIATE_COLUMN` | string | `condition` | Column ComBat keeps in its model (see above); `none` = no covariate. |
| `XDIM`, `YDIM` | integer | `8`, `8` | SOM grid (cyCombine default); the correction is done per node. Smaller = coarser correction, for very different batches. |
| `NORM_METHOD` | `scale` \| `rank` \| `qnorm` | `scale` | Batch-wise normalisation used only for clustering. The authors recommend `scale` within one study, `rank` across studies / instruments, and advise against `qnorm`. |
| `MAX_CELLS_PER_SAMPLE` | string | `auto` | Events per sample corrected and written (seeded); 0 = all. `auto` keeps cells × markers within `MAX_VALUES`, at most 50000 per sample. |
| `MAX_VALUES` | integer | `30000000` | Memory budget of `auto`, in corrected values (cells × markers); the evaluation uses 40% of it. 30 M peaked at 5.3 GB. |
| `MARKERS` | string | (all `type=marker`) | Comma-separated markers or channels to cluster and correct. |
| `SEED` | integer | `42` | Seed for the subsample and both SOMs. |

ComBat is run parametric (cyCombine's default), method `ComBat`, no reference batch.

## Results on the demo data

Full demo data from `fcs-prep` (`data/work/prep`: 24 files, 10,212,112 events, 25 markers;
no PeacoQC yet), default parameters: 50,000 events per sample corrected (1,175,055 cells;
batch01_tcu2 has only 25,055), `COVARIATE_COLUMN=condition` (every batch has 2 stimulated
and 2 unstimulated samples, so the covariate is not confounded), SOM 8x8 (64 nodes used),
`NORM_METHOD=scale`. All 24 files corrected.

| Metric (cyCombine evaluation, 480,000 cells) | Value |
|---|---|
| EMD reduction | **0.756** (9,008 of 24,000 node x marker x batch-pair EMDs were >= 2 before or after) |
| MAD score | 0.022 (median change of within-batch spread; small = biology kept) |

Mean EMD between batches per marker, before -> after (`emd_by_marker.csv`): GZMB 4.89 ->
1.58, CD14 4.88 -> 0.85, CCR7 3.18 -> 0.85, CD45RA 2.86 -> 0.89, HLADR 2.79 -> 0.71,
CD3 2.10 -> 0.58, CD19 0.85 -> 0.37; every marker improves. The same with 25,000 events
per sample gives EMD reduction 0.753 and MAD score 0.022, so the subsample size barely
changes the correction. Caveat for interpretation: ComBat has no notion of donor, and
the stimulated samples are different donors in different batches (see Limitations), so
some of what is removed in them is donor biology. Full outputs:
`workflow-immune-cytometry/data/work/cycombine/`.

## Runtime

Memory grows with cells × markers, not cells alone, so the default cap is a budget:
`MAX_CELLS_PER_SAMPLE=auto` takes `min(50000, MAX_VALUES / (markers × samples))` per sample,
and the evaluation `min(20000, 0.4 × MAX_VALUES / (markers × samples))`. On the spectral demo
(24 samples × 25 markers) that is 50000 and 20000, as before. On the CyTOF demo (48 samples ×
42 markers, Kowli et al. 2025, 20k events per file) it is 14880 and 5952: 235 s, R-reported peak
3.0 GiB; with 50000 per sample the evaluation alone was killed at 6 GB.

`computeTypes: [standard]`, 2 vCPU / 8 GiB, 2 h timeout.

Measured on the full demo data (24 files, 10.2 M events, 1.3 GB), linux/amd64 under
emulation on an Apple-silicon Mac (Docker Desktop, `--cpus 4`), cgroup memory (includes
the page cache of the 1.3 GB of input read twice):

| `MAX_CELLS_PER_SAMPLE` | cells corrected | wall time | peak memory |
|---|---|---|---|
| 25,000 | 600,000 | 100 s | 4.3 GB |
| **50,000 (default)** | 1,175,055 | **128 s** (read 9 s, SOM 32 s, ComBat 17 s, write 7 s, evaluation 43 s, plots 13 s) | **5.3 GB** |
| 100,000 | 2,318,290 | - | killed at 3.9 GB by the shared 8 GB Docker VM (other jobs running), not measurable here |
| 0 (all events) | 10,212,112 | not run | extrapolated ~20 GB |

Memory grows by about 1.7 GB per million corrected cells (slope between the two
measured points; the pooled table, its batch-normalised copy, the per-node splits that
ComBat works on and the corrected copy all exist at once), on top of ~3.3 GB that is
mostly page cache. So the default fits 8 GiB; correcting every event of the demo data
(10.2 M cells) would need roughly 20-24 GB (`MAX_CELLS_PER_SAMPLE=0` with a 30 GiB
task), which is why the default subsamples. The evaluation is capped at 20,000 cells per
sample and does not grow further.

The image is 579 MB (R 4.6.1, flowCore, sva, kohonen, cyCombine and dependencies).

## Limitations

- The output files are subsampled (default 50,000 events per sample): downstream steps
  see fewer events than the input had. Rare populations (< 0.1 %) become thin.
- Correction assumes every SOM node is the same cell type in every batch and that, after
  the covariate, the remaining differences between batches are technical. On the demo
  data the stimulated samples come from donor 199C in batches 01 and 05 and donor 198C in
  02, 03, 04 and 06 (FCS `$TUBENAME`), so part of the donor difference will be removed as
  "batch".
- ComBat removes per-batch mean and variance differences within a node; non-linear
  distortions (a compressed top end) are only partly corrected (CytoNorm's splines handle
  those).
- Memory grows linearly with the number of corrected cells (all held at once).

## Running it locally

```sh
pip install -r requirements-dev.txt
make fixture DEMO=../workflow-immune-cytometry/data/zenodo-15388817   # 20,000 events per file
make test                                                             # unit tests (synthetic data)
make run                                                              # docker compose over data/fixture
```

Or by hand:

```sh
docker build --platform linux/amd64 -t processor-cycombine .
docker run --rm --platform linux/amd64 -e MAX_CELLS_PER_SAMPLE=50000 \
  -v "$PWD/data/fixture:/data/input:ro" -v "$PWD/data/output:/data/output" \
  processor-cycombine
```

`make test` runs `tests/run_tests.R` with a local R that has the packages, else in the
image. The tests build a synthetic 3-batch fixture and check channels, keywords and
non-marker values, that replicate medians align, the subsample (count, same events,
keyword, reproducible with `SEED`), no-covariate mode, unlinked inputs and per-sample
failures.

## Licences

This app's code is Apache-2.0 (`LICENSE`). The image bundles R and R packages under their
own licences (cyCombine MIT; sva Artistic-2.0; R and kohonen GPL; flowCore's backend
cytolib AGPL-3.0), so the image as a whole is distributed under GPL/AGPL terms; see
`THIRD_PARTY.md`.

## Citations

- Pedersen CB, Dam SH, Barnkob MB, Leipold MD, Purroy N, Rassenti LZ, Kipps TJ, Nguyen J,
  Lederer JA, Gohil SH, Wu CJ, Olsen LR. cyCombine allows for robust integration of
  single-cell cytometry datasets within and across technologies. *Nat Commun*
  13(1):1698 (2022). doi:10.1038/s41467-022-29383-5. PMID 35361793.
- Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data
  using empirical Bayes methods. *Biostatistics* 8(1):118-127 (2007).
  doi:10.1093/biostatistics/kxj037. PMID 16632515. (ComBat)
- Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. The sva package for removing batch
  effects and other unwanted variation in high-throughput experiments. *Bioinformatics*
  28(6):882-883 (2012). doi:10.1093/bioinformatics/bts034. PMID 22257669.
- Wehrens R, Kruisselbrink J. Flexible Self-Organizing Maps in kohonen 3.0. *Journal of
  Statistical Software* 87(7) (2018). doi:10.18637/jss.v087.i07.
- Hahne F, LeMeur N, Brinkman RR, Ellis B, Haaland P, Sarkar D, Spidlen J, Strain E,
  Gentleman R. flowCore: a Bioconductor package for high throughput flow cytometry.
  *BMC Bioinformatics* 10:106 (2009). doi:10.1186/1471-2105-10-106. PMID 19358741.
- Demo data: Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow
  cytometry dataset. Zenodo (2025). doi:10.5281/zenodo.15388817 (CC BY 4.0).

Software: cyCombine https://github.com/biosurf/cyCombine at commit `fcd42dbb` (version 0.3.0).
