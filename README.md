# CommunityCarbon — Northern Peatlands SOC Workflow

Mapping soil organic carbon near **Fort Severn, Ontario** (Hudson Bay Lowlands)
from eight community-collected soil cores, supported by open spatial data and
published reference layers.

The short version of what this workflow found: **eight cores in two clusters
cannot support a covariate-driven carbon map of this landscape.** They can
support a careful accounting of what the cores measure and a stratified
estimate with honest uncertainty. Both are produced. A prediction map is
attempted, tested against a null model, and released only if it passes.

Full findings: **[`outputs/REPORT.qmd`](outputs/REPORT.qmd)** (technical) and
**[`outputs/COMMUNITY_REPORT.qmd`](outputs/COMMUNITY_REPORT.qmd)** (plain
language). Both are Quarto documents with **no executable code chunks** —
every number is computed by the R scripts and written in as finished text — so
they render with the Quarto CLI alone:

```bash
quarto render outputs/COMMUNITY_REPORT.qmd    # no knitr, no R packages needed
```
Community-facing summary and shareable figures: **[`outputs/COMMUNITY_BRIEF.md`](outputs/COMMUNITY_BRIEF.md)** and `outputs/figures/community_*.png`.

---

## Quick start

```bash
Rscript run_all.R               # 01,02,05,09,13,14 → 07,08,10,11,12 — base R
Rscript run_all.R --gee         # adds 03, 04, 06     — needs rgee + EE auth
Rscript tests/test_functions.R  # 193 unit tests
Rscript tests/test_pipeline.R   # 131 integration tests
```

**Scripts 01, 02, 05 and 07–12 require no R packages at all** — not even
for spatial output. That is deliberate: the scientifically load-bearing steps —
quality control, stock computation, stratified estimation, the validation
ledger, and the GeoTIFFs — should run on any R installation without a
package-install ordeal, because a community-science pipeline that partners
cannot re-run is not reproducible in any useful sense. `R/geotiff.R` writes
valid GeoTIFFs and GeoJSON directly, so GDAL and `terra` are never required to
produce a map. Only the Earth Engine steps (03, 04) and the model that consumes
them (06) need `rgee`, and the pipeline completes and reports without them.

## Spatial outputs

`09_spatial_products.R` writes to `outputs/spatial/` and `outputs/figures/`:

| File | Status |
|------|--------|
| `cores.geojson` | OBSERVATION — the 8 cores with stocks, flags, depth support |
| `dist_to_nearest_core_m.tif` | GEOMETRY — exact, no inference |
| `nearest_core_index.tif` | DIAGNOSTIC — Thiessen assignment |
| `DIAGNOSTIC_soc_nearest_core_kgm2.tif` | DIAGNOSTIC — **not** a carbon map |
| `PRODUCT1b_soc_within_credible_radius_kgm2.tif` | SAMPLE-BASED, MASKED |
| `09_*.png` | quick-look figures |

The masked product is the only core-derived surface this workflow will call a
product. Taking the credible radius as the median core-to-core spacing
(2.5 km), **the eight cores constrain 8.4% of the 1,417 km² study area** — 119
km². The rest is NoData because nothing in the dataset constrains it. That
number is the argument for the next field season rather than for a better
interpolator.

## Bayesian maps (PRODUCT 4)

`11_bayesian_map.R` treats a published carbon map as the **prior** and the
cores as **observations**, in a conjugate Gaussian update (`R/bayes.R`). A
core's pull on a pixel fades with distance (Gaussian kernel, length scale =
median core spacing, truncated at 7.5 km) and precisions add across cores, so
sample size enters automatically. Outputs: posterior mean, sd,
`core_info_fraction` (the share of the answer coming from cores rather than the
prior — the most honest layer in the set) and `shift_from_prior`.

**Depth is handled explicitly, because Li et al. is full peat-column carbon
(~86 kg/m² over ~184 cm) and the cores are 0–30 cm (3–14 kg/m²).** Those differ
~6× because they measure different things. So the column is split into its top
30 cm and everything below, only the top is updated, and the deep part is added
back untouched. The split fraction is the single largest assumption and is
configurable, reported, and varied in sensitivity.

**The headline result is not the map update — it is a corroboration.** The
cores move the full-column prior by at most 0.21%, which is correct: 30 cm is a
small slice of a metres-deep column. But a measured 0–30 cm stock of
10.02 kg C/m² against Li's 86 kg C/m² implies **11.6%** of the column sits in
the top 30 cm, against **16.3%** under uniform density. Peat compacts with
depth, so *less* than proportional is exactly what physics predicts — and these
cores show bulk density rising 4.6× from surface to 30 cm. **The community
measurement independently corroborates the published regional map at this
site.**

## Earth Engine outputs go to Drive *and* local disk

`R/gee_io.R` exports every Earth Engine raster to Google Drive and downloads it
to `outputs/gee/`. The Drive copy is for sharing; the local copy is what the
reports embed, what `terra` can read, and — the part that changes the analysis
— what `11_bayesian_map.R` looks for.

**Running `04` upgrades the Bayesian maps automatically.** With no reference
raster on disk, `11` falls back to a spatially constant regional prior and
prefixes its outputs `DEMO_`. The downloads land under exactly the filenames
`11` expects, so re-running it afterwards produces the real Li-prior product
with no further action.

A failed download never loses work: the Drive export is started first and
succeeds on its own, downloading is a separate recoverable step, and each
product's fate is recorded in an export manifest (`drive_and_local`,
`drive_only`, `drive_pending`, `cached`). Downloaded files are verified with
`terra` — dimensions, CRS, and whether the raster is entirely NoData, which a
successful download can still leave you with.

## Coast-following AOI (script 13)

An axis-aligned box fits a coastal study area badly — it either clips the
sampled transect or spends most of its area over open water. `13` builds a
rectangle **aligned to the trend of the shoreline**, and derives that trend
from the cores' own principal axis rather than drawing it by eye: **142.9°
from north**, 8° off the NW–SE Hudson Bay coastal trend, with 80% of the core
variance on that axis. With a 30 km along-axis and 12 km across-axis buffer the
AOI is **83.6 × 33.0 km = 2,758 km²**, which is 42% of its own axis-aligned
envelope — so the rotation genuinely earns its keep.

Outputs: `aoi_coast_oriented.geojson`, `aoi_coast_oriented_mask.tif`, and
`aoi_for_earthengine.js` (a `ee.Geometry.Polygon` literal, so the AOI used
locally and in GEE are provably the same shape). Override the bearing via
`CFG$aoi_oriented$bearing_deg`; the script falls back to an axis-aligned box if
the cores are too clustered for a direction to mean anything.

The coast-relative coordinates it produces expose something the lat/lon view
hides: **all peat cores sit 3.6–5.1 km to one side of the transect axis and all
mineral cores 1.0–3.9 km to the other.** Peat-versus-mineral is also
one-side-versus-the-other, which compounds the campaign confounding.

## National context (script 14)

The AAFC National Pedon Database (9,017 pedons, 349 complete to 30 cm after QC)
places Fort Severn nationally. Two results:

- **The nearest existing core with a usable 0–30 cm stock is 700 km away, and
  there are zero within 500 km.** The community cores are not merely few —
  they are *all there is*. Every published carbon map covering Fort Severn is,
  in this neighbourhood, extrapolating from hundreds of kilometres away.
- **The Fort Severn peat/mineral ordering appears nationally, but weakly.**
  Organic-horizon soils median 6.23 kg C/m² over 0–30 cm vs 7.15 for mineral
  (Wilcoxon p = 0.031; a random organic core exceeds a random mineral one 43%
  of the time). The ordering matches, but IQRs overlap heavily and organic
  soils are simply far more variable (IQR 9.0 vs 4.5). Read as corroboration,
  not proof. What survives: **a 0–30 cm accounting does not rank peatlands
  above mineral soils** — depth is what distinguishes them.

The community cores sit at the 24th–89th national percentile, i.e. squarely
within the national range, which is itself a check on data quality.

### Peat depth: a substantive field finding

Detecting the peat/mineral contact from the organic-matter profile gives peat
depths of **15.2 cm (PM-02-A)** and **25.2 cm (PM-03-A)**, with PM-01-A still
in peat at 14.5 cm. Against a regional mean of 184 ± 48 cm, these are **peat
margins or shallow fen, not the deep peat plateaus** that hold most of the
region's carbon — which is both a limit on what these cores can say about a
full-column product and a useful conservation finding in its own right, since
thin peat is the vulnerable kind.

---

## Layout

```
config.R                      single shared config: paths, CRS, seed, thresholds,
                              depth windows, GEE asset ids. Sourcing it also
                              loads R/, so scripts need one bootstrap block.
R/
  utils.R                     logging, I/O, small pure helpers
  qc.R                        hard gates and value flags (pure predicates)
  depth.R                     depth semantics and window integration (pure)
  stats.R                     bootstrap, cross-validation, variogram diagnostic (pure)
  gee_io.R                    Earth Engine export to Drive AND local disk
  bayes.R                     conjugate Bayesian update, kernels, depth splitting (pure)
  geotiff.R                   base-R GeoTIFF and GeoJSON writers (no GDAL)
scripts/
  01_ingest_qc.R              read, recompute from first principles, reconcile, flag
  02_depth_harmonise.R        like-for-like depth windows, stocks, lower-bound flags
  03_gee_covariates.R         rgee: buffered covariate stack + extraction [needs GEE]
  04_gee_reference_audit.R    rgee: audit units/depth of published layers [needs GEE]
  05_stratified_estimate.R    PRODUCT 1 and 2; the case against kriging
  06_covariate_model.R        PRODUCT 3 + honest cross-validation        [needs GEE]
  07_validation_ledger.R      what is validated, what is not, where leakage was possible
  08_report.R                 assembles outputs/REPORT.qmd (technical)
  09_spatial_products.R       GeoTIFFs + GeoJSON, no GDAL required
  10_community_figures.R      shareable PNG figures + COMMUNITY_BRIEF.md
  11_bayesian_map.R           PRODUCT 4: published map as prior, cores as data
  12_community_report.R       plain-language outputs/COMMUNITY_REPORT.qmd
  13_aoi_boundary.R           coast-following oriented AOI (GeoJSON + mask + GEE)
  14_npdb_context.R           national context from the AAFC Pedon Database
  15_landcover_carbon.R       carbon per ecosystem type                 [needs GEE]
  16_embedding_clusters.R     k-means on AlphaEarth embeddings          [needs GEE]
tests/                        unit + integration tests
data/raw/                     the 22-row core CSV, verbatim headers
data/derived/                 script outputs
outputs/tables|figures|gee/   tables, and GEE rasters when retrieved
```

---

### Earth Engine covariate extraction note

`03_gee_covariates.R` expands the mapping AOI around the actual core locations
with a configurable 50 km buffer before compositing and sampling. It also
unmasks the sampling stack with a no-data sentinel so Earth Engine cannot drop
an entire core just because one optical/radar/reference band is masked; the
sentinel is converted back to `NA` in the CSV, and a per-band coverage table is
written to `outputs/gee/03_covariate_coverage.csv`.

## Products, and what each one is

| # | Product | Scientific status |
|---|---------|-------------------|
| 1 | Stratum mean stock | **Sample-based estimate** (design-based) |
| 2 | Stratum-mean assignment map | **Class-mean assignment** — *not* interpolation, *not* prediction |
| 3 | Covariate prediction | **Prediction / extrapolation** — released only if it beats a null |
| R | SoilGrids / Sothe / Li | **Reference** (external products, audited before use) |

No kriged or IDW surface is produced. Eight cores give 28 pairs; zero lag bins
reach the ~30 pairs conventionally needed before a bin is informative, and the
two campaigns form two spatial clusters that differ in landscape type, so any
apparent spatial correlation would largely be the peat/mineral contrast
reappearing as a distance effect. `05` demonstrates this rather than asserting
it.

---

## Four things quality control found

All four would have been invisible had the supplied columns been trusted.

1. **Longitudes were in the wrong hemisphere.** All 22 rows carried positive
   longitudes; Fort Severn is at −87.63. As supplied the cores plot ~7,540 km
   away in western Siberia. This is a hard gate: the run stops unless the
   repair is explicitly authorised, and every repaired row is flagged.

2. **`Organic Carbon Density (g/cm^2)` is off by one dimension.** It reconciles
   exactly as `bulk density × SOC/100`, which is g/cm**³** — a concentration,
   not an areal density. A prior pipeline summed it down each core; for the
   peat cores, where segment thickness varies, that error ranges from 7.7× to
   14.5× and cannot be undone by rescaling.

3. **The two campaigns used different OM→carbon factors** — 0.580 (van
   Bemmelen) for the peat cores, 0.554 (≈1/1.8) for the mineral ones. `SOC` is
   two differently-derived quantities sharing a column, and the peat cores
   carry the factor that overstates peat carbon. Campaign, landscape, year,
   operator and conversion factor are perfectly confounded.

4. **`Depth` is segment thickness, not cumulative depth** — derived from the
   data (it decreases within two cores) rather than assumed.

---

## Depth handling

Two windows, never interchangeable:

- **Common support, 0–14.5 cm** — the deepest window every core fully covers,
  derived from the data. No extrapolation; all eight cores directly comparable.
  This is the response variable for any model.
- **Reference, 0–30 cm** — matches the published layers. Six cores reach it;
  PM-01-A (14.5 cm) and PM-02-A (20.7 cm) do not, so their stocks are **lower
  bounds**, flagged and never rescaled to look complete.

**A 0–30 cm map of this landscape is not a map of its carbon.** Peat here has
bulk density of 0.05–0.27 g/cm³ against 0.75–2.18 for the mineral soils, so
over the common-support window the peat cores actually hold *less* carbon
(3.08 kg/m²) than the mineral cores (5.26 kg/m²). The peatland advantage lives
in depth: Hudson Bay Lowlands peat is metres thick and nearly all its carbon
sits below the window every 0–30 cm product reports.

---

## Cross-validation: leave-one-core-out, never leave-one-segment-out

The most consequential choice in the workflow. There are 22 segments, which
looks like a workable sample size — but every segment of a core sits at one
coordinate pair, so at 30 m resolution they all sample the *same pixel* of
every covariate, and they share a campaign, operator and lab batch. The
effective sample size for a covariate model is **8**, not 22. Segment-level
splitting would produce impressive statistics measuring the model's ability to
recall cores it had already seen. No script here does it; `06` demonstrates the
size of the inflation on a controlled example.

Leave-one-**stratum**-out is also run — train on mineral, predict peat — since
mapping means predicting ground unlike the training points. Every model is
scored against a mean-only null, and `06` writes no raster when the model
cannot beat it.

Validation against SoilGrids and Sothe et al. is **quasi-external**: both
predate these cores, but both were trained on profile databases that plausibly
include Hudson Bay Lowlands profiles from the same population. They have seen
the population without seeing the cores.

---

## Reference layers (`Prior_data` catalogue)

Asset IDs, units and depth bases come from the project's `Prior_data` GEE
reference library, ingested into `config.R`. `04_gee_reference_audit.R` still
verifies each one at runtime rather than trusting the catalogue.

**Li et al. 2025** — `McMasterCarbon30mkgm2version1`: the `30` is the **pixel
size**, not a depth. The product is **full peat-column** carbon, surface to the
base of the peat, with a published HBL mean of **86 ± 35 kg C/m² over a mean
peat depth of 184 ± 48 cm**. Comparing it to a 0–30 cm stock is a category
error — 30 cm is 16% of the mean peat depth. `04` confirms the pixel grid is
~30 m, checks the AOI mean against the published mean in standard-deviation
units, and refuses to combine layers whose depth support differs. (Earlier
project notes labelled this asset "0–100 cm", which disagrees with the paper;
the runtime check is what settles it.)

**Sothe et al. 2022** ships at **two depths** — 0–30 cm and 0–1 m — with
identical units. `04` audits both and bars the 0–1 m layer from the 0–30 cm
comparison, checking that the deeper layer exceeds the shallower one as nesting
requires.

**SoilGrids 2.0** is built here from `soc_mean` × `bdod_mean` rather than the
packaged `ocs_mean`, so the integration matches the one applied to the cores:
`(soc/10 g/kg) × (bdod/100 g/cm³) × thickness_cm / 100 = kg C/m²`.

**GWL_FCS30** (Zhang et al. 2023) is the preferred stratum layer over ESA
WorldCover, which reads treed bog on peat and upland forest on mineral soil
alike — precisely the distinction these two campaigns represent.

---

## Community figure pack

The simplified communication workflow now ends with `10_community_figures.R`,
which produces a plain-language brief and four PNGs for meetings, slides and
printouts:

1. `community_01_core_locations.png` — where the peat/wetland and
   mineral/forest cores are located.
2. `community_02_core_carbon_stocks.png` — what the new cores measured in the
   shared 0-14.5 cm layer.
3. `community_03_compare_context.png` — why shallow core values should not be
   compared directly with full peat-column stocks.
4. `community_04_conservation_messages.png` — the conservation takeaways.

These figures answer four community questions directly: carbon stock around
Fort Severn; what the new peat/mineral samples taught us; how the values
compare with shallow products and full peat-column estimates; and why the
conservation priority is keeping peat wet and measuring deeper peat stores.

## Reproducibility

Seed `20260727`, set in `config.R` and used by every stochastic step. QC flags,
never drops — implausible values are surfaced with named flags and the row is
kept, because dropping data is the analyst's decision, not the pipeline's.
`outputs/REPORT.md` is generated entirely from the derived tables.
