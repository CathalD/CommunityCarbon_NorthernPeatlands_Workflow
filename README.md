# CommunityCarbon — Northern Peatlands SOC Workflow

Mapping soil organic carbon near **Fort Severn, Ontario** (Hudson Bay Lowlands)
from eight community-collected soil cores, supported by open spatial data and
published reference layers.

The short version of what this workflow found: **eight cores in two clusters
cannot support a covariate-driven carbon map of this landscape.** They can
support a careful accounting of what the cores measure and a stratified
estimate with honest uncertainty. Both are produced. A prediction map is
attempted, tested against a null model, and released only if it passes.

Full findings: **[`outputs/REPORT.md`](outputs/REPORT.md)**.

---

## Quick start

```bash
Rscript run_all.R               # 01, 02, 05, 09, 07, 08 — base R only, ~5 s
Rscript run_all.R --gee         # adds 03, 04, 06        — needs rgee + EE auth
Rscript tests/test_functions.R  # 100 unit tests
Rscript tests/test_pipeline.R   # 63 integration tests
```

**Scripts 01, 02, 05, 07, 08 and 09 require no R packages at all** — not even
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
scripts/
  01_ingest_qc.R              read, recompute from first principles, reconcile, flag
  02_depth_harmonise.R        like-for-like depth windows, stocks, lower-bound flags
  03_gee_covariates.R         rgee: covariate stack + extraction         [needs GEE]
  04_gee_reference_audit.R    rgee: audit units/depth of published layers [needs GEE]
  05_stratified_estimate.R    PRODUCT 1 and 2; the case against kriging
  06_covariate_model.R        PRODUCT 3 + honest cross-validation        [needs GEE]
  07_validation_ledger.R      what is validated, what is not, where leakage was possible
  08_report.R                 assembles outputs/REPORT.md
tests/                        unit + integration tests
data/raw/                     the 22-row core CSV, verbatim headers
data/derived/                 script outputs
outputs/tables|figures|gee/   tables, and GEE rasters when retrieved
```

---

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

## Reproducibility

Seed `20260727`, set in `config.R` and used by every stochastic step. QC flags,
never drops — implausible values are surfaced with named flags and the row is
kept, because dropping data is the analyst's decision, not the pipeline's.
`outputs/REPORT.md` is generated entirely from the derived tables.
