# CommunityCarbon MVP

A simplified, 7-step version of the full workflow one level up in this repo
(`../scripts/`, `../R/`, untouched). Same scientific backbone -- cores →
carbon stocks → satellite-covariate model → landscape prediction → Bayesian
fusion with a published prior -- built from mature R packages instead of
~10,800 lines of custom research-grade code. See
`.claude/brainstorm-mvp.md` in the repo root for the full design rationale
if you want it; this README is just "how do I run it."

**This pipeline is iterative.** The first run fuses a published carbon map
(Li et al. 2025) with your 8 cores. Every run after that fuses your last
run's own output with whatever new cores you've added since -- the map gets
better each field season without ever being rebuilt from zero. See step 6.

## Setup (once)

```r
source("mvp/install_packages.R")
```

Then, for step 02 only, set up Earth Engine once interactively:

```r
rgee::ee_install()
rgee::ee_Initialize(project = "your-gcp-project-id")
```

(The project ID defaults to the one in `config.R`; override with the
`EE_PROJECT` environment variable if you're using your own.)

## Running it

Run each step **one at a time**, from the repo root, and check its printed
output before moving on:

```bash
Rscript mvp/R/01_clean_and_stocks.R
Rscript mvp/R/01b_plot_profiles.R      # two carbon-depth profile PNGs
Rscript mvp/R/02_covariates.R          # needs Earth Engine auth, takes longest
Rscript mvp/R/03_training_data.R
Rscript mvp/R/04_train_model.R
Rscript mvp/R/05_predict_and_compare.R
Rscript mvp/R/06_bayesian_update.R
Rscript mvp/R/07_version_and_export.R
```

## The study area

`mvp/data/aoi.geojson` is the coast-following polygon that defines the mapped
extent (~5,400 km²). Step 02 uses it if present and only falls back to a
buffered hull around the cores if it's missing, so **edit that file to change
the mapping extent** — no code change needed. Step 02 stops with an error if
any core falls outside it, since a core outside the AOI silently loses all its
covariates.

Then the external-comparison add-on:

```bash
Rscript mvp/R/08_external_ingest.R       # no GEE needed
Rscript mvp/R/09_external_ecosystem.R    # needs GEE, ~11,500 points, chunked
Rscript mvp/R/10_comparison_outputs.R    # figures + GeoPackage
```

Once you trust the whole chain, `Rscript mvp/run_all.R` runs everything in
order.

Every step reads from and writes to `mvp/outputs/current/` — that folder is
always "the latest state." Step 7 archives a full copy of it into
`mvp/outputs/versions/carbon_map_v1/` (then `v2/`, `v3/`, ...) with a
`metadata.csv` describing what changed.

## What to check at each step

| Step | What "it worked" looks like |
|---|---|
| 01 | Prints a table of 8 cores with `stock_kgm2` roughly 3–15 kg/m². `cores_clean.geojson` opens in QGIS/geojson.io. |
| 01b | Two PNGs in `outputs/current/figures/`. Peat cores should show low carbon *density* but be organic-rich; the harmonised plot should show PM-01-A and PM-02-A ending early as dashed lines. |
| 02 | Prints a pre-export check confirming every band has a value at all 8 cores, then downloads `predictors.tif` (9 bands), `prior_li2025.tif`, `prior_sothe_0_30.tif`. Slowest step — several minutes. |
| 03 | `training_data.csv` has 8 rows, one column per covariate, no unexpected `NA`s. |
| 04 | Prints leave-one-core-out RMSE and R². **Expect a weak score** — 8 cores is a small sample and the script says so rather than hiding it. This is expected, not broken. |
| 05 | `carbon_prediction.tif` and `carbon_uncertainty.tif` cover the AOI; `prediction_residuals.csv` shows 8 rows of observed-minus-prior. |
| 06 | `carbon_posterior_mean.tif` should look like the prior almost everywhere, with a visible pull toward the observed value near each core. `carbon_difference_from_prior.tif` should be ~0 far from every core. |
| 07 | `outputs/versions/carbon_map_v1/metadata.csv` exists with today's date and 8 cores. |

## Running it again next field season

1. Add new rows to `../data/raw/community_soil_cores.csv` (or point
   `CFG$file_cores_raw` at a new file with the old + new cores combined).
2. In `config.R`, set:
   ```r
   bayes = list(prior_source = "previous_version", ...)
   ```
3. Run steps 01–07 again. Step 5 will now load `carbon_map_v1`'s posterior
   as the prior instead of Li et al., and step 6 fuses in only the new
   cores' residuals against it. Step 7 writes `carbon_map_v2/`.

## The external comparison (steps 08–10, an add-on)

Steps 01–07 produce the map. Steps 08–10 answer a different question — how do
these eight cores sit against everything else that has been measured — and
nothing in 01–07 depends on them.

Four databases, harmonised onto one schema in step 08:

| Dataset | Profiles | Layers | What it is |
|---|---|---|---|
| CanPeat | 1,217 | 37,072 | Canadian peatlands |
| NPDB | 9,017 | 48,372 | AAFC National Pedon Database — the only source with a mineral/organic flag |
| Janousek | 1,284 | 23,018 | US Pacific + Gulf tidal wetlands |
| WOSIS Canada | 29 | 124 | Canadian subset only |

**The global WOSIS table is deliberately excluded.** It is 14,596 profiles, 92%
United States, and would swamp every comparison with irrelevant geography.

**Two stock columns, never interchangeable.** `stock_kgm2_total` covers
whatever depth each profile actually reached (CanPeat's median is 243 cm; ours
is 30 cm). `stock_kgm2_0_30` covers 0–30 cm only, integrated from the layers
with partial layers apportioned by overlap. Comparing our 30 cm total against a
243 cm total is a category error, so both exist and every figure says which it
uses.

For the GeoPackage this matters for symbology: size-by-`stock_kgm2_total` makes
our cores (3–14) nearly invisible beside CanPeat (median 119), not because they
hold less carbon but because they are a tenth of the depth. Use
`stock_kgm2_0_30` for a like-for-like size comparison.

### The finding

Inside the Hudson & James Bay Lowlands, across all four databases combined:
**86 external profiles, and every one of them is organic.** Zero mineral soil
profiles. The nearest profile anywhere with a complete, comparable 0–30 cm
stock is **700 km away**. Step 08 recomputes these numbers on every run and
writes them to `narrative_stats.csv`, so anything quoted from them is traceable
rather than remembered.

## Known simplifications (intentional, for this first pass)

These are all documented as deliberate MVP trade-offs in the design report
(`.claude/brainstorm-mvp.md`), not oversights:

- **One depth window** (however deep each core actually reached), not the
  original's dual common-support/reference-window system. Cores under
  30 cm are flagged `stock_is_lower_bound = TRUE`, never rescaled.
- **One carbon-conversion factor**, not the dual as-supplied/harmonised-peat
  scenario comparison.
- **One random forest**, cross-validated honestly (leave-one-core-out), not
  a 5-candidate tournament with pre-registered selection rules and
  permutation-null significance testing.
- **A smaller covariate stack** (elevation, slope, Sentinel-1 VV/VH,
  Sentinel-2 NDVI, JRC water occurrence) — distance-to-coast and the
  wetland/land-cover strata layers are dropped for now.
- **A simple buffered convex hull AOI**, not the coast-oriented,
  PCA-rotated rectangle.
- **A spatially constant prior uncertainty** (35 kg C/m², Li et al.'s
  published regional figure) rather than a per-pixel uncertainty raster.
- **No QA/QC ceremony** — the three real data errors (longitude sign, and
  the two the original found in unit handling) are one-line `dplyr`
  checks; the 5-table audit trail is gone.

Add any of these back deliberately once the pipeline is running smoothly
and you understand every line — that was the explicit brief for this first
pass, not a limitation to work around.

## Start each run in a clean R session

If you `source()` these from an R session that also has the original
workflow's `R/` files loaded, you'll see masking warnings like
`yardstick ... masked _by_ '.GlobalEnv': mae, metric_set, rmse` and
`terra ... masked _by_ '.GlobalEnv': draw, zonal`. `.GlobalEnv` wins those
conflicts, so the old workflow's functions would shadow the packages' —
harmless as the code stands (it calls `rmse_vec()`, not `rmse()`), but it's
a trap waiting to spring. Restart R before a run.

## If a step fails

Each script's header comment lists its `INPUT`/`OUTPUT` files. Almost every
failure will be one of:

- A missing package → rerun `mvp/install_packages.R`.
- A missing input file → the previous step didn't finish or didn't write
  where this step expects. Check `outputs/current/` for what's actually
  there.
- Earth Engine auth (step 02 only) → `rgee::ee_Initialize()` again.

### Earth Engine gotcha worth knowing

`ERROR in Earth Engine servers: Exported bands must have compatible data
types; found inconsistent types: Float32 and Float64.`

Earth Engine will not export a multi-band image whose bands have different
data types, and a covariate stack naturally mixes them: `ee$Terrain$slope()`
returns Float32, `median()` and `normalizedDifference()` return Float64, and
JRC water occurrence is uint8. Step 02 calls `$toFloat()` on the finished
stack to make them uniform. **If you add another covariate band later, the
`$toFloat()` at the end of the `predictors <- ...` chain is what keeps the
export working** — don't drop it.

Paste the error and the last few lines of console output back and we'll
fix that one step without touching the others.
