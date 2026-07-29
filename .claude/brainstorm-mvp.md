# Brainstorm: Simplified Community Carbon Mapping MVP

## Context
Full workflow: 18 scripts, Earth Engine integration, national database context, Bayesian mapping, complex spatial analysis.

**Practitioner MVP Goal:** Core carbon accounting & basic spatial output without GIS expertise.

---

## MVP SCOPE: What to Keep, What to Cut

### ✅ KEEP (Core Value)
- Raw data ingestion (8 cores, 22 segments)
- Quality control (catch obvious errors, flag suspects)
- Carbon stock calculation (bulk density × SOC% × depth)
- Depth standardization (common support window)
- Stratified estimate (peat vs. mineral mean + uncertainty)
- Core locations map (GeoJSON)
- Community summary figures (4 simple PNG)

### ❌ CUT (Complexity Without Proportional Value for Amateur)
- Earth Engine covariates (needs auth, 50 km buffer logic)
- Predictive/kriging surfaces (8 cores insufficient anyway)
- Bayesian prior update (needs Li et al. peat-column data)
- Random forest covariate models
- National database context (NPDB integration)
- Embedding clusters analysis
- Coast-following AOI (axis-aligned box sufficient)
- Advanced validation ledger
- Technical Quarto reports

---

## MVP OUTPUT STRUCTURE

| Output | File | Description | Package |
|--------|------|-------------|---------|
| **1. Clean core data** | `cores_clean.csv` | 8 cores with QC flags, carbon stocks at 0–14.5 cm | Base R / data.table |
| **2. Segment depth table** | `segments_standardized.csv` | 22 segments at harmonized depths | Base R |
| **3. Spatial locations** | `core_locations.geojson` | 8 points with stocks, strata, flags | `sf` or `geojsonsf` |
| **4. Stratified estimate** | `stratum_summary.csv` | Mean ± SD, n, CI by soil type | Base R |
| **5. Spatial raster mask** | `study_area_mask.tif` | AOI boundary (simple polygon) | `terra` or `raster` |
| **6. Community brief** | `BRIEF.md` | Plain-language summary | Markdown |
| **7. Figures** | `fig_01_core_locations.png` etc. | 4 simple PNGs (locations, stocks, stratum comparison) | `ggplot2` or `base` |

---

## STEP-BY-STEP MVP PIPELINE

### Step 1: DATA INGESTION & QC

**What:** Read CSV, recompute derived columns, flag implausible values  
**Package:** `data.table` (fast CSV, no dependencies) OR Base R  
**Functions:**
- `fread()` / `read.csv()` — load raw cores
- Check column names match expected
- Validate: longitude ∈ [-88, -87], latitude ∈ [55, 56] (Fort Severn bounds)
- Flag: bulk density > 2.5 g/cm³ (soil ceiling)
- Flag: depth < 0 or > 50 cm
- Flag: missing required columns (lat, lon, bulk_density, OM, SOC%, depth)

**Output:** `cores_qc.csv` (8 rows, all kept, flagged)

---

### Step 2: CARBON STOCK CALCULATION

**What:** Compute stock from measurements; detect unit errors  
**Formula:**  
```
stock_kgm2 = bulk_density_g_cm3 × SOC_pct / 100 × depth_cm / 100
```

**Package:** Base R (vectorized arithmetic)  
**Functions:**
- Apply formula to each segment
- Compute campaign-specific OM→carbon factor (if different years supplied)
- Reconcile supplied vs. computed SOC (flag large discrepancies)
- Store both supplied and computed values

**Output:** `segments_with_stocks.csv` (22 rows, with `stock_kgm2` column)

---

### Step 3: DEPTH HARMONIZATION

**What:** Find common depth window all cores reach; rescale shorter cores as lower bounds  
**Package:** Base R  
**Functions:**
- Read actual depth reached by each core from data
- Find minimum depth reached = common support window (here: 14.5 cm)
- Aggregate segments to that window per core
- For cores shallower than 30 cm: mark `depth_30cm_is_lower_bound = TRUE`

**Output:** `cores_harmonized.csv` (8 rows, with stocks at 0–14.5 cm + 0–30 cm flag)

---

### Step 4: STRATIFIED ESTIMATE

**What:** Mean & SD of carbon stock, by soil type or campaign  
**Package:** Base R or `tidyverse` (dplyr if available, else aggregate)  
**Functions:**
- Group by campaign/stratum (peat vs. mineral)
- Calculate: mean, sd, se, n, 95% CI per group
- Boot-strap CI optional (Base R `sample()` loop if needed, ~1000 iterations)
- Output: summary table, one row per stratum

**Output:** `stratum_summary.csv` (2 rows: peat, mineral with stats)

---

### Step 5: SPATIAL: CORE LOCATIONS

**What:** Simple GIS layer of 8 core points with stocks & attributes  
**Package:** `sf` (modern, lightweight) OR `geojsonsf` (GeoJSON focused) OR manual JSON  
**Functions:**
- Create sf object from lat/lon columns
- Set CRS = EPSG:4326 (WGS84, standard; no projection needed for visual)
- Attributes: core_id, stock_0_14.5cm, stock_0_30cm, stratum, flags
- Export: GeoJSON (human-readable, opens in any web map)

**Output:** `core_locations.geojson`

---

### Step 6: SPATIAL: STUDY AREA MASK (OPTIONAL, SIMPLE)

**What:** Bounding box or convex hull around cores as a reference polygon  
**Package:** `sf` or manual polygon in GeoJSON  
**Functions:**
- Compute convex hull of core coordinates
- OR use a simple rectangular buffer (e.g., 20 km around cores)
- Export as GeoTIFF raster (1 = inside AOI, 0 = outside) OR GeoJSON polygon

**Output:** `aoi_simple.tif` OR `aoi_simple.geojson`

---

### Step 7: FIGURES (BASE GRAPHICS OR GGPLOT2)

**What:** Four publication-ready community figures  
**Package:** Base R `plot()` + `png()` OR `ggplot2` + `ggsave()`  
**Functions:**

#### Figure 1: Core Locations Map
```
plot(lon, lat, col=c("orange","blue")[as.factor(stratum)], 
     main="Soil Core Locations", pch=19, cex=2)
legend(..., c("Peat", "Mineral"))
```
**Output:** `fig_01_core_locations.png`

#### Figure 2: Carbon Stocks by Core
```
barplot(cores_harmonized$stock_0_14.5cm, 
        names=cores_harmonized$core_id,
        col=c("orange","blue")[...],
        ylab="Carbon (kg C/m²)", main="0–14.5 cm Stocks")
```
**Output:** `fig_02_core_stocks.png`

#### Figure 3: Stratum Comparison
```
boxplot(stock_0_14.5cm ~ stratum, 
        main="Peat vs. Mineral\n(0–14.5 cm)", ylab="Carbon (kg C/m²)")
```
**Output:** `fig_03_stratum_comparison.png`

#### Figure 4: Depth Profile (One Core)
```
plot(cumsum(stock_per_segment) ~ depth, 
     main="Carbon Accumulation with Depth (Example Core)", 
     type="l", xlab="Depth (cm)", ylab="Cumulative C (kg/m²)")
```
**Output:** `fig_04_depth_profile.png`

---

### Step 8: COMMUNITY BRIEF

**What:** Plain-language markdown summary answering 4 key questions  
**Package:** None (write markdown directly in R via `cat()` or `writeLines()`)  
**Functions:**
- Template: 4 headings
  1. What is the carbon stock around [location]?
  2. What did the new data teach us?
  3. How do these results compare to other areas?
  4. Conservation perspective
- Populate with numbers from steps 1–4
- Add data quality notes
- Cite methods simply

**Output:** `BRIEF.md` (shareable, no code, 1-2 page markdown)

---

## DEPENDENCIES FOR MVP

### Minimal Set (No External Packages)
- **Base R only:** Read CSV, QC, calculations, basic stats, markdown generation
- Scripts: ~3 main R files (~200 lines each)
- Execution: `Rscript main.R` (no package install ordeal)

### Lightweight Additions (Optional, Improves Output)
| Need | Package | Why | Install |
|------|---------|-----|---------|
| Fast CSV | `data.table::fread()` | 10× faster on large files | `install.packages("data.table")` |
| Modern maps | `sf` | Standard for GIS layers; GeoJSON export | `install.packages("sf")` |
| Better figures | `ggplot2` | Publication-ready graphics | `install.packages("ggplot2")` |
| Bootstrap stats | Base R `sample()` | Already included | (none) |
| Spatial raster | `terra` | AOI mask as GeoTIFF | `install.packages("terra")` |

---

## SOIL/SEDIMENT PACKAGES FOR EXTENSION

If practitioner later wants carbon estimation tools:

| Package | Function | Use Case |
|---------|----------|----------|
| `soilDB::fetchOSD()` | Fetch USDA soil series | Compare against regional soil types |
| `soilDB::fetchSDA()` | NRCS Soil Data Access | Lookup known profiles at site |
| `isoSpikeGrowth::estimateC()` | Estimate C from SOM | If direct C measurement unavailable |
| Custom `depth_integrate()` | Integrate profile to any depth | General purpose depth standardization |

---

## ALTERNATIVE ARCHITECTURE: PYTHON

If practitioner prefers Python (growing in soil science):

| Task | R Package | Python Equiv |
|------|-----------|--------------|
| Data I/O | Base R / data.table | `pandas` |
| Spatial | `sf` | `geopandas` |
| Stats | Base R | `numpy` / `scipy` |
| Figures | `ggplot2` | `matplotlib` / `plotly` |
| GeoTIFF | `terra` | `rasterio` |
| Soil DB | `soilDB` | `open-soil-data` or custom API |

**Minimum Python stack:** `pandas`, `geopandas`, `matplotlib`

---

## EXECUTION SUMMARY

```
main.R
├── 01_ingest_qc.R         [reads CSV, validates, flags]
├── 02_stocks_depth.R      [calculates stocks, harmonizes depth]
├── 03_stratum_summary.R   [groupwise statistics + CI]
├── 04_spatial.R           [writes GeoJSON + simple mask]
├── 05_figures.R           [4 PNGs + BRIEF.md]
└── run_all.R              [wrapper: sources all above]

OUTPUT TREE:
├── cores_qc.csv
├── segments_standardized.csv
├── stratum_summary.csv
├── core_locations.geojson
├── aoi_simple.geojson
├── fig_01_core_locations.png
├── fig_02_core_stocks.png
├── fig_03_stratum_comparison.png
├── fig_04_depth_profile.png
└── BRIEF.md
```

---

## WHAT'S NOT IN MVP (But Could Be Added Later)

1. **Kriging / IDW interpolation** — Not justified for 8 cores; explain why instead
2. **Bayesian prior update** — Requires external full-column peat-depth product
3. **Earth Engine integration** — Adds complexity; local data often sufficient
4. **Machine learning model** — Cross-validation over 8 cores is overfit; better as "lessons learned"
5. **National context** — Out of scope for practitioner; focus on local story
6. **Advanced visualization** — Leaflet maps nice but not essential; static PNG sufficient

---

## MVP FEASIBILITY CHECKLIST

- [ ] Can be run in base R (no external packages required)
- [ ] Produces valid GIS layers (GeoJSON, optionally GeoTIFF)
- [ ] Includes honest uncertainty (CI, not false precision)
- [ ] Flags data quality (keeps all rows, marks suspects)
- [ ] Outputs community-readable brief + figures
- [ ] Reproducible from a single `Rscript main.R` call
- [ ] Execution time < 10 seconds
- [ ] Disk footprint < 10 MB (no large rasters)
- [ ] Practitioner can modify parameters (e.g., depth windows) with 1 config edit

---

## NOTES FOR MVP DESIGN

1. **Depth window choice:** Use 0–14.5 cm (deepest all cores reach). Explain *why* not 0–30 cm in brief.
2. **Stratum definition:** Use campaign (peat vs. mineral). Simple, scientifically defensible.
3. **Null model:** Just report the mean ± SD. No interpolation.
4. **Confidence interval:** Bootstrap or t-CI, both ~2 lines of code. Report it.
5. **QC philosophy:** Catch and flag, never drop. Analyst decides.
6. **Outputs are consumable:** CSV (Excel), GeoJSON (web + GIS), PNG (email/print), Markdown (GitHub/web).
