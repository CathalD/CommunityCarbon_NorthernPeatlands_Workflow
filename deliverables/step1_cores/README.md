# Fort Severn community soil cores — GIS layers

Prepared 29 July 2026 for the Fort Severn First Nation governing council.

## What this is

The complete record of 8 soil cores collected by community members across two field seasons, as GIS point layers. Nothing has been removed. Values that looked doubtful are kept and marked with a flag, so you can see exactly what was measured and what we thought about it.

## Files

| File | What it is |
|------|------------|
| `fort_severn_cores.gpkg` | GeoPackage, opens in ArcGIS and QGIS. Two layers: `cores_qc` (one point per core) and `segments_raw` (one point per sampled slice). |
| `fort_severn_cores_shp.zip` | The same data as Shapefiles, zipped for uploading to Google Earth Engine. |
| `*.csv` | Plain tables. Open in Excel or anything else. |
| `*.geojson` | Web-mapping format. |
| `FIELD_DICTIONARY.csv` | Every column explained, with its units. |

**Coordinate system: EPSG:3978 (NAD83 / Canada Atlas Lambert).** Latitude and longitude are also included in every file.

## The two layers

**`cores_qc`** — one point per core (8 points). Use this for maps and for totals.

**`segments_raw`** — one point per slice cut from a core (22 points). Several slices share a location because they came from the same hole, one below the other. Use this to look at how carbon changes with depth.

## Two numbers for carbon, and why

Cores stopped at different depths, between 14.5 and 30.7 cm. To compare them fairly you must compare the same depth in each.

- **`soc_shal`** — carbon over the top 14.5 cm, the deepest layer *every* core reached. Every core has a real number here, so this is the fair comparison.
- **`soc_30cm`** — carbon over the top 30 cm, the depth published maps use. Two cores never got that deep. For those, the number is a **minimum**, not a measurement, and `soc30_lb = 1` marks them.

Never mix the two columns in one total.

## What the flags mean

| Flag | Meaning |
|------|---------|
| `LON_SIGN_REPAIRED` | Longitude sign was wrong for the stated study area and was repaired. The supplied value would place this core in the eastern hemisphere. |
| `BD_ABOVE_MAX` | Dry bulk density exceeds the plausible ceiling for soil; approaches parent-rock density. Suspect a sample-volume or compaction error. |
| `SEGMENT_INDEX_ABSENT` | Sample Id carries no trailing integer, so segment order had to be assigned rather than read. |
| `DUPLICATE_LAB_VALUES` | This sample shares an identical organic matter AND organic carbon pair with a sample from a DIFFERENT core. Possible transcription or copy-paste error. |
| `SINGLE_SEGMENT_CORE` | Core is represented by a single segment, so within-core variability cannot be estimated for it. |

## Two things you should know about the data

**1. The longitudes were recorded without a minus sign.** As supplied, the coordinates placed these cores in Siberia, about 7,540 km away. This was corrected on loading, and every corrected row is marked `lon_fixed = 1`. After correction the cores sit 2.7 to 19.4 km from Fort Severn, which confirms the correction was right.

**2. The depth ordering is an assumption that needs a field check.** The `Depth` column records the THICKNESS of each slice, not its depth below the surface. We established this from the data itself: read as depth-below-surface the values would decrease going down two cores, which is impossible. Slice order comes from the number at the end of the sample name.

> **Please confirm against the field notes** that slice `PM-03-A-1` is the surface slice and `PM-03-A-4` the deepest, and likewise for the other cores. Every depth in these files, and every carbon total, depends on that ordering being right.

## How carbon was calculated

```
carbon in a slice  =  bulk density  x  carbon %  /  100  x  slice thickness
     (kg C/m2)          (g/cm3)                              (cm)

core total  =  sum of all its slices
```

One caution on the carbon percentages: the two field seasons used different laboratory conversion factors to turn organic matter into organic carbon (0.58 for the 2024 peat cores, 0.554 for the 2025 mineral cores). This means a peat-versus-mineral comparison is partly a laboratory difference, not only a difference in the ground. Using one method for both, or measuring carbon directly on a few samples, would settle it.

## Questions

Every number in these files is produced by open scripts from the original field data, and can be regenerated from scratch. If a value looks wrong, it can be traced back to the exact measurement it came from.
