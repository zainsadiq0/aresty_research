# Cinder Cone Morphometry

Automated extraction of volcanic cone shape parameters — basal diameter, crater
diameter, height, steepness, flatness, and volume — from DEM elevation profiles,
using slope-based break-in-slope detection rather than fixed profile endpoints.

Developed as Aresty undergraduate research at Rutgers University under
Prof. K. G. Bemis, applied to cinder cones in the San Francisco Volcanic Field,
Arizona.

---

## The problem

Measuring a cinder cone sounds trivial and isn't. Every shape parameter depends
on one decision: **where does the cone end and the surrounding terrain begin?**

The naive approach takes the endpoints of a profile line as the cone's base. That
answer is set by where the analyst happened to stop drawing, not by the landform,
and it overestimates basal diameter, sometimes badly. Everything downstream
(steepness, flatness, volume) inherits the error.

The harder problem is that real cones don't sit on flat ground. On the cone
measured here, the two ends of a single profile differ in elevation by up to
**123 m**. Averaging them into one "base elevation" silently discards a real
property of the landform.

## Method

The pipeline went through three iterations, each driven by a failure of the last.

### 1. Profile endpoints — rejected

Basal diameter taken as `dist[-1] - dist[0]`. Fails as described above: the
measurement reports the analyst's line length, not the cone.

### 2. Profile curvature — rejected

The toe of a cone is a *bend*, not a high or low point, so it is invisible to a
peak finder on elevation but appears as a spike in the second derivative. This
works, but proved fragile: basal diameter moved **170 m** on one profile from
changing the smoothing window alone (15 → 21 points), because over-smoothing
erases the very bend being detected. The controlling parameter — a window width —
also has no physical meaning to defend.

### 3. Slope threshold — current

A cone flank sits near the angle of repose; the surrounding plain does not. On
this cone the separation is roughly fivefold:

| | Slope |
|---|---|
| Surrounding terrain | 2.6 – 7.1° |
| Cone flank (median) | 25.9 – 27.9° |
| Cone flank (max) | 38.4 – 40.6° |

Toe detection walks outward from each crater rim until slope stays below a
threshold for three consecutive points. The persistence requirement matters: a
single noisy DEM cell partway up a flank can dip below threshold and stop a
naive single-crossing test far above the true base.

Three things make this better than curvature:

- **More stable.** Basal diameter moves ~40–75 m across thresholds of 6–12°,
  versus 170 m for a smoothing-window change under the curvature method.
- **Physically defensible.** "10°, chosen because it lies in the empty gap between
  two distinct slope populations" is a methods-section sentence. "A 15-point
  smoothing window" is not.
- **Visible.** The cone is unmistakable on a slope raster while it blends into
  surrounding terrain on a hillshade — the method matches how the landform
  actually presents.

### Base elevation

Rather than averaging the two toe elevations, a first-order surface is fitted by
least squares to the terrain *outside* both detected toes and evaluated at the
cone centre. On a tilted regional surface this recovers the elevation the ground
would have had without the cone. Per-flank heights are reported alongside the
single value so asymmetry stays visible rather than being averaged away.

---

## Validation

Measured cone: S P Crater (S P Mountain), San Francisco Volcanic Field, Arizona.
USGS 1/3 arc-second (~10 m) DEM, four radial profiles.

**Against published values:**

| | This method | Published | Source |
|---|---|---|---|
| Crater diameter | 421.8 m (mean) | 400 m | Ulrich & Bailey (1987) |
| Crater depth | ~110 m | ~120 m | Ulrich & Bailey (1987) |
| Summit elevation | 2138.7 m | 2140 m | NGVD 29 benchmark |

Crater diameter runs ~5% high, which is expected — crater measurements are more
sensitive to DEM resolution than cone measurements because the feature is
smaller.

**Internal consistency checks:**

- **Crater floor elevation** recovered independently from four transects:
  2013.6 / 2016.5 / 2015.7 / 2017.6 m — a 4.0 m spread. The profiles genuinely
  intersect the same crater.
- **Median flank slope** across the same four transects: 25.9 – 27.9°, consistent
  with the angle of repose for unconsolidated scoria.

**A finding, not an error:** published sources give S P Crater a height of ~250 m.
This method returns a mean height of 208 m above the fitted regional surface — but
per-flank heights on the downhill sides come out at 268.7, 252.0, and 241.1 m.
The published figure is maximum relief measured from the low side of the plain;
the value here is mean height above a fitted surface. On a cone this asymmetric
the two quantities are genuinely different, and reporting a single height obscures
that.

**Hand check:** all measurements were independently reproduced by manual
inspection of the raw elevation table. Three of four basal diameters matched to
the decimal. Remaining differences traced to two specific causes — smoothed-array
rim elevations shaving peaks by ~4 m, and the fitted-surface-versus-toe-average
base definition — rather than to general imprecision.

---

## Usage

Requires ArcGIS Pro's bundled Python environment (arcpy, numpy, scipy, pandas,
matplotlib are all present by default). Run from the Pro Python window
(View → Python).

1. Generate a Stack Profile in ArcGIS Pro from your DEM and profile lines
   (Stack Profile tool, 3D Analyst).
2. Set `TABLE_NAME` and `OUTPUT_DIR` at the top of the script.
3. Run.

Outputs to `OUTPUT_DIR`:

| File | Contents |
|---|---|
| `profile_chart.png` | Elevation profiles, all lines overlaid |
| `diagnostic_profileN.png` | Per-profile elevation + slope with detected toe and rim marked |
| `cone_table.png` | Formatted morphometry table |
| `cone_morphology.csv` | Same table as CSV |

### Parameters

| Parameter | Default | Notes |
|---|---|---|
| `SLOPE_THRESHOLD` | 10.0° | Working range ~6–12°. Above ~12° the threshold starts cutting into the lower flank. |
| `PERSIST_PTS` | 3 | Consecutive sub-threshold points required to declare the toe. |
| `SMOOTH_WINDOW` | 15 | Odd. Still required — slope is a derivative. |
| `RIM_PROMINENCE` | 10 m | Minimum prominence for crater rim detection. |
| `USE_FITTED_BASE` | True | False falls back to toe-elevation averaging. |

The script prints a threshold sensitivity sweep (6/8/10/12°) on every run, so
parameter stability is reported rather than assumed.

### Self-diagnostics

Rather than silently returning numbers, the script warns when its own assumptions
break:

- a toe walk that never found terrain below threshold (the profile line does not
  extend past the cone)
- terrain outside a detected toe within 80% of the threshold (toe position poorly
  constrained)
- flank heights differing by more than 40 m (the cone sits on a sloping surface —
  report per-flank values)

---

## Where this method works, and where it doesn't

Optimal on a young, unbreached, radially symmetric scoria cone with a
well-defined summit crater, flanks near the angle of repose, on a planar surface,
isolated from neighbouring edifices, and large relative to DEM cell size.

It degrades or fails on:

- **Craterless cones.** Rim detection requires two prominent maxima; a single-peak
  profile is skipped rather than mismeasured.
- **Eroded cones.** Flanks at 15–18° leave almost no gap above the threshold, and
  the 6–12° working window collapses. Relevant across the San Francisco Volcanic
  Field, which spans ~6 Myr of cone ages.
- **Breached cones.** The breach side has no rim peak, and any lava apron extends
  the apparent toe arbitrarily.
- **Elongated cones.** The frustum volume formula assumes one basal radius.
  Fissure-aligned cones need the planform approach below.
- **Coalesced cones.** The toe on a shared flank never reaches flat terrain.
- **Small cones.** At 10 m spacing the default smoothing window spans ~150 m. Below
  roughly 600–800 m basal diameter the window swallows the flank; reduce it
  proportionally.

## Roadmap

- **Standalone CSV input**, so the pipeline runs without an ArcGIS Pro licence.
- **2D raster delineation.** Threshold the slope raster, region-group, extract the
  cone as a polygon. Gives planform geometry (equivalent circular diameter,
  best-fit ellipse major/minor axes, circularity) and — more importantly — volume
  by integrating elevation above a fitted trend surface cell by cell, dropping the
  symmetric-frustum assumption entirely.
- **Raw-window rim elevation.** Rim *location* from the smoothed curve, rim
  *elevation* from raw data in a ±2 point window, removing a known ~4 m
  peak-shaving bias.
- Batch processing across multiple cones.

## Repository structure

```
cinder-cone-morphometry/
├── README.md
├── cone_morphometry.py          # main pipeline
├── requirements.txt
├── data/
│   └── sp_crater_profiles.csv   # sample stack profile output
├── figures/                     # generated charts and diagnostics
└── docs/
    └── method_notes.md          # extended derivation and parameter rationale
```

## References

Ulrich, G.E., and Bailey, N.G., 1987, *Geologic map of the SP Mountain part of
the San Francisco Volcanic Field, north-central Arizona*: U.S. Geological Survey
Miscellaneous Field Studies Map MF-1956. https://doi.org/10.3133/mf1956

U.S. Geological Survey Astrogeology Science Center, 2021, *A geologic field guide
to S P Mountain and its lava flow, San Francisco Volcanic Field, Arizona*.
https://www.usgs.gov/publications/a-geologic-field-guide-s-p-mountain-and-its-lava-flow-san-francisco-volcanic-field

*A review of techniques for characterising scoria cone morphologies*, Frontiers in
Earth Science, 2025. https://doi.org/10.3389/feart.2025.1667680

Elevation data: USGS 3D Elevation Program (3DEP), 1/3 arc-second DEM.

## Acknowledgments

Research conducted through the Aresty Research Center, Rutgers University, under
the supervision of Prof. K. G. Bemis, whose critique of the endpoint-based
approach prompted the move to derivative-based edge detection.

AI assistance (Claude) was used for code implementation and refactoring. The
methodology, the progression from endpoints to curvature to slope, the fitted
base surface, and the validation approach was developed by the author in
consultation with Prof. Bemis.

## License

MIT
