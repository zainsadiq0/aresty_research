# Cone Planform Morphometry

Automated 2D morphometry of volcanic scoria cones from a DEM, built for research on the San Francisco Volcanic Field (SFVF), Arizona.

The pipeline outlines a cone's base from a slope raster and measures the whole edifice as one object: planform area and shape, height above a fitted regional surface, crater size and depth, volume, and slope distribution. Shape parameters follow Grosse et al. (2012) and the MORVOLC conventions where applicable, adapted for scoria cones.

It complements a profile-based approach (measuring cones along transects). Outline shape, elongation direction and an asymmetry-free volume can only be measured from the full 2D surface.

## How it works

1. **Slope.** Computes slope for every cell with the Horn (1981) method, the same algorithm as the ArcGIS Pro Slope tool.
2. **Outline.** Marks cells steeper than a threshold (default 10°), removes small speckle, fills the crater hole (the flat crater floor falls below the threshold), and keeps the connected region containing a seed point. The 10° default sits between typical SFVF plain slopes (about 1–6°) and cone flanks near the angle of repose (about 27°).
3. **Planform metrics.** Area, perimeter, equivalent diameter (Wco), best-fit ellipse with compass azimuth, and the ellipticity, dissection and irregularity indices.
4. **Base surface.** Fits a least-squares plane through a collar of terrain just outside the outline, so cones on tilted ground are measured against the regional surface rather than an averaged toe elevation.
5. **Height and volume.** Height is summit minus the base surface under the cone's centroid. Volume is integrated cell by cell above the base, with no cone or frustum shape assumed.
6. **Crater.** Located two ways, because they differ measurably:
   - *Closed depression* (MORVOLC's lowest closed contour): conservative, but a breached rim makes it small.
   - *Rim crest*: walks outward from the depression to the ring of highest mean elevation.
7. **Slope by height.** Mean and median slope in 20 normalised height bands.

## Files

| File | Purpose |
|---|---|
| `run_planform.py` | Settings (CONFIG), the analysis driver, and output writing. **This is the file you run.** |
| `dem_io.py` | Loads the DEM (ArcGIS raster, ESRI ASCII grid, GeoTIFF or NumPy array) and converts map coordinates to grid cells. |
| `terrain.py` | Slope, fitted base surface, and volume above it. |
| `delineation.py` | Cone outline from the slope raster, and crater detection. |
| `metrics.py` | Planform size and shape metrics, and slope-by-height bands. |
| `export.py` | GeoJSON outline and the four-panel overview figure. |

All six files must be in the same folder.

## Requirements

- Python 3.9+
- numpy, scipy, matplotlib, scikit-image
- arcpy (optional): to read rasters directly from an ArcGIS Pro geodatabase
- rasterio (optional): to read GeoTIFFs outside ArcGIS

### ArcGIS Pro setup

ArcGIS Pro's default Python environment includes numpy, scipy and matplotlib but not scikit-image, and it cannot be modified. Clone it once:

1. **Project → Package Manager**, then the gear icon next to *Active Environment*.
2. Clone `arcgispro-py3`, select the clone as the active environment, and restart Pro.
3. In **Package Manager → Add Packages**, install `scikit-image`.

## Preparing the DEM

The code assumes **linear units in metres**. USGS 3DEP DEMs are usually distributed in geographic coordinates (degrees), and running on them gives wrong slopes, areas and volumes without raising an error.

In ArcGIS Pro, run **Project Raster** on the elevation DEM:

- Output coordinate system: `NAD 1983 UTM Zone 12N` (covers the SFVF)
- Resampling: `Bilinear`
- Cell size: `10` for a 1/3 arc-second DEM
- Save the output inside a file geodatabase (no file extension)

Clip the DEM with a few hundred metres of surrounding terrain on every side of the cone, so the base surface has ground to fit, and keep neighbouring cones out of the extent where possible.

## Configuration

Edit the CONFIG block at the top of `run_planform.py`:

| Setting | Meaning | Default |
|---|---|---|
| `DEM_SOURCE` | Path to the projected DEM. A geodatabase raster path must end with the raster name, not `.gdb`. | required |
| `OUTPUT_DIR` | Folder for results (created if missing). | required |
| `SEED` | `(x, y)` map coordinate inside the target cone. `None` uses the highest cell. | `None` |
| `SLOPE_THRESHOLD` | Degrees; steeper cells count as cone. | `10.0` |
| `MIN_FEATURE_CELLS` | Steep patches smaller than this are removed as noise. | `25` |
| `RING_WIDTH_CELLS` | Width of the terrain collar used for the base fit (10 cells = 100 m at 10 m resolution). | `10` |
| `TREND_ORDER` | Base surface: `1` tilted plane, `2` quadratic. | `1` |
| `SINK_TOLERANCE` | Minimum depth (m) for a closed depression to count as a crater. | `1.0` |
| `MAX_PLAUSIBLE_DIAMETER` | Warns if the outline is wider than this (m), which usually means it merged with a neighbour. | `5000.0` |
| `CELL_SIZE` | Only needed for `.npy` input. | `None` |

Use raw strings for Windows paths, e.g. `r"C:\data\project.gdb\colton_utm"`.

## Running

**In the ArcGIS Pro Python window** (enter one line at a time):

```python
import sys, importlib
sys.path.insert(0, r"C:\path\to\this\repo")
import run_planform
run_planform.main()
```

To re-run after changing CONFIG in the same session:

```python
importlib.reload(run_planform)
run_planform.main()
```

**From a terminal** (with the requirements installed):

```bash
python run_planform.py
```

## Outputs

| File | Contents |
|---|---|
| `cone_planform.csv` | Every measured parameter. |
| `slope_by_height.csv` | Mean and median slope per height band. |
| `cone_basal_outline.geojson` | Basal outline polygon in the DEM's map coordinates. |
| `planform_overview.png` | Hillshade with outline and crater, slope map, relief above the base, and slope vs. height. |

**Always check the outline in `planform_overview.png` before using the numbers.** If the outline spreads onto a neighbouring cone, lava flow or scarp, the diameter, area and volume are overestimated.

### Measured parameters

- **Planform:** area, perimeter, equivalent diameter (Wco), ellipse major/minor axes and azimuth
- **Shape indices:** ellipticity, dissection, irregularity, circularity
- **Relief:** summit (highest cell and 99.9th percentile), base elevation, height, max/min relief, toe elevation range, regional dip, base-fit RMS
- **Crater:** equivalent diameter (closed depression and rim crest), major/minor axes, floor elevation, depth (closed and to rim crest), area
- **Volume:** clipped (above base only) and net
- **Slope:** flank mean, median and maximum
- **Ratios:** H/Wco, H/Wco from maximum relief, Wcr/Wco (Porter 1972; Wood 1980)

## Example: Colton Crater

10 m USGS DEM projected to UTM Zone 12N, slope threshold 10°.

| Parameter | Value |
|---|---|
| Equivalent base diameter (Wco) | 2,244 m |
| Height above fitted base | 245 m |
| Crater diameter, closed depression / rim crest | 849 m / 1,355 m |
| Crater depth, closed / to rim crest | 131 m / 201 m |
| Volume, clipped / net | 2.93 × 10⁸ m³ / 2.76 × 10⁸ m³ |
| Median flank slope | 20.1° |
| H/Wco | 0.109 |

These are preliminary. The outline currently includes two attached lobes on the cone's south side, which inflates Wco, area and volume, and Colton's rim is breached to the north, so the rim-crest crater values are the more representative ones.

## Validation

- **Synthetic cone.** A cone of known dimensions (elliptical base, tilted regional surface, correlated noise) was run through the pipeline. The two crater definitions bracket the true diameter at about −6% and +6%.
- **Cross-platform check.** The same Colton DEM processed in ArcGIS Pro and in standalone Python agreed within 0.1% for area, height, crater and volume. Perimeter differed by about 4% because of small differences between the two projected grids, so perimeter-based indices are the least stable outputs.

## Known limitations

- Adjacent landforms steeper than the threshold merge into the outline. An option to exclude areas with a user-supplied polygon is planned.
- With `SEED = None`, the highest cell is used, which selects the wrong cone if a taller neighbour is in the extent.
- NoData cells are filled with the mean elevation before slope is computed, creating artificial steep edges near NoData areas.
- Shape indices are computed for the basal outline only, not averaged over every flank contour as in Grosse et al. (2012), so values are not directly comparable with theirs.
- The GeoJSON outline has no CRS tag and is offset by half a cell.
- The GeoTIFF reader does not convert NoData values to NaN; use a geodatabase raster or ASCII grid until this is fixed.

## References

- Grosse, P., van Wyk de Vries, B., Euillades, P. A., Kervyn, M., and Petrinovic, I. A. (2012). Systematic morphometric characterization of volcanic edifices using digital elevation models. *Geomorphology*, 136, 114–131.
- Horn, B. K. P. (1981). Hill shading and the reflectance map. *Proceedings of the IEEE*, 69(1), 14–47.
- Porter, S. C. (1972). Distribution, morphology, and size frequency of cinder cones on Mauna Kea volcano, Hawaii. *Geological Society of America Bulletin*, 83, 3607–3612.
- Wood, C. A. (1980). Morphometric evolution of cinder cones. *Journal of Volcanology and Geothermal Research*, 7, 387–413.

## Author

Zain Sadiq, Rutgers University. Developed as part of Aresty Undergraduate Research under Prof. Karen Bemis. Written with assistance from Claude (Anthropic).
