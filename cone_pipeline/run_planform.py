#this is the script you run to run all other scripts for cone_planform
#docstring section
"""
cone_planform.py — 2D planform morphometry of a volcanic cone from a DEM.

This script delineates the whole basal outline from a slope raster and
measures the edifice as a single object.

Parameters like ellipticity, irregularity and planform area
cannot be obtained from profiles at all.

Parameters follow Grosse et al. (2012) / MORVOLC where applicable, adapted for
scoria cones.

INPUT SECTION :
    * ArcGIS Pro raster (set DEM_SOURCE to the layer name)
    * ESRI ASCII grid .asc   (Pro: Raster to ASCII — no extra dependencies)
    * GeoTIFF .tif           (requires rasterio)
    * NumPy .npy             (set CELL_SIZE manually)

OUTPUT — printed table, CSV, polygon as .asc mask + GeoJSON, 4-panel figure.
import csv
import importlib
import math
import os

import numpy as np
from scipy import ndimage
from skimage import measure

import dem_io, terrain, delineation, metrics, export

# While developing in ArcGIS Pro: force Python to re-read edited files
for _m in (dem_io, terrain, delineation, metrics, export):
    importlib.reload(_m)

from dem_io import load_dem, rowcol_from_xy
from terrain import slope_degrees, fit_base_surface, volume_above_base
from delineation import delineate, find_crater
from metrics import planform_metrics, slope_by_elevation
from export import mask_to_geojson, overview_figure
"""
import csv
import importlib
import math
import os

import numpy as np
from scipy import ndimage
from skimage import measure

import dem_io, terrain, delineation, metrics, export

# While developing in ArcGIS Pro: force Python to re-read edited files
for _m in (dem_io, terrain, delineation, metrics, export):
    importlib.reload(_m)

from dem_io import load_dem, rowcol_from_xy
from terrain import slope_degrees, fit_base_surface, volume_above_base
from delineation import delineate, find_crater
from metrics import planform_metrics, slope_by_elevation
from export import mask_to_geojson, overview_figure

# CONFIGURATION SECTION - USER NEEDS TO EDIT 3 VARIABLES
#DEM_SOURCE
#OUTPUT_DIR
#SEED

#CHANGE DEM_SOURCE WHEN USING OWN DATA
DEM_SOURCE = r"C:\Users\zs424\Box\aresty_research\aresty_week4\aresty_week4.gdb\colton_utm"

#CHANGE OUTPUT_DIR WHEN USING OWN DATA
OUTPUT_DIR = r"C:\Users\zs424\Box\aresty_research\outputs"

CELL_SIZE = None        

SLOPE_THRESHOLD = 10.0    # degrees — same parameter as the profile pipeline
SEED = None               # (x, y) map coords inside the cone. None = highest cell.

MIN_FEATURE_CELLS = 25    # speckle smaller than this is removed before labelling
RING_WIDTH_CELLS = 10     # width of the collar outside the cone used for the base fit
TREND_ORDER = 1           # 1 = tilted plane (Grosse's polynomial base), 2 = quadratic
SINK_TOLERANCE = 1.0      # m — minimum depth for a closed depression to count as crater

MAX_PLAUSIBLE_DIAMETER = 5000.0   # m — guards against merging with a neighbour

#----------------------------------
#DRIVER SECTION
#-----------------------------------
def analyse(dem, cell_size, origin, seed_xy=None):
    slope = slope_degrees(dem, cell_size)

    if seed_xy is None:
        seed_rc = np.unravel_index(np.nanargmax(dem), dem.shape)
    else:
        seed_rc = rowcol_from_xy(seed_xy, origin, cell_size, dem.shape)

    mask = delineate(slope, seed_rc, SLOPE_THRESHOLD, MIN_FEATURE_CELLS)
    pf = planform_metrics(mask, cell_size)

    if pf["equivalent_diameter_m"] > MAX_PLAUSIBLE_DIAMETER:
        print(f"  WARNING - delineated region is {pf['equivalent_diameter_m']:.0f} m across. "
              f"It has probably merged with an adjacent high-slope feature "
              f"(neighbouring cone, flow front, scarp). Check the outline against "
              f"the hillshade before trusting any of these numbers.")

    base, dip, rms = fit_base_surface(dem, mask, cell_size, RING_WIDTH_CELLS, TREND_ORDER)
    vol_clip, vol_net = volume_above_base(dem, base, mask, cell_size)

    # `summit` is the single highest cell, which is an extreme-value statistic
    # over thousands of noisy cells and therefore biases HIGH — the mirror of
    # the profile pipeline's smoothed-rim bias, which runs low. The 99.9th
    # percentile is the robust alternative; on a synthetic cone it halves the
    # height error. Both are reported so the bias direction stays visible.
    in_cone = np.where(mask, dem, np.nan)
    summit = float(np.nanmax(in_cone))
    summit_robust = float(np.nanpercentile(in_cone, 99.9))
    rr, cc = np.nonzero(mask)
    cen_r, cen_c = rr.mean(), cc.mean()
    base_at_centre = float(base[int(round(cen_r)), int(round(cen_c))])
    height = summit - base_at_centre

    edge = mask & ~ndimage.binary_erosion(mask, np.ones((3, 3), bool))
    toe_lo = float(np.nanmin(np.where(edge, dem, np.nan)))
    toe_hi = float(np.nanmax(np.where(edge, dem, np.nan)))

    cr = find_crater(dem, mask, SINK_TOLERANCE)
    if cr is not None:
        crater = cr["mask"]
        crater_area = float(crater.sum()) * cell_size ** 2
        crater_d = 2.0 * math.sqrt(crater_area / math.pi)
        cprops = measure.regionprops(crater.astype(int))[0]
        crater_major = float(cprops.axis_major_length) * cell_size
        crater_minor = float(cprops.axis_minor_length) * cell_size
        rim_area = float(cr["rim_region"].sum()) * cell_size ** 2
        crater_d_rim = 2.0 * math.sqrt(rim_area / math.pi)
        crater_depth = cr["depth_closed"]
        crater_depth_rim = cr["depth_to_rim"]
        crater_floor = cr["floor_elev"]
    else:
        crater = None
        crater_area = crater_d = crater_major = crater_minor = float("nan")
        crater_d_rim = crater_depth = crater_depth_rim = crater_floor = float("nan")

    bands = slope_by_elevation(dem, slope, base, mask)
    flank = slope[mask & (slope > SLOPE_THRESHOLD)]

    Wco = pf["equivalent_diameter_m"]
    results = dict(pf)
    results.update({
        "summit_elev_m": summit,
        "summit_elev_robust_m": summit_robust,
        "base_elev_at_centre_m": base_at_centre,
        "height_m": height,
        "height_robust_m": summit_robust - base_at_centre,
        "max_relief_m": summit - toe_lo,
        "min_relief_m": summit - toe_hi,
        "toe_elev_range_m": toe_hi - toe_lo,
        "regional_dip_deg": dip,
        "base_fit_rms_m": rms,
        "crater_area_m2": crater_area,
        "crater_equiv_diameter_m": crater_d,
        "crater_equiv_diameter_rim_m": crater_d_rim,
        "crater_major_m": crater_major,
        "crater_minor_m": crater_minor,
        "crater_floor_elev_m": crater_floor,
        "crater_depth_m": crater_depth,
        "crater_depth_to_rim_m": crater_depth_rim,
        "volume_m3": vol_clip,
        "volume_net_m3": vol_net,
        "flank_slope_mean_deg": float(np.mean(flank)) if flank.size else float("nan"),
        "flank_slope_median_deg": float(np.median(flank)) if flank.size else float("nan"),
        "flank_slope_max_deg": float(np.max(flank)) if flank.size else float("nan"),
        # Literature-standard ratios (Porter 1972; Wood 1980)
        "H_over_Wco": height / Wco if Wco else float("nan"),
        "H_over_Wco_max_relief": (summit - toe_lo) / Wco if Wco else float("nan"),
        "Wcr_over_Wco": crater_d / Wco if Wco else float("nan"),
    })
    return results, mask, crater, slope, base, bands


def main():
    dem, cell_size, origin, _ = load_dem(DEM_SOURCE, CELL_SIZE)
    print(f"DEM {dem.shape[0]}x{dem.shape[1]} cells at {cell_size:.2f} m")

    res, mask, crater, slope, base, bands = analyse(dem, cell_size, origin, SEED)

    groups = [
        ("Planform", ["planform_area_m2", "perimeter_m", "equivalent_diameter_m",
                      "ellipse_major_m", "ellipse_minor_m", "ellipse_azimuth_deg"]),
        ("Shape indices", ["ellipticity_index", "dissection_index",
                           "irregularity_index", "circularity"]),
        ("Relief", ["summit_elev_m", "summit_elev_robust_m", "base_elev_at_centre_m",
                    "height_m", "height_robust_m",
                    "max_relief_m", "min_relief_m", "toe_elev_range_m",
                    "regional_dip_deg", "base_fit_rms_m"]),
        ("Crater", ["crater_equiv_diameter_m", "crater_equiv_diameter_rim_m",
                    "crater_major_m", "crater_minor_m", "crater_floor_elev_m",
                    "crater_depth_m", "crater_depth_to_rim_m", "crater_area_m2"]),
        ("Volume", ["volume_m3", "volume_net_m3"]),
        ("Slope", ["flank_slope_mean_deg", "flank_slope_median_deg",
                   "flank_slope_max_deg"]),
        ("Standard ratios", ["H_over_Wco", "H_over_Wco_max_relief", "Wcr_over_Wco"]),
    ]
    for title, keys in groups:
        print(f"\n{title}")
        for k in keys:
            v = res[k]
            s = f"{v:.4g}" if abs(v) < 1e6 else f"{v:.4e}"
            print(f"  {k:<28} {s:>14}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    import csv
    with open(os.path.join(OUTPUT_DIR, "cone_planform.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "value"])
        for k, v in res.items():
            w.writerow([k, v])

    mask_to_geojson(mask, cell_size, origin,
                    os.path.join(OUTPUT_DIR, "cone_basal_outline.geojson"))
        
    overview_figure(dem, slope, mask, crater, base, res, bands, cell_size,
                    os.path.join(OUTPUT_DIR, "planform_overview.png"),
                    SLOPE_THRESHOLD)

    if bands:
        with open(os.path.join(OUTPUT_DIR, "slope_by_height.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(bands[0].keys()))
            w.writeheader()
            w.writerows(bands)

    print(f"\nWrote CSV, GeoJSON outline and overview figure to {OUTPUT_DIR}")
    print("Check the outline against the hillshade panel before using any number.")


if __name__ == "__main__":
    main()
