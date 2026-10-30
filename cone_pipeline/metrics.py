"""Planform shape metrics and slope-by-height bands."""
import math
import numpy as np
from skimage import measure

#planform metrics
def planform_metrics(mask, cell_size):
    """Area, perimeter, equivalent diameter, best-fit ellipse, shape indices.

    Index definitions (Grosse et al. 2012):
        ellipticity  ei = pi * (major/2)^2 / A        1.0 for a circle, a/b for an ellipse
        dissection   di = P / (2 * sqrt(pi * A))      1.0 for a circle, >1 as outline convolutes
        irregularity ii = di(outline) - (di(ellipse) - 1)
                                                      1.0 for a circle or any ellipse
        where di(ellipse) is the dissection index of an ellipse whose ellipticity
        index equals the outline's (ei_ellipse = ei_outline), per Grosse et al.
        (2012, eq. on p. 121; after Grosse et al. 2009).

    Grosse et al. compute ei and ii for every closed flank contour and report the
    averages (eiAVG, iiAVG). Here they are computed for the basal outline only.
    """
    area = float(mask.sum()) * cell_size ** 2

    # Crofton perimeter is far less biased than counting staircase cell edges,
    # which overestimates by up to ~27% on a diagonal boundary.
    perim = float(measure.perimeter_crofton(mask, directions=4)) * cell_size

    props = measure.regionprops(mask.astype(int))[0]
    major = float(props.axis_major_length) * cell_size
    minor = float(props.axis_minor_length) * cell_size
    # skimage measures orientation counter-clockwise from the row axis, which
    # points SOUTH on a north-up raster. Negating it gives a compass azimuth
    # (deg clockwise from north, mod 180). Checked on synthetic ellipses at
    # 0-150 deg. Elongation azimuth is what tests fissure/vent alignment.
    orient = (-math.degrees(float(props.orientation))) % 180.0

    equiv_d = 2.0 * math.sqrt(area / math.pi)
    ei = math.pi * (major / 2.0) ** 2 / area if area > 0 else float("nan")
    di = perim / (2.0 * math.sqrt(math.pi * area)) if area > 0 else float("nan")

    # Dissection of an ellipse with the SAME ellipticity index as the outline.
    # For an ellipse ei = a/b, so b = a/ei; di is scale-free, so a is arbitrary.
    # Ramanujan perimeter approximation.
    a = 1.0
    b = a / ei if (ei and math.isfinite(ei) and ei > 0) else a
    h = ((a - b) ** 2) / ((a + b) ** 2) if (a + b) > 0 else 0.0
    p_ell = math.pi * (a + b) * (1 + 3 * h / (10 + math.sqrt(4 - 3 * h)))
    a_ell = math.pi * a * b
    di_ell = p_ell / (2.0 * math.sqrt(math.pi * a_ell)) if a_ell > 0 else float("nan")

    return {
        "planform_area_m2": area,
        "perimeter_m": perim,
        "equivalent_diameter_m": equiv_d,
        "ellipse_major_m": major,
        "ellipse_minor_m": minor,
        "ellipse_azimuth_deg": orient,
        "ellipticity_index": ei,
        "dissection_index": di,
        "irregularity_index": di - (di_ell - 1.0),
        "circularity": 4.0 * math.pi * area / perim ** 2 if perim > 0 else float("nan"),
    }

# slope
def slope_by_elevation(dem, slope, base, mask, nbins=20):
    """Mean and median slope in normalised-height bands (MORVOLC-style profile)."""
    h = np.where(mask, dem - base, np.nan)
    hmax = np.nanmax(h)
    if not np.isfinite(hmax) or hmax <= 0:
        return None
    frac = h / hmax
    edges = np.linspace(0, 1, nbins + 1)
    rows = []
    for i in range(nbins):
        sel = mask & (frac >= edges[i]) & (frac < edges[i + 1]) & np.isfinite(slope)
        if sel.sum() < 5:
            continue
        rows.append({
            "h_norm_mid": float((edges[i] + edges[i + 1]) / 2),
            "n_cells": int(sel.sum()),
            "slope_mean_deg": float(np.mean(slope[sel])),
            "slope_median_deg": float(np.median(slope[sel])),
        })
    return rows

