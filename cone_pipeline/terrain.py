"""Slope, fitted base surface, and volume above it."""
import math
import numpy as np
from scipy import ndimage

#surface derivatives
def slope_degrees(dem, cell_size):
    """Horn (1981) 3x3 slope — the same algorithm ArcGIS Pro's Slope tool uses."""
    kx = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], float) / (8.0 * cell_size)
    ky = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]], float) / (8.0 * cell_size)
    filled = np.where(np.isfinite(dem), dem, np.nanmean(dem))
    dzdx = ndimage.correlate(filled, kx, mode="nearest")
    dzdy = ndimage.correlate(filled, ky, mode="nearest")
    return np.degrees(np.arctan(np.hypot(dzdx, dzdy)))

# base surface and volume
def fit_base_surface(dem, mask, cell_size, ring_cells, order):
    """Least-squares trend surface through the collar of terrain outside the cone.

    This is the 2D generalisation of the profile pipeline's fitted regional line,
    and corresponds to MORVOLC's polynomial base. A tilted or gently warped
    regional surface is removed rather than averaged away.
    """
    grown = ndimage.binary_dilation(mask, np.ones((3, 3), bool), iterations=int(ring_cells))
    ring = grown & ~mask & np.isfinite(dem)
    if ring.sum() < 30:
        raise RuntimeError("collar outside the cone is too small for a base fit — "
                           "the DEM may be clipped too tightly around the edifice")

    rr, cc = np.nonzero(ring)
    x = cc * cell_size
    y = rr * cell_size
    z = dem[ring]

    if order == 1:
        A = np.column_stack([np.ones_like(x), x, y])
    else:
        A = np.column_stack([np.ones_like(x), x, y, x * x, x * y, y * y])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)

    R, C = np.indices(dem.shape)
    X = C * cell_size
    Y = R * cell_size
    if order == 1:
        surf = coef[0] + coef[1] * X + coef[2] * Y
        dip = math.degrees(math.atan(math.hypot(coef[1], coef[2])))
    else:
        surf = (coef[0] + coef[1] * X + coef[2] * Y
                + coef[3] * X * X + coef[4] * X * Y + coef[5] * Y * Y)
        dip = math.degrees(math.atan(math.hypot(coef[1], coef[2])))

    resid = z - surf[ring]
    return surf, dip, float(np.sqrt(np.mean(resid ** 2)))


def volume_above_base(dem, base, mask, cell_size):
    """Integrate elevation above the base surface, cell by cell.

    No symmetry assumed and no frustum. `clipped` counts only material above the
    base; `net` allows below-base cells to subtract, which flags a poor fit if the
    two diverge much.
    """
    diff = np.where(mask & np.isfinite(dem), dem - base, 0.0)
    cell_area = cell_size ** 2
    return float(np.clip(diff, 0, None).sum() * cell_area), float(diff.sum() * cell_area)