"""Loading DEMs and converting map coordinates to grid cells."""
import os
import numpy as np

def load_dem(source, cell_size=None):
    """Return (dem, cell_size, origin_xy, nodata_mask)."""
    ext = os.path.splitext(str(source))[1].lower()

    if ext == ".asc":
        return _load_ascii_grid(source)

    if ext == ".npy":
        if cell_size is None:
            raise ValueError("CELL_SIZE must be set for .npy input")
        dem = np.load(source).astype(float)
        return dem, float(cell_size), (0.0, 0.0), np.zeros(dem.shape, bool)

    if ext in (".tif", ".tiff"):
        import rasterio
        with rasterio.open(source) as src:
            dem = src.read(1).astype(float)
            nodata = (dem == src.nodata) if src.nodata is not None else np.zeros(dem.shape, bool)
            return dem, float(src.transform.a), (src.transform.c, src.transform.f), nodata

    # Otherwise assume an ArcGIS Pro layer name
    import arcpy
    r = arcpy.Raster(source)
    dem = arcpy.RasterToNumPyArray(r, nodata_to_value=np.nan).astype(float)
    origin = (r.extent.XMin, r.extent.YMax)
    return dem, float(r.meanCellWidth), origin, ~np.isfinite(dem)


def _load_ascii_grid(path):
    hdr = {}
    with open(path) as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
    dem = np.loadtxt(path, skiprows=6).astype(float)
    nodata = dem == hdr.get("nodata_value", -9999.0)
    dem[nodata] = np.nan
    cs = hdr["cellsize"]
    # xllcorner/yllcorner is the lower-left; convert to upper-left origin
    origin = (hdr.get("xllcorner", 0.0),
              hdr.get("yllcorner", 0.0) + dem.shape[0] * cs)
    return dem, cs, origin, nodata


def rowcol_from_xy(xy, origin, cell_size, shape):
    col = int(round((xy[0] - origin[0]) / cell_size))
    row = int(round((origin[1] - xy[1]) / cell_size))
    row = min(max(row, 0), shape[0] - 1)
    col = min(max(col, 0), shape[1] - 1)
    return row, col
