"""Writing the outline to GeoJSON and drawing the overview figure."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage import measure

#export
def mask_to_geojson(mask, cell_size, origin, path):
    """Largest contour of the mask as a GeoJSON polygon in map coordinates."""
    padded = np.pad(mask.astype(float), 1)
    contours = measure.find_contours(padded, 0.5)
    if not contours:
        return None
    c = max(contours, key=len) - 1
    xs = origin[0] + c[:, 1] * cell_size
    ys = origin[1] - c[:, 0] * cell_size
    ring = [[float(x), float(y)] for x, y in zip(xs, ys)]
    if ring[0] != ring[-1]:
        ring.append(ring[0])
    gj = {"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"name": "cone_basal_outline"},
         "geometry": {"type": "Polygon", "coordinates": [ring]}}]}
    with open(path, "w") as f:
        json.dump(gj, f)
    return path

def overview_figure(dem, slope, mask, crater, base, stats, bands, cell_size, path, threshold):
    fig, ax = plt.subplots(2, 2, figsize=(13, 11))

    ls = matplotlib.colors.LightSource(azdeg=315, altdeg=45)
    shade = ls.hillshade(np.nan_to_num(dem, nan=float(np.nanmin(dem))),
                         vert_exag=1.5, dx=cell_size, dy=cell_size)
    ax[0, 0].imshow(shade, cmap="gray")
    ax[0, 0].contour(mask.astype(float), levels=[0.5], colors="#1D9E75", linewidths=2)
    if crater is not None:
        ax[0, 0].contour(crater.astype(float), levels=[0.5], colors="#D85A30",
                         linewidths=1.5, linestyles="--")
    ax[0, 0].set_title("Hillshade — basal outline (solid), crater (dashed)")
    ax[0, 0].axis("off")

    im = ax[0, 1].imshow(slope, cmap="magma", vmin=0, vmax=40)
    ax[0, 1].contour(mask.astype(float), levels=[0.5], colors="#1D9E75", linewidths=2)
    plt.colorbar(im, ax=ax[0, 1], shrink=0.8, label="Slope (deg)")
    ax[0, 1].set_title(f"Slope — threshold {threshold:.0f} deg")
    ax[0, 1].axis("off")

    relief = np.where(mask, dem - base, np.nan)
    im2 = ax[1, 0].imshow(relief, cmap="viridis")
    plt.colorbar(im2, ax=ax[1, 0], shrink=0.8, label="Height above base (m)")
    ax[1, 0].set_title("Relief above fitted base surface")
    ax[1, 0].axis("off")

    if bands:
        hn = [b["h_norm_mid"] for b in bands]
        ax[1, 1].plot([b["slope_mean_deg"] for b in bands], hn, "o-", label="mean")
        ax[1, 1].plot([b["slope_median_deg"] for b in bands], hn, "s--", label="median")
        ax[1, 1].set_xlabel("Slope (deg)")
        ax[1, 1].set_ylabel("Normalised height")
        ax[1, 1].set_title("Slope vs. height")
        ax[1, 1].grid(alpha=0.3)
        ax[1, 1].legend()

    plt.tight_layout()
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)