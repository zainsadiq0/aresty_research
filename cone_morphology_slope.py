"""
Author: Zain Sadiq
Script made to measure cones based on slope (Runs in ArcPY)
Work under Karen Bemis' research on Volcano Morphology

READ ------ This script is highly accurate and suitable for scoria cone shapes,
however more tests must be run on asymettric and "difficult" cones.


Base surface:
  Rather than averaging the two toe elevations, this fits a straight line to
  the terrain outside the detected toes and measures height above that fitted
  regional surface at the cone centre. Per-flank heights are also reported so
  asymmetry stays visible instead of being averaged away.

Outputs to OUTPUT_DIR:
  profile_chart.png        clean elevation profile, all lines overlaid
  diagnostic_profileN.png  per-line check: elevation + slope with picks marked
  cone_table.png           formatted morphology table
  cone_morphology.csv      same table as CSV
"""

import math
from collections import defaultdict

import arcpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter

# ----------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------
TABLE_NAME = "stackprofile"
OUTPUT_DIR = r"C:\path\to\output"

SLOPE_THRESHOLD = 10.0   # degrees. Working range ~6-12. See notes below.
PERSIST_PTS = 3          # consecutive points below threshold to call it "off the cone"
SMOOTH_WINDOW = 15       # odd. Still needed - slope is a derivative.
RIM_PROMINENCE = 10      # metres, for crater rim detection

USE_FITTED_BASE = True   # fit regional surface from terrain outside the toes
SHOW_DIAGNOSTICS = True
# ----------------------------------------------------------------------
# NOTE ON SLOPE_THRESHOLD
#   Above ~12 deg the threshold starts cutting into the lower flank itself
#   and base diameter collapses. Below ~5 deg it starts chasing noise on the
#   plain. Check the printed "outer terrain slope" diagnostic: the threshold
#   should sit comfortably between that value and the flank slope.
# ----------------------------------------------------------------------


def load_profiles(table_name):
    """Pull LINE_ID / distance / elevation straight from the open project."""
    profiles = defaultdict(list)
    with arcpy.da.SearchCursor(table_name,
                               ["LINE_ID", "FIRST_DIST", "FIRST_Z"]) as cur:
        for line_id, dist, elev in cur:
            if dist is None or elev is None:
                continue
            profiles[line_id].append((dist, elev))
    for pts in profiles.values():
        pts.sort(key=lambda p: p[0])
    return profiles


def smooth(elevs, window):
    n = len(elevs)
    w = min(window, n if n % 2 else n - 1)
    if w < 5:
        return elevs.astype(float)
    if w % 2 == 0:
        w -= 1
    return savgol_filter(elevs, window_length=w, polyorder=3)


def find_rims(elevs_s):
    """Two most prominent local maxima, returned left-to-right."""
    peaks, props = find_peaks(elevs_s, prominence=RIM_PROMINENCE)
    if len(peaks) < 2:
        return None
    top2 = np.sort(peaks[np.argsort(props["prominences"])[-2:]])
    return int(top2[0]), int(top2[1])


def find_toes(slope_deg, rim1, rim2, thr, persist):
    """Walk outward from each rim until slope stays below threshold.

    Returns (left_idx, right_idx, left_ok, right_ok). The *_ok flags are
    False when the walk ran off the end of the profile without ever finding
    terrain below threshold - i.e. the line never reached the plain.
    """
    n = len(slope_deg)

    left_idx, left_ok = 0, False
    for i in range(rim1, persist - 1, -1):
        if np.all(slope_deg[i - persist:i] < thr):
            left_idx, left_ok = i, True
            break

    right_idx, right_ok = n - 1, False
    for i in range(rim2, n - persist):
        if np.all(slope_deg[i:i + persist] < thr):
            right_idx, right_ok = i, True
            break

    return left_idx, right_idx, left_ok, right_ok


def fit_regional_surface(dists, elevs_s, left, right):
    """Least-squares line through the terrain outside the cone.

    In profile space a tilted regional plane appears as a sloping line, so a
    linear fit to the outer segments recovers the surface the cone sits on.
    Returns a callable, or None if there aren't enough outside points.
    """
    outside = np.concatenate([np.arange(0, left), np.arange(right + 1, len(dists))])
    if len(outside) < 6:
        return None
    coef = np.polyfit(dists[outside], elevs_s[outside], 1)
    return np.poly1d(coef)


def analyse(dists, elevs):
    elevs_s = smooth(elevs, SMOOTH_WINDOW)
    slope_deg = np.degrees(np.arctan(np.abs(np.gradient(elevs_s, dists))))

    rims = find_rims(elevs_s)
    if rims is None:
        return None, "fewer than 2 rim peaks - lower RIM_PROMINENCE"
    rim1, rim2 = rims

    left, right, left_ok, right_ok = find_toes(
        slope_deg, rim1, rim2, SLOPE_THRESHOLD, PERSIST_PTS)

    base_d = dists[right] - dists[left]
    crater_d = dists[rim2] - dists[rim1]
    rim_z = (elevs_s[rim1] + elevs_s[rim2]) / 2
    centre = (dists[rim1] + dists[rim2]) / 2

    # Base elevation: fitted regional surface, or fall back to toe average
    surface = fit_regional_surface(dists, elevs_s, left, right) if USE_FITTED_BASE else None
    if surface is not None:
        base_z = float(surface(centre))
        base_method = "fitted"
        regional_dip = math.degrees(math.atan(abs(surface.coefficients[0])))
    else:
        base_z = (elevs_s[left] + elevs_s[right]) / 2
        base_method = "toe-average"
        regional_dip = float("nan")

    height = rim_z - base_z
    h_left = rim_z - elevs_s[left]
    h_right = rim_z - elevs_s[right]

    R, c = base_d / 2, crater_d / 2
    if R <= c or height <= 0:
        return None, "degenerate geometry (crater >= base, or height <= 0)"

    stats = {
        "Base diameter (m)": base_d,
        "Top/crater diameter (m)": crater_d,
        "Height (m)": height,
        "Height, left flank (m)": h_left,
        "Height, right flank (m)": h_right,
        "Steepness (h/(r-c))": height / (R - c),
        "Flatness (c/r)": c / R,
        "Volume (m^3)": (math.pi * height / 3) * (R ** 2 + R * c + c ** 2),
    }
    meta = {
        "left": left, "right": right, "rim1": rim1, "rim2": rim2,
        "smooth": elevs_s, "slope": slope_deg,
        "left_ok": left_ok, "right_ok": right_ok,
        "surface": surface, "base_z": base_z, "base_method": base_method,
        "regional_dip": regional_dip,
        "outer_slope_left": float(np.median(slope_deg[:max(3, left)])) if left > 2 else float("nan"),
        "outer_slope_right": float(np.median(slope_deg[right:])) if right < len(dists) - 3 else float("nan"),
        "flank_slope": float(np.median(slope_deg[slope_deg > SLOPE_THRESHOLD])),
        "asymmetry": abs(h_left - h_right),
    }
    return stats, meta


def check(line_id, meta):
    """Return a list of warning strings for this profile."""
    w = []
    if not meta["left_ok"]:
        w.append("left arm never dropped below the slope threshold - the line "
                 "does not reach flat terrain on that side")
    if not meta["right_ok"]:
        w.append("right arm never dropped below the slope threshold - the line "
                 "does not reach flat terrain on that side")
    for side in ("left", "right"):
        v = meta[f"outer_slope_{side}"]
        if not math.isnan(v) and v > SLOPE_THRESHOLD * 0.8:
            w.append(f"terrain outside the {side} toe averages {v:.1f} deg, close "
                     f"to the {SLOPE_THRESHOLD:.0f} deg threshold - toe position "
                     f"is not well constrained on that side")
    if meta["asymmetry"] > 40:
        w.append(f"flank heights differ by {meta['asymmetry']:.0f} m - the cone "
                 f"sits on a sloping surface; report per-flank values")
    return w


def diagnostic_plot(line_id, dists, elevs, meta, out_path):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    zs = meta["smooth"]
    L, R, r1, r2 = meta["left"], meta["right"], meta["rim1"], meta["rim2"]

    axes[0].plot(dists, elevs, lw=0.8, alpha=0.35, label="raw")
    axes[0].plot(dists, zs, lw=1.5, label="smoothed")
    if meta["surface"] is not None:
        axes[0].plot(dists, meta["surface"](dists), ls=":", lw=1.2,
                     label="fitted regional surface")
    axes[0].plot(dists[[L, R]], zs[[L, R]], "v", ms=12, label="toe")
    axes[0].plot(dists[[r1, r2]], zs[[r1, r2]], "^", ms=12, label="rim")
    axes[0].set_ylabel("Elevation (m)")
    axes[0].set_title(f"Profile {line_id} - slope-based detection "
                      f"(threshold {SLOPE_THRESHOLD:.0f} deg)")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].plot(dists, meta["slope"], lw=1)
    axes[1].axhline(SLOPE_THRESHOLD, ls="--", lw=1,
                    label=f"threshold {SLOPE_THRESHOLD:.0f} deg")
    axes[1].axvline(dists[L], ls=":", lw=1)
    axes[1].axvline(dists[R], ls=":", lw=1)
    axes[1].set_xlabel("Distance (m)")
    axes[1].set_ylabel("Slope (deg)")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def presentation_chart(profiles, out_path):
    fig, ax = plt.subplots(figsize=(11, 6))
    for line_id, pts in sorted(profiles.items()):
        ax.plot([p[0] for p in pts], [p[1] for p in pts],
                lw=1.4, label=f"Profile {line_id}")
    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Elevation (m)")
    ax.set_title("Cone elevation profiles")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def format_table(df):
    def fmt(row):
        if "Volume" in row.name:
            return row.map(lambda v: f"{v:.2e}".replace("e+0", "x10^"))
        if "Steepness" in row.name or "Flatness" in row.name:
            return row.map(lambda v: f"{v:.3f}")
        return row.map(lambda v: f"{v:.1f}")
    return df.apply(fmt, axis=1)


def table_image(df_fmt, out_path):
    fig, ax = plt.subplots(figsize=(9.5, 3.0))
    ax.axis("off")
    tbl = ax.table(cellText=df_fmt.values, rowLabels=df_fmt.index,
                   colLabels=df_fmt.columns, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1, 1.6)
    for j in range(len(df_fmt.columns)):
        tbl[(0, j)].set_text_props(weight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def sensitivity(profiles, thresholds=(6, 8, 10, 12)):
    """Re-run toe detection across thresholds so you can report stability."""
    print("\n--- THRESHOLD SENSITIVITY (base diameter, m) ---")
    header = "  thr  " + "  ".join(f"P{ln:>6}" for ln in sorted(profiles))
    print(header)
    for thr in thresholds:
        row = []
        for line_id, pts in sorted(profiles.items()):
            d = np.array([p[0] for p in pts], float)
            z = np.array([p[1] for p in pts], float)
            zs = smooth(z, SMOOTH_WINDOW)
            sl = np.degrees(np.arctan(np.abs(np.gradient(zs, d))))
            rims = find_rims(zs)
            if rims is None:
                row.append(float("nan"))
                continue
            L, R, _, _ = find_toes(sl, rims[0], rims[1], thr, PERSIST_PTS)
            row.append(d[R] - d[L])
        print(f"  {thr:3.0f}  " + "  ".join(f"{v:7.1f}" for v in row))
    print("  (a spread of a few tens of metres is expected; hundreds means the "
          "toe is not well constrained)")


def main():
    profiles = load_profiles(TABLE_NAME)
    if not profiles:
        print("No profile data loaded. Check TABLE_NAME and field names.")
        return

    results, all_warnings = {}, []

    for line_id, pts in sorted(profiles.items()):
        dists = np.array([p[0] for p in pts], dtype=float)
        elevs = np.array([p[1] for p in pts], dtype=float)

        stats, meta = analyse(dists, elevs)
        if stats is None:
            print(f"Profile {line_id}: SKIPPED - {meta}")
            continue

        print(f"Profile {line_id}: {len(dists)} pts, spacing "
              f"~{np.median(np.diff(dists)):.2f} m, length {dists[-1]:.1f} m")
        print(f"  outer terrain {meta['outer_slope_left']:.1f} / "
              f"{meta['outer_slope_right']:.1f} deg | flank median "
              f"{meta['flank_slope']:.1f} deg | regional dip "
              f"{meta['regional_dip']:.1f} deg")
        print(f"  base {stats['Base diameter (m)']:.1f} m, crater "
              f"{stats['Top/crater diameter (m)']:.1f} m, height "
              f"{stats['Height (m)']:.1f} m ({meta['base_method']}) "
              f"[L {stats['Height, left flank (m)']:.1f} / "
              f"R {stats['Height, right flank (m)']:.1f}]")

        for w in check(line_id, meta):
            print(f"  WARNING - {w}")
            all_warnings.append(f"Profile {line_id}: {w}")

        results[line_id] = stats
        if SHOW_DIAGNOSTICS:
            diagnostic_plot(line_id, dists, elevs, meta,
                            rf"{OUTPUT_DIR}\diagnostic_profile{line_id}.png")

    if not results:
        print("\nNo profiles produced measurements.")
        return

    df = pd.DataFrame(results)
    df.columns = [f"Profile {c}" for c in df.columns]
    df["Average"] = df.mean(axis=1)
    df_fmt = format_table(df)
    print("\n" + df_fmt.to_string() + "\n")

    presentation_chart(profiles, rf"{OUTPUT_DIR}\profile_chart.png")
    table_image(df_fmt, rf"{OUTPUT_DIR}\cone_table.png")
    df_fmt.to_csv(rf"{OUTPUT_DIR}\cone_morphology.csv")
    print(f"Saved chart, table, and CSV to {OUTPUT_DIR}")

    sensitivity(profiles)

    if all_warnings:
        print("\n--- WARNINGS ---")
        for w in all_warnings:
            print(w)


if __name__ == "__main__":
    main()
