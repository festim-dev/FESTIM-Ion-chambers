"""Digitise the experimental points from the figures in the TFG 2026 slides.

Source: N. J. Sovare, W. T. Shmayda, "Modelling the Performance of Process Monitors
from First Principles", 50th Tritium Focus Group, 19 Oct 2026. The figures in
figures/ are the raster images embedded in that deck, unmodified.

Only the experimental markers are extracted (no model or theory series). Markers are
found by colour. Circular markers are located by fitting a disc of the known marker
radius, which tolerates partial occlusion by lines, other markers and the axes;
other markers use the blob centroid. Every CSV row carries a quality flag:

    ok        marker fully or nearly fully visible
    occluded  less than 60 % of the marker visible; position less certain
    split     marker separated from touching same-colour markers
    clipped   marker cut by the plot frame

Run from anywhere:  python data/experimental/digitise_slides.py
Writes the CSVs next to this file and check overlays to figures/checks/.
"""

import csv
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi
from scipy.cluster.vq import kmeans2

HERE = Path(__file__).parent
FIG = HERE / "figures"
CHECKS = FIG / "checks"


# --------------------------------------------------------------------------- tools
def axis_map(p0, v0, p1, v1, log=False):
    """Pixel -> data map from two reference ticks."""
    if log:
        v0, v1 = np.log10(v0), np.log10(v1)
    k = (v1 - v0) / (p1 - p0)
    return lambda p: 10 ** (v0 + k * (p - p0)) if log else v0 + k * (p - p0)


def load(name):
    return np.asarray(Image.open(FIG / name).convert("RGB")).astype(int)


def colour_mask(im, rgb, tol, box, exclude=()):
    m = np.abs(im - np.array(rgb)).max(2) <= tol
    keep = np.zeros_like(m)
    x0, x1, y0, y1 = box
    keep[y0:y1, x0:x1] = True
    for ex0, ex1, ey0, ey1 in exclude:
        keep[ey0:ey1, ex0:ex1] = False
    return m & keep, keep


def blob_markers(im, rgb, tol, area, box, exclude=(), close=2, min_frac=0.2):
    """Non-circular markers (triangles, small dots): blob centroids.

    `close` bridges the gap left where a line is drawn over a marker. Blobs much
    larger than one marker are split by k-means with k = area / nominal area.
    """
    m, keep = colour_mask(im, rgb, tol, box, exclude)
    if close:
        m = ndi.binary_closing(m, np.ones((2 * close + 1,) * 2), border_value=0) & keep
    lab, n = ndi.label(m)
    out = []
    for i in range(1, n + 1):
        yy, xx = np.nonzero(lab == i)
        if len(xx) < min_frac * area:
            continue
        k = max(1, int(round(len(xx) / area)))
        if k == 1:
            q = "ok" if len(xx) >= 0.6 * area else "occluded"
            out.append((xx.mean(), yy.mean(), q))
        else:
            c, _ = kmeans2(np.c_[xx, yy].astype(float), k, seed=1, minit="++")
            out += [(cx, cy, "split") for cx, cy in c]
    return out


def fit_disc(xx, yy, white, r, cx0, cy0):
    """Centre of a radius-r disc maximising (blob px inside) - (background px inside).

    Pixels hidden under other markers or lines count neither way, so occlusion and
    clipping by the frame do not bias the centre.
    """
    H, W = white.shape
    oy, ox = np.mgrid[-int(r) - 1:int(r) + 2, -int(r) - 1:int(r) + 2]
    d = ox**2 + oy**2 <= r**2
    ox, oy = ox[d], oy[d]
    s_r = int(np.ceil(r))
    best = (-np.inf, cx0, cy0)
    for dy in np.arange(-s_r, s_r + 0.01, 0.5):
        for dx in np.arange(-s_r, s_r + 0.01, 0.5):
            cx, cy = cx0 + dx, cy0 + dy
            px, py = np.round(cx + ox).astype(int), np.round(cy + oy).astype(int)
            ok = (px >= 0) & (px < W) & (py >= 0) & (py < H)
            s = ((xx - cx) ** 2 + (yy - cy) ** 2 <= r**2).sum() - white[py[ok], px[ok]].sum()
            if s > best[0]:
                best = (s, cx, cy)
    return best[1], best[2]


def disc_markers(im, rgb, tol, r, box, exclude=(), close=0, open_=0, min_area=10):
    """Circular markers of known radius r (pixels)."""
    m, keep = colour_mask(im, rgb, tol, box, exclude)
    if open_:  # strip thin lines drawn in the marker colour
        yy, xx = np.mgrid[-open_:open_ + 1, -open_:open_ + 1]
        m = ndi.binary_opening(m, structure=xx**2 + yy**2 <= open_**2)
    if close:
        m = ndi.binary_closing(m, np.ones((2 * close + 1,) * 2), border_value=0) & keep
    white = im.min(2) >= 235
    full = np.pi * r**2
    x0, x1, y0, y1 = box
    lab, n = ndi.label(m)
    out = []
    for i in range(1, n + 1):
        yy, xx = np.nonzero(lab == i)
        if len(xx) < min_area:
            continue
        if xx.max() - xx.min() + 1 <= 2.6 * r:
            cx, cy = fit_disc(xx, yy, white, r, xx.mean(), yy.mean())
            out.append([cx, cy, len(xx), "ok" if len(xx) >= 0.6 * full else "occluded"])
            continue
        found, a0 = [], len(xx)
        while len(xx) >= min_area:  # overlapping markers in a row: greedy, left first
            sel = xx <= xx.min() + 2 * r
            cx, cy = fit_disc(xx, yy, white, r, xx.min() + r, yy[sel].mean())
            found.append([cx, cy, len(xx), "split"])
            rest = (xx - cx) ** 2 + (yy - cy) ** 2 > (r + 1.5) ** 2
            xx, yy = xx[rest], yy[rest]
        if len(found) == 1:  # wide only because of a stray fragment
            found[0][3] = "ok" if a0 >= 0.6 * full else "occluded"
        out += found
    out.sort(key=lambda p: -p[2])  # drop leftover fragments duplicating a marker
    kept = []
    for p in out:
        if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 > r**2 for q in kept):
            if p[1] + r > y1 or p[0] - r < x0:
                p[3] = "clipped"
            kept.append(p)
    return [(cx, cy, q) for cx, cy, _, q in kept]


def overlay(name, pts, r=4):
    img = Image.open(FIG / name).convert("RGB")
    d = ImageDraw.Draw(img)
    for cx, cy, q in pts:
        col = (255, 0, 255) if q == "ok" else (255, 140, 0)
        d.ellipse([cx - r, cy - r, cx + r, cy + r], outline=col)
        d.point((cx, cy), fill=(0, 0, 0))
    CHECKS.mkdir(exist_ok=True)
    img.resize((img.width * 2, img.height * 2), Image.NEAREST).save(CHECKS / name)


def write(fname, header, rows):
    with open(HERE / fname, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in rows:
            w.writerow([f"{v:.4g}" if isinstance(v, float) else v for v in row])
    print(f"{fname}: {len(rows)} points")


# ------------------------------------------------------------------------ figures
def calibration_1L_N2(name, out, yname):
    """Slides 12 and 20: 1 L monitor, N2 carrier, 0-50 uCi/m3, linear axes."""
    im = load(name)
    fx, fy = axis_map(98, 0, 553, 50), axis_map(414, 0, 68, 50)
    box, legend = (82, 575, 60, 426), [(84, 315, 62, 178)]
    series = [  # (flow LPM, bias V, colour, nominal marker area px)
        (5.0, -100, (255, 0, 0), 21), (5.0, -25, (0, 128, 0), 49),
        (7.5, -25, (0, 0, 255), 49), (7.5, +25, (0, 0, 0), 21)]
    rows, allpts = [], []
    for flow, bias, rgb, area in series:
        pts = blob_markers(im, rgb, 60 if rgb != (0, 128, 0) else 50, area, box, legend)
        allpts += pts
        rows += [(flow, bias, float(fx(cx)), float(fy(cy)), q) for cx, cy, q in pts]
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    overlay(name, allpts)
    write(out, ["flow_LPM", "bias_V", "concentration_uCi_m3", yname, "quality"], rows)


def voltage_dependence():
    """Slide 12: detected activity vs bias magnitude at three concentrations, log-log."""
    name = "slide12_voltage_dependence.png"
    im = load(name)
    fx, fy = axis_map(155, 1, 656, 100, log=True), axis_map(518, 1e-2, 156, 1, log=True)
    pts = disc_markers(im, (255, 0, 0), 60, 6.8, (121, 864, 88, 641), close=2, min_area=30)
    overlay(name, pts, r=6.8)
    nominal = np.array([1, 3, 10, 30, 60, 100, 125, 150, 175, 200, 300, 400, 500])
    labels = ["9.7 mCi/m3", "105 mCi/m3", "1.74 Ci/m3"]  # from the curve labels
    rows = []
    V = np.array([float(fx(cx)) for cx, _, _ in pts])
    Vn = nominal[np.abs(np.log(V[:, None] / nominal[None, :])).argmin(1)]
    for v in nominal:  # three points per voltage, lowest = 9.7 mCi/m3
        grp = sorted([(float(fy(cy)), V[i], q) for i, (cx, cy, q) in enumerate(pts) if Vn[i] == v])
        assert len(grp) == 3, (v, grp)
        rows += [(lab, int(v), Vd, y, q) for lab, (y, Vd, q) in zip(labels, grp)]
    rows.sort(key=lambda r: (labels.index(r[0]), r[1]))
    write("voltage_dependence_slide12.csv",
          ["series", "bias_magnitude_V_nominal", "bias_magnitude_V_digitised",
           "detected_activity_Ci_m3", "quality"], rows)


def pressure_dependence_Ar():
    """Slide 12: detected activity vs argon pressure, 9.45 mCi/m3."""
    name = "slide12_Ar_pressure_dependence.png"
    im = load(name)
    fx, fy = axis_map(102, 0, 553, 800), axis_map(427, 1e-3, 70, 1e-2, log=True)
    # the model line has the marker colour: an opening removes it first
    pts = disc_markers(im, (31, 119, 180), 40, 4.2, (81, 576, 59, 427),
                       exclude=[(385, 572, 370, 425)], open_=2)
    overlay(name, pts, r=4.2)
    rows = sorted((float(fx(cx)), float(fy(cy)), q) for cx, cy, q in pts)
    write("Ar_pressure_dependence_slide12.csv",
          ["pressure_Torr", "detected_activity_Ci_m3", "quality"], rows)


def wire_cage_pressures():
    """Slide 14: 20 cc wire cage, 100-1000 Torr, current vs concentration, log-log.

    The ten series overlap almost completely and 1000 Torr is drawn on top. Only
    markers with at least 20 visible pixels are kept: all of 1000 Torr, 100 Torr where
    it shows below the stack, and 200/300 Torr at the lowest concentration.
    """
    name = "slide14_20cc_pressure_series.png"
    im = load(name)
    fx, fy = axis_map(348.5, 100, 814.5, 1e5, log=True), axis_map(630.5, 1e-7, 172.5, 1e-2, log=True)
    colours = {100: (66, 133, 244), 200: (234, 67, 53), 300: (251, 188, 4),
               400: (52, 168, 83), 500: (255, 109, 1), 600: (70, 189, 198),
               700: (123, 170, 247), 800: (240, 123, 114), 900: (252, 208, 79),
               1000: (113, 194, 135)}
    rows, allpts = [], []
    for p, rgb in colours.items():
        pts = disc_markers(im, rgb, 24, 8.6, (212, 830, 150, 708), min_area=20)
        allpts += pts
        rows += [(p, float(fx(cx)), float(fy(cy)), q) for cx, cy, q in pts]
    overlay(name, allpts, r=8.6)
    rows.sort(key=lambda r: (r[0], r[1]))
    write("wire_cage_20cc_pressures_slide14.csv",
          ["pressure_Torr", "concentration_Ci_m3", "current_A", "quality"], rows)


def wire_cage_Ar():
    """Slide 21: 20 cc wire cage in argon, current vs concentration, log-log."""
    name = "slide21_20cc_Ar.png"
    im = load(name)
    fx = axis_map(446.5, 100, 1222.5, 1e5, log=True)
    fy = axis_map(644.5, 1e-9, 132.5, 1e-4, log=True)
    pts = disc_markers(im, (66, 133, 244), 20, 8.6, (226, 1240, 40, 708), min_area=8)
    overlay(name, pts, r=8.6)
    rows = sorted((float(fx(cx)), float(fy(cy)), q) for cx, cy, q in pts)
    write("wire_cage_20cc_Ar_slide21.csv", ["concentration_Ci_m3", "current_A", "quality"], rows)


if __name__ == "__main__":
    calibration_1L_N2("slide12_1L_N2_calibration.png", "calibration_1L_N2_slide12.csv",
                      "current_fA")
    calibration_1L_N2("slide20_1L_N2_calibration_model.png", "calibration_1L_N2_slide20.csv",
                      "indicated_activity_uCi_m3")
    voltage_dependence()
    pressure_dependence_Ar()
    wire_cage_pressures()
    wire_cage_Ar()
