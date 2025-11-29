#!/usr/bin/env python3
"""
Auto-generate rectangular ROI polygons by clustering FOV positions.

Intended for CosMx datasets where multiple tissue slices/blocks were stitched
into one dataset, forming several spatial "islands". This script clusters FOVs
into K groups and writes K ROI CSV files that are directly consumable by
geneSCOPE::createSCOPE_cosmx (expects a file with columns 'x','y').

Coordinate source priority (no heavy dependencies like pyarrow):
  1) flatFiles/<dataset>/*_fov_positions_file.csv.gz (uses x_global_mm/y_global_mm if present)
  2) If 1) not found, you can pass a custom CSV via --fov-csv with columns containing
     x/y in either µm (recommended) or pixels (then provide --pixel-size-um).

Clustering methods:
  - kmeans (default): robust when islands are well separated; requires a target K
  - graph: build a radius graph using eps = 1.5 * median nearest-neighbor distance and
           take connected components (DBSCAN-like, no sklearn). Use --method graph.

Outputs:
  - <out_dir>/roi_cluster_01.csv ... K files. Each contains a closed rectangle polygon
    (five rows: four corners plus a repeated first point) with columns 'x','y' in µm.
  - <out_dir>/roi_manifest.csv mapping of FOV -> cluster_id and bounding boxes.

Example:
  python genescope/scripts/auto_rois_from_fovs.py \
    --root _External__esophageal2__Resegmentation__13_08_2025_22_44_20_700 \
    --k 6 --out-dir AutoROIs --method kmeans --margin-um 50

Processing logic (overview):
  1) Parse CLI arguments: dataset root, clustering method/params, ROI padding, optional manifest/preview paths.
  2) Load FOV anchor positions (µm): prefer an explicit --fov-csv; otherwise read flatFiles/<dataset>/*_fov_positions_file.csv.gz
     using x_global_mm/y_global_mm when present or x_global_px/y_global_px with --pixel-size-um.
  3) Cluster FOV coordinates:
       - kmeans (default) with lightweight numpy implementation, requires target K (>0).
       - graph mode builds a radius graph (eps = 1.5× median nearest-neighbor distance unless overridden) and takes connected components;
         optional fallback to kmeans if component count ≠ K unless --no-fallback.
       - Tiny islands with < --min-fovs can be reassigned to nearest kept component.
  4) Optionally enforce a minimum aspect ratio per ROI via axis binning and splitting (--min-xy-fov-ratio).
  5) For each cluster (or split part), compute bounding box, apply margin (--margin-um), and write a closed rectangle polygon CSV
     (x,y in µm) plus an ROI manifest mapping FOVs to cluster IDs and box extents.
  6) Unless disabled with --no-preview, render a PNG preview showing FOV points and ROI rectangles (saved to --preview-png or <out-dir>/roi_preview.png).
"""

from __future__ import annotations

import argparse
import csv
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import gzip
import math
import numpy as np
import pandas as pd


@dataclass
class FOVPoint:
    fov_id: str
    x_um: float
    y_um: float


def _load_fov_positions_auto(root: Path, pixel_size_um: float | None, explicit_csv: Path | None) -> List[FOVPoint]:
    """Load FOV anchors as (x_um, y_um) from known sources without pyarrow.

    Priority:
      - if explicit_csv provided, use it;
      - else look for flatFiles/<dataset>/<dataset>_fov_positions_file.csv.gz
    """
    # 1) explicit CSV
    if explicit_csv is not None:
        return _load_fov_positions_from_csv(explicit_csv, pixel_size_um)

    # 2) CosMx flatFiles fov positions
    flat = root / "flatFiles"
    if flat.exists():
        # pick the first dataset directory
        ds = next((p for p in flat.iterdir() if p.is_dir()), None)
        if ds is not None:
            cands = sorted(ds.glob(f"{ds.name}_fov_positions_file.csv.gz"))
            if cands:
                return _load_fov_positions_from_csv(cands[0], pixel_size_um)

    raise FileNotFoundError(
        "Could not locate FOV positions. Provide --fov-csv or ensure flatFiles/<dataset>/*_fov_positions_file.csv.gz exists."
    )


def _load_fov_positions_from_csv(path: Path, pixel_size_um: float | None) -> List[FOVPoint]:
    # Support .gz and plain CSV; try pandas for convenience
    is_gz = path.suffix == ".gz"
    open_fn = gzip.open if is_gz else open
    df = pd.read_csv(open_fn(path, "rt"))
    cols = {c.lower(): c for c in df.columns}

    # FOV id
    fov_col = None
    for k in ("fov", "fov_id", "fovid", "fovindex"):
        if k in cols:
            fov_col = cols[k]
            break
    if fov_col is None:
        # fallback: try to synthesize from row index
        fov_ids = [f"FOV{int(i+1):05d}" for i in range(len(df))]
    else:
        raw = df[fov_col]
        fov_ids = [f"FOV{int(v):05d}" if _is_number(v) else str(v) for v in raw]

    # Prefer x/y in mm (already global stage frame), else use px with pixel_size
    if "x_global_mm" in cols and "y_global_mm" in cols:
        x = pd.to_numeric(df[cols["x_global_mm"]], errors="coerce").to_numpy(dtype=float) * 1000.0
        y = pd.to_numeric(df[cols["y_global_mm"]], errors="coerce").to_numpy(dtype=float) * 1000.0
    elif "x_global_px" in cols and "y_global_px" in cols and pixel_size_um is not None:
        x = pd.to_numeric(df[cols["x_global_px"]], errors="coerce").to_numpy(dtype=float) * float(pixel_size_um)
        y = pd.to_numeric(df[cols["y_global_px"]], errors="coerce").to_numpy(dtype=float) * float(pixel_size_um)
    else:
        raise ValueError(
            "CSV does not contain x_global_mm/y_global_mm; and x_global_px/y_global_px present but --pixel-size-um not provided."
        )

    out: List[FOVPoint] = []
    for fid, xi, yi in zip(fov_ids, x, y):
        if np.isfinite(xi) and np.isfinite(yi):
            out.append(FOVPoint(fid, float(xi), float(yi)))
    if not out:
        raise RuntimeError("No valid FOV positions parsed from CSV.")
    return out


def _is_number(x) -> bool:
    try:
        float(x)
        return True
    except Exception:
        return False


def _kmeans(points: np.ndarray, k: int, n_init: int = 10, max_iter: int = 200, rng: np.random.Generator | None = None) -> np.ndarray:
    """Lightweight K-Means (numpy only). Returns labels for each row in points."""
    if rng is None:
        rng = np.random.default_rng()
    best_inertia = math.inf
    best_labels = None
    for _ in range(n_init):
        # k-means++ seeding (simplified)
        centers = _kmeanspp(points, k, rng)
        labels = None
        for _ in range(max_iter):
            # assign
            d2 = _cdist2(points, centers)
            new_labels = np.argmin(d2, axis=1)
            if labels is not None and np.array_equal(new_labels, labels):
                break
            labels = new_labels
            # update
            new_centers = np.vstack([points[labels == j].mean(axis=0) if np.any(labels == j) else centers[j] for j in range(k)])
            centers = new_centers
        inertia = float(np.sum(np.min(_cdist2(points, centers), axis=1)))
        if inertia < best_inertia:
            best_inertia = inertia
            best_labels = labels
    assert best_labels is not None
    return best_labels


def _kmeanspp(points: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    n = points.shape[0]
    centers = np.empty((k, points.shape[1]), dtype=float)
    # first center random
    idx = rng.integers(0, n)
    centers[0] = points[idx]
    # subsequent centers by D^2 weighting
    closest_d2 = _cdist2(points, centers[0:1]).ravel()
    for i in range(1, k):
        probs = closest_d2 / closest_d2.sum()
        idx = rng.choice(n, p=probs)
        centers[i] = points[idx]
        closest_d2 = np.minimum(closest_d2, _row_d2(points, centers[i]))
    return centers


def _cdist2(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    # squared Euclidean distance
    # (a^2 + b^2 - 2ab)
    a2 = np.sum(a * a, axis=1)[:, None]
    b2 = np.sum(b * b, axis=1)[None, :]
    ab = a @ b.T
    return a2 + b2 - 2.0 * ab


def _row_d2(a: np.ndarray, row: np.ndarray) -> np.ndarray:
    diff = a - row[None, :]
    return np.sum(diff * diff, axis=1)


def _graph_components(points: np.ndarray, eps: float) -> np.ndarray:
    """Connected components using radius graph with threshold eps (µm). Returns labels [0..C-1]."""
    n = points.shape[0]
    labels = -np.ones(n, dtype=int)
    adj = _radius_neighbors(points, eps)
    cid = 0
    for i in range(n):
        if labels[i] != -1:
            continue
        # BFS
        queue = [i]
        labels[i] = cid
        while queue:
            u = queue.pop()
            for v in adj[u]:
                if labels[v] == -1:
                    labels[v] = cid
                    queue.append(v)
        cid += 1
    return labels


def _radius_neighbors(points: np.ndarray, eps: float) -> List[List[int]]:
    n = points.shape[0]
    adj: List[List[int]] = [[] for _ in range(n)]
    # naive O(n^2) is fine for n ~ few hundreds
    for i in range(n):
        # vectorized distances from i to all j>i
        d2 = _row_d2(points, points[i])
        within = np.where((d2 <= eps * eps) & (np.arange(n) != i))[0]
        for j in within:
            adj[i].append(j)
            adj[j].append(i)
    return adj


def _estimate_axis_step(values: np.ndarray) -> float:
    if values.size < 2:
        return 0.0
    sorted_vals = np.sort(values)
    diffs = np.diff(sorted_vals)
    diffs = diffs[diffs > 1e-3]
    if diffs.size == 0:
        return 0.0
    hi = float(np.percentile(diffs, 75))
    trimmed = diffs[diffs <= hi]
    ref = trimmed if trimmed.size > 0 else diffs
    return float(np.median(ref))


def _axis_bins(points: np.ndarray, idx: np.ndarray, axis: int, axis_step: float) -> List[np.ndarray]:
    """Group FOV indices by axis-aligned bins using a dataset-wide step estimate."""
    if idx.size == 0:
        return []
    coords = points[idx, axis]
    if coords.size == 0:
        return []
    order = np.argsort(coords)
    sorted_idx = idx[order]
    sorted_vals = coords[order]
    tol = max(5.0, axis_step * 0.5) if axis_step > 0 else 5.0
    bins: List[np.ndarray] = []
    start = 0
    last = sorted_vals[0]
    for i in range(1, sorted_vals.size):
        if sorted_vals[i] - last > tol:
            bins.append(sorted_idx[start:i])
            start = i
        last = sorted_vals[i]
    bins.append(sorted_idx[start:])
    return bins


def _enforce_min_ratio(points: np.ndarray, idx: np.ndarray, min_ratio: float, step_x: float, step_y: float) -> List[np.ndarray]:
    """Split ROI indices so that min(nx, ny)/max(nx, ny) >= min_ratio."""
    if min_ratio <= 0 or idx.size == 0:
        return [idx]
    pending: deque[np.ndarray] = deque([idx])
    out: List[np.ndarray] = []
    while pending:
        cur = pending.popleft()
        bins_x = _axis_bins(points, cur, axis=0, axis_step=step_x)
        bins_y = _axis_bins(points, cur, axis=1, axis_step=step_y)
        count_x = len(bins_x) if bins_x else (1 if cur.size else 0)
        count_y = len(bins_y) if bins_y else (1 if cur.size else 0)
        short = min(count_x, count_y)
        long = max(count_x, count_y)
        if short == 0 or long == 0 or long == 1 or short / long >= min_ratio:
            out.append(cur)
            continue
        if count_x >= count_y:
            axis_bins = bins_x
            other_count = count_y
            axis_name = "x"
        else:
            axis_bins = bins_y
            other_count = count_x
            axis_name = "y"
        max_bins_allowed = int(math.floor(other_count / min_ratio))
        if max_bins_allowed < 1:
            max_bins_allowed = 1
        if max_bins_allowed >= len(axis_bins):
            print(
                f"[warn] Could not enforce min-ratio {min_ratio} for ROI (axis={axis_name}, bins={len(axis_bins)}, other={other_count})."
            )
            out.append(cur)
            continue
        start = 0
        while start < len(axis_bins):
            group = axis_bins[start : start + max_bins_allowed]
            new_idx = np.concatenate(group).astype(int, copy=False)
            pending.append(new_idx)
            start += max_bins_allowed
    return out


def _write_roi_csv(path: Path, xmin: float, xmax: float, ymin: float, ymax: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x", "y"])
        # closed rectangle polygon (repeat first point at end)
        w.writerow([xmin, ymin])
        w.writerow([xmax, ymin])
        w.writerow([xmax, ymax])
        w.writerow([xmin, ymax])
        w.writerow([xmin, ymin])


def main() -> None:
    ap = argparse.ArgumentParser(description="Cluster FOV positions into spatial islands and write ROI rectangles (µm)")
    ap.add_argument("--root", required=True, type=Path, help="CosMx dataset root (contains flatFiles/ and Parquets)")
    ap.add_argument("--k", type=int, default=6, help="Number of ROI clusters; set 0 to 'auto' (one ROI per connected island in graph mode)")
    ap.add_argument("--method", choices=["kmeans", "graph"], default="kmeans", help="Clustering method (graph recommended for stitched islands)")
    ap.add_argument("--fov-csv", type=Path, default=None, help="Optional path to FOV positions CSV(.gz) if not using flatFiles default")
    ap.add_argument("--pixel-size-um", type=float, default=None, help="Only needed when CSV lacks mm columns and has pixels")
    ap.add_argument("--margin-um", type=float, default=0.0, help="Pad each ROI bbox by this many micrometers")
    ap.add_argument("--out-dir", type=Path, default=Path("AutoROIs"), help="Output directory for ROI CSV files")
    ap.add_argument("--manifest", type=Path, default=None, help="Optional path for a FOV->cluster manifest CSV")
    ap.add_argument("--preview-png", type=Path, default=None, help="Optional path to save a preview PNG showing FOVs and ROI boxes; default is <out-dir>/roi_preview.png")
    ap.add_argument("--no-preview", action="store_true", help="Disable preview image generation")
    ap.add_argument("--min-fovs", type=int, default=2, help="Ignore tiny islands with < this many FOVs (graph mode)")
    ap.add_argument(
        "--min-xy-fov-ratio",
        type=float,
        default=0.0,
        help=(
            "Shape guardrail per ROI. Values <=1 enforce min(nx,ny)/max(nx,ny) >= value; "
            ">1 are interpreted as max(nx,ny)/min(nx,ny) <= value. Set 0 to disable."
        ),
    )
    ap.add_argument("--eps-um", type=float, default=None, help="Override graph radius (µm) used to connect neighboring FOVs; larger merges islands")
    ap.add_argument("--no-fallback", action="store_true", help="When --method graph and components != K, do NOT fallback to kmeans (keep components)")
    args = ap.parse_args()

    if args.min_xy_fov_ratio < 0:
        ap.error("--min-xy-fov-ratio must be >= 0.")

    ratio_arg = float(args.min_xy_fov_ratio)
    if ratio_arg == 0:
        ratio_threshold = 0.0
        ratio_desc = "disabled"
    elif ratio_arg <= 1.0:
        ratio_threshold = ratio_arg
        ratio_desc = f"min/max >= {ratio_arg}"
    else:
        ratio_threshold = 1.0 / ratio_arg
        ratio_desc = f"max/min <= {ratio_arg}"
        print(
            f"[info] --min-xy-fov-ratio={ratio_arg} interpreted as max(nx,ny)/min(nx,ny) <= {ratio_arg} "
            f"(min/max >= {ratio_threshold:.4f})."
        )

    root = args.root.resolve()
    pts = _load_fov_positions_auto(root, args.pixel_size_um, args.fov_csv)
    arr = np.array([(p.x_um, p.y_um) for p in pts], dtype=float)
    grid_step_x = _estimate_axis_step(arr[:, 0])
    grid_step_y = _estimate_axis_step(arr[:, 1])

    if args.method == "kmeans":
        k = max(1, int(args.k))
        labels = _kmeans(arr, k)
    else:
        # Graph islands: eps = 1.5 * median nearest-neighbor distance
        if arr.shape[0] < 2:
            labels = np.zeros(arr.shape[0], dtype=int)
        else:
            if args.eps_um is not None and args.eps_um > 0:
                eps = float(args.eps_um)
            else:
                # nearest neighbor distance for each point
                d2 = _cdist2(arr, arr)
                np.fill_diagonal(d2, np.inf)
                nn = np.sqrt(np.min(d2, axis=1))
                eps = 1.5 * float(np.median(nn))
            labels = _graph_components(arr, eps)
        comps = sorted(set(labels.tolist()))
        # Optionally drop tiny components
        if args.min_fovs > 1:
            sizes = {c: int(np.sum(labels == c)) for c in comps}
            keep = [c for c in comps if sizes[c] >= args.min_fovs]
            if keep and len(keep) != len(comps):
                # reassign dropped points to nearest kept centroid
                kept_centers = np.vstack([arr[labels == c].mean(axis=0) for c in keep])
                kept_ids = np.array(keep)
                for c in comps:
                    if c in keep:
                        continue
                    idx = (labels == c)
                    if not np.any(idx):
                        continue
                    d2k = _cdist2(arr[idx], kept_centers)
                    labels[idx] = kept_ids[np.argmin(d2k, axis=1)]
                comps = keep
        # K harmonization
        if int(args.k) == 0:
            # auto: keep graph components as-is
            pass
        else:
            if len(set(labels.tolist())) != int(args.k):
                if args.no_fallback:
                    # keep components and warn via print
                    print(f"[warn] graph found {len(set(labels))} components != K={args.k}; keeping components (no-fallback)")
                else:
                    labels = _kmeans(arr, int(args.k))

    # Build bboxes per cluster
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cluster_ids = sorted(set(labels.tolist()))
    roi_segments: List[tuple[int, np.ndarray]] = []
    split_events: List[tuple[int | str, int]] = []
    for lab in cluster_ids:
        idx = np.where(labels == lab)[0]
        parts = _enforce_min_ratio(arr, idx, ratio_threshold, grid_step_x, grid_step_y)
        if ratio_threshold > 0 and len(parts) > 1:
            split_events.append((lab, len(parts)))
        for seg in parts:
            if seg.size == 0:
                continue
            roi_segments.append((lab, seg))

    if split_events:
        for lab, count in split_events:
            print(
                f"[info] cluster {lab} split into {count} ROI(s) to satisfy ratio constraint ({ratio_desc})."
            )

    # optional manifest
    man_rows: List[dict] = []
    roi_boxes = {}

    for cid, (lab, idx) in enumerate(roi_segments, start=1):
        xs = arr[idx, 0]
        ys = arr[idx, 1]
        xmin = float(xs.min()); xmax = float(xs.max())
        ymin = float(ys.min()); ymax = float(ys.max())
        if args.margin_um > 0:
            xmin -= args.margin_um; ymin -= args.margin_um
            xmax += args.margin_um; ymax += args.margin_um
        out_csv = out_dir / f"roi_cluster_{cid:02d}.csv"
        _write_roi_csv(out_csv, xmin, xmax, ymin, ymax)
        roi_boxes[cid] = (xmin, xmax, ymin, ymax)

        for i in idx:
            man_rows.append({
                "fov_id": pts[i].fov_id,
                "cluster_id": cid,
                "x_um": float(arr[i, 0]),
                "y_um": float(arr[i, 1]),
                "roi_xmin": xmin,
                "roi_xmax": xmax,
                "roi_ymin": ymin,
                "roi_ymax": ymax,
            })

    manifest_path = args.manifest if args.manifest is not None else (out_dir / "roi_manifest.csv")
    pd.DataFrame(man_rows).to_csv(manifest_path, index=False)
    print(f"[OK] Wrote {len(roi_segments)} ROI CSVs to {out_dir}")
    print(f"[OK] Manifest: {manifest_path}")

    # Preview figure (optional)
    if not args.no_preview:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import Rectangle

            # Downsample points for readability if needed
            n = arr.shape[0]
            step = max(1, n // 8000)
            arr_ds = arr[::step]
            labels_ds = labels[::step]

            fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
            sc = ax.scatter(arr_ds[:, 0], arr_ds[:, 1], c=labels_ds, s=6, alpha=0.6, cmap="tab20")

            # draw ROI rectangles
            for cid, (xmin, xmax, ymin, ymax) in roi_boxes.items():
                w = xmax - xmin
                h = ymax - ymin
                rect = Rectangle((xmin, ymin), w, h, fill=False, lw=1.5, ec="k")
                ax.add_patch(rect)
                ax.text(xmin + 4, ymax - 4, f"ROI {cid:02d}", fontsize=8, color="k", ha="left", va="top",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))

            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("X (µm)")
            ax.set_ylabel("Y (µm)")
            ax.set_title("ROI islands preview (FOV positions + ROI boxes)")
            ax.grid(alpha=0.2, lw=0.3)
            preview_path = args.preview_png if args.preview_png is not None else (out_dir / "roi_preview.png")
            fig.tight_layout()
            fig.savefig(preview_path)
            plt.close(fig)
            print(f"[OK] Preview saved to {preview_path}")
        except Exception as e:
            print(f"[warn] Could not generate preview image: {e}")


if __name__ == "__main__":
    main()
auto_rois_from_fovs.py
