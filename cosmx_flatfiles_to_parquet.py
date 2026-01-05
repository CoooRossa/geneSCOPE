#!/usr/bin/env python3
"""
Convert CosMx AtomX flatFiles (CSV.GZ) to geneSCOPE-ready Parquet files and
ensure all coordinates are in micrometers (µm), matching geneSCOPE defaults.

Outputs (written under --root):
  - transcripts.parquet      (x_location,y_location in µm; feature_name, qv, nucleus_distance)
  - cells.parquet            (x_centroid,y_centroid in µm; per (fov,cell_id) centroid from polygons)
  - segmentation_boundaries.parquet / cell_boundaries.parquet / nucleus_boundaries.parquet
    (if present, converted in-place to µm with a backup *.px.parquet)

Usage:
  python genescope/scripts/cosmx_flatfiles_to_parquet.py \
      --root <CosMxRunRoot> \
      --pixel-size-um 0.120280945 \
      --build-transcripts --build-cells --convert-segmentation --overwrite

Notes:
  - tx_file without header is supported via positional mapping:
      col6->x_global_px, col7->y_global_px, col9->target, col8->nucleus_distance
  - polygons without header is supported:
      col1->fov, col2->cell_id, col6->x_global_px, col7->y_global_px
  - Existing Parquet files will be overwritten only with --overwrite (backups
    for segmentation are always created as *.px.parquet before conversion).

Processing logic (overview):
  1) Parse CLI flags (root, pixel size, build/convert toggles, segmentation mode)
     and locate the dataset directory under flatFiles/ (first subfolder).
  2) build_transcripts: find *_tx_file.csv(.gz); auto-detect missing header;
     stream in chunks; map x/y from px→µm (or mm→µm fallback); carry nucleus_distance
     either as px converted by --pixel-size-um or as µm when --nucleus-distance-unit=um;
     write Parquet with feature_name, qv (NaN→+inf) and nucleus_distance (NaN→0).
  3) build_cells_from_polygons: read polygon CSV.GZ in chunks; accumulate per (fov,cell)
     centroid in pixel space; scale to µm and write cells.parquet.
  4) build_segmentation_from_polygons: stream polygons and emit segmentation_boundaries.parquet.
     In compat mode (default) it only rescales coordinates; in enhanced mode it also
     breaks long cross-FOV lines (path_id/vertex_order), with optional subsampling and
     per-chunk path joining governed by --seg-* flags.
  5) convert_segmentation_to_um: if existing cell/nucleus/segmentation Parquets already exist,
     back them up as *.px.parquet and rescale coordinate columns to µm, adding label_id when absent.
  6) main() orchestrates optional steps; overwrites outputs only when --overwrite is supplied and
     always leaves a px backup for segmentation conversions.
"""
from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
from typing import Optional, Tuple, Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def find_dataset_dir(flatfiles_dir: Path) -> Path:
    cands = [p for p in flatfiles_dir.iterdir() if p.is_dir()]
    if not cands:
        raise FileNotFoundError(f"No dataset folder under {flatfiles_dir}")
    return sorted(cands)[0]


def is_headerless_csv_gz(path: Path) -> bool:
    """Heuristic: treat as 'headerless' unless first line contains expected header keys.

    CosMx tx_file can have a first data row containing strings (e.g., gene names),
    so presence of letters is not a reliable indicator of a header.
    """
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        first = (f.readline() or "").strip().lower()
    tokens = [t.strip() for t in first.split(",")]
    expected = {"x_global_px", "y_global_px", "target", "x_mm", "y_mm", "x_global_mm", "y_global_mm"}
    return not any(t in expected for t in tokens)


def _parquet_writer(out_path: Path, schema: pa.schema) -> pq.ParquetWriter:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    return pq.ParquetWriter(out_path, schema, compression="zstd")


def build_transcripts(root: Path, dataset_dir: Path, out_path: Path, pixel_um: float, overwrite: bool) -> None:
    if out_path.exists() and not overwrite:
        print(f"[SKIP] {out_path.name} exists (use --overwrite to rebuild)")
        return

    # find *_tx_file.csv(.gz)
    tx_path = None
    for cand in [dataset_dir / f"{dataset_dir.name}_tx_file.csv.gz", dataset_dir / f"{dataset_dir.name}_tx_file.csv"]:
        if cand.exists():
            tx_path = cand
            break
    if tx_path is None:
        raise FileNotFoundError("Could not find *_tx_file.csv(.gz) under flatFiles")

    headerless = tx_path.suffix == ".gz" and is_headerless_csv_gz(tx_path)
    print(f"[INFO] Building {out_path.name} from {tx_path.name}; headerless={headerless}")

    # define schema
    schema = pa.schema([
        pa.field("x_location", pa.float64()),
        pa.field("y_location", pa.float64()),
        pa.field("feature_name", pa.string()),
        pa.field("qv", pa.float64()),
        pa.field("nucleus_distance", pa.float64()),
    ])
    writer = _parquet_writer(out_path, schema)

    chunksize = 2_000_000
    read_kwargs = dict(compression="gzip") if tx_path.suffix == ".gz" else {}

    if headerless:
        # positional mapping: col6, col7, col9, (col8 as ndist), (col10 as qv if numeric)
        usecols = [5, 6, 7, 8, 9]  # 1-based: 6,7,8,9,10
        names = [f"v{i+1}" for i in range(max(usecols) + 1)]
        dtype_map = {"v10": str, "v9": str}
        for chunk in pd.read_csv(
            tx_path,
            header=None,
            names=names,
            usecols=usecols,
            chunksize=chunksize,
            low_memory=False,
            dtype=dtype_map,
            **read_kwargs,
        ):
            xg = pd.to_numeric(chunk.iloc[:, 0], errors="coerce")  # col6
            yg = pd.to_numeric(chunk.iloc[:, 1], errors="coerce")  # col7
            nd = pd.to_numeric(chunk.iloc[:, 2], errors="coerce")  # col8
            tgt = chunk.iloc[:, 3].astype(str)                      # col9
            qv = pd.to_numeric(chunk.iloc[:, 4], errors="coerce")  # col10 (may be str)

            df = pd.DataFrame({
                "x_location": xg * pixel_um,
                "y_location": yg * pixel_um,
                "feature_name": tgt,
                "qv": qv.fillna(np.inf).astype(float),
                "nucleus_distance": nd.fillna(0).astype(float),
            })
            # Optional: convert nucleus_distance to µm if it's in px
            if getattr(build_transcripts, "_nd_unit", "px") == "px":
                df["nucleus_distance"] = (df["nucleus_distance"].clip(lower=0)) * pixel_um
            else:
                df["nucleus_distance"] = df["nucleus_distance"].clip(lower=0)
            table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)
            writer.write_table(table)
    else:
        # named columns (use low_memory=False to avoid dtype surprises)
        read_kwargs2 = dict(read_kwargs)
        read_kwargs2.update(dict(low_memory=False))
        for chunk in pd.read_csv(tx_path, chunksize=chunksize, **read_kwargs2):
            low = {c.lower(): c for c in chunk.columns}

            def get_series(*ks):
                for k in ks:
                    if k in low:
                        return pd.to_numeric(chunk[low[k]], errors="coerce")
                return pd.Series(np.nan, index=chunk.index)

            def get_text(*ks):
                for k in ks:
                    if k in low:
                        return chunk[low[k]].astype(str)
                return pd.Series(["" for _ in range(len(chunk))])

            x_px = get_series("x_global_px", "x_px")
            y_px = get_series("y_global_px", "y_px")
            x_mm = get_series("x_global_mm", "x_mm")
            y_mm = get_series("y_global_mm", "y_mm")

            # final µm coordinates
            x_um = np.where(~x_px.isna(), x_px * pixel_um, x_mm * 1000.0)
            y_um = np.where(~y_px.isna(), y_px * pixel_um, y_mm * 1000.0)

            tgt = get_text("target", "gene", "feature_name")
            qv = get_series("qv")
            nd = get_series("nucleus_distance")

            df = pd.DataFrame({
                "x_location": pd.to_numeric(x_um, errors="coerce"),
                "y_location": pd.to_numeric(y_um, errors="coerce"),
                "feature_name": tgt,
                "qv": qv.fillna(np.inf).astype(float),
                "nucleus_distance": nd.fillna(0).astype(float),
            })
            if getattr(build_transcripts, "_nd_unit", "px") == "px":
                df["nucleus_distance"] = (df["nucleus_distance"].clip(lower=0)) * pixel_um
            else:
                df["nucleus_distance"] = df["nucleus_distance"].clip(lower=0)
            table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)
            writer.write_table(table)

    writer.close()
    print(f"[OK] Wrote {out_path}")


def build_cells_from_polygons(dataset_dir: Path, out_path: Path, pixel_um: float, overwrite: bool) -> None:
    if out_path.exists() and not overwrite:
        print(f"[SKIP] {out_path.name} exists (use --overwrite)")
        return

    poly_path = None
    for cand in [dataset_dir / f"{dataset_dir.name}-polygons.csv.gz", dataset_dir / f"{dataset_dir.name}_polygons.csv.gz"]:
        if cand.exists():
            poly_path = cand
            break
    if poly_path is None:
        raise FileNotFoundError("Could not find *-polygons.csv.gz under flatFiles")

    headerless = is_headerless_csv_gz(poly_path)
    print(f"[INFO] Building {out_path.name} from {poly_path.name}; headerless={headerless}")

    chunksize = 4_000_000
    read_kwargs = dict(compression="gzip")

    # accumulators for centroid = mean(x,y) per (fov, cell_id)
    sums = {}
    counts = {}

    if headerless:
        usecols = [0, 1, 5, 6]  # fov, cellid, x_px, y_px (1-based: 1,2,6,7)
        names = [f"v{i+1}" for i in range(max(usecols) + 1)]
        for chunk in pd.read_csv(poly_path, header=None, names=names, usecols=usecols, chunksize=chunksize, **read_kwargs):
            fov = pd.to_numeric(chunk.iloc[:, 0], errors="coerce").astype("Int64")
            cid = chunk.iloc[:, 1].astype(str)
            xg = pd.to_numeric(chunk.iloc[:, 2], errors="coerce")
            yg = pd.to_numeric(chunk.iloc[:, 3], errors="coerce")
            for fo, ce, x, y in zip(fov, cid, xg, yg):
                if pd.isna(fo) or pd.isna(x) or pd.isna(y):
                    continue
                key = (int(fo), ce)
                sx, sy = sums.get(key, (0.0, 0.0))
                c = counts.get(key, 0)
                sums[key] = (sx + float(x), sy + float(y))
                counts[key] = c + 1
    else:
        # named columns (accept several variants)
        for chunk in pd.read_csv(poly_path, chunksize=chunksize, **read_kwargs):
            low = {c.lower(): c for c in chunk.columns}
            def getcol(*ks):
                for k in ks:
                    if k in low:
                        return chunk[low[k]]
                return None
            fov = pd.to_numeric(getcol("fov"), errors="coerce").astype("Int64")
            cid = getcol("cellid", "cell_id").astype(str)
            xg = pd.to_numeric(getcol("x_global_px", "x_px", "x"), errors="coerce")
            yg = pd.to_numeric(getcol("y_global_px", "y_px", "y"), errors="coerce")
            for fo, ce, x, y in zip(fov, cid, xg, yg):
                if pd.isna(fo) or pd.isna(x) or pd.isna(y):
                    continue
                key = (int(fo), ce)
                sx, sy = sums.get(key, (0.0, 0.0))
                c = counts.get(key, 0)
                sums[key] = (sx + float(x), sy + float(y))
                counts[key] = c + 1

    if not counts:
        raise RuntimeError("No polygon points accumulated; cannot compute centroids.")

    rows = []
    for (fo, ce), (sx, sy) in sums.items():
        c = counts[(fo, ce)]
        rows.append({
            "cell_id": ce,
            "fov": int(fo),
            "x_centroid": (sx / c) * pixel_um,
            "y_centroid": (sy / c) * pixel_um,
        })
    df = pd.DataFrame(rows)
    table = pa.Table.from_pandas(df, preserve_index=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path, compression="zstd")
    print(f"[OK] Wrote {out_path} (n={len(df)})")


def build_segmentation_from_polygons(
    dataset_dir: Path,
    out_path: Path,
    pixel_um: float,
    overwrite: bool,
    *,
    # Compatibility mode: faithfully replicates the original script output (only coordinate scaling; no polyline break enhancement)
    compat: bool = True,
    # The following parameters apply only when compat=False to enhance polyline segmentation for visualization (avoid cross-FOV links)
    break_mult: float = 10.0,
    break_px: Optional[float] = None,
    mode: str = "cell",
    join_chunks: bool = True,
    subsample: int = 1,
) -> None:
    """Rebuild segmentation boundaries from polygons CSV.GZ into Parquet in µm.

    Adds polyline structure to avoid cross-FOV long lines when plotting by providing
    per-cell path breaks and vertex order.

    Columns:
      - cell_id, label_id, fov
      - vertex_x (µm), vertex_y (µm)
      - path_id (int32): increments at discontinuities within the same (fov,cell)
      - vertex_order (int32): index within (fov,cell,path)
    """
    if out_path.exists() and not overwrite:
        print(f"[SKIP] {out_path.name} exists (use --overwrite)")
        return

    poly_path = None
    for cand in [dataset_dir / f"{dataset_dir.name}-polygons.csv.gz", dataset_dir / f"{dataset_dir.name}_polygons.csv.gz"]:
        if cand.exists():
            poly_path = cand
            break
    if poly_path is None:
        raise FileNotFoundError("Could not find *-polygons.csv.gz under flatFiles")

    headerless = is_headerless_csv_gz(poly_path)
    print(f"[INFO] Building {out_path.name} from {poly_path.name}; headerless={headerless}")
    # Compatibility mode: schema matches the original script (no path_id / vertex_order)
    if compat:
        schema = pa.schema([
            pa.field("cell_id", pa.string()),
            pa.field("label_id", pa.string()),
            pa.field("vertex_x", pa.float64()),
            pa.field("vertex_y", pa.float64()),
            pa.field("fov", pa.int64()),
        ])
    else:
        if mode not in {"cell", "fov"}:
            raise ValueError("segmentation build: mode must be 'cell' or 'fov'")
        subsample = max(1, int(subsample))
        if subsample > 1:
            print(f"[INFO] Segmentation vertex subsampling: keep every {subsample} point (fast plotting mode)")
        schema = pa.schema([
            pa.field("cell_id", pa.string()),
            pa.field("label_id", pa.string()),
            pa.field("vertex_x", pa.float64()),
            pa.field("vertex_y", pa.float64()),
            pa.field("fov", pa.int64()),
            pa.field("path_id", pa.int32()),
            pa.field("vertex_order", pa.int32()),
        ])
    writer = _parquet_writer(out_path, schema)

    chunksize = 4_000_000
    read_kwargs = dict(compression="gzip")

    # Compatibility mode: read/write exactly like the original script, no segmentation enhancement
    if compat:
        if headerless:
            usecols = [0, 1, 5, 6]
            names = [f"v{i+1}" for i in range(max(usecols) + 1)]
            for chunk in pd.read_csv(poly_path, header=None, names=names, usecols=usecols, chunksize=chunksize, **read_kwargs):
                fov = pd.to_numeric(chunk.iloc[:, 0], errors="coerce").astype("Int64")
                cid = chunk.iloc[:, 1].astype(str)
                xg = pd.to_numeric(chunk.iloc[:, 2], errors="coerce")
                yg = pd.to_numeric(chunk.iloc[:, 3], errors="coerce")
                df = pd.DataFrame({
                    "cell_id": cid,
                    "label_id": cid,
                    "vertex_x": xg * pixel_um,
                    "vertex_y": yg * pixel_um,
                    "fov": fov.astype("Int64"),
                }).dropna(subset=["vertex_x", "vertex_y", "fov"])
                table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)
                writer.write_table(table)
        else:
            for chunk in pd.read_csv(poly_path, chunksize=chunksize, **read_kwargs):
                low = {c.lower(): c for c in chunk.columns}
                def getcol(*ks):
                    for k in ks:
                        if k in low:
                            return chunk[low[k]]
                    return None
                fov = pd.to_numeric(getcol("fov"), errors="coerce").astype("Int64")
                cid = getcol("cellid", "cell_id").astype(str)
                xg = pd.to_numeric(getcol("x_global_px", "x_px", "x"), errors="coerce")
                yg = pd.to_numeric(getcol("y_global_px", "y_px", "y"), errors="coerce")
                df = pd.DataFrame({
                    "cell_id": cid,
                    "label_id": cid,
                    "vertex_x": xg * pixel_um,
                    "vertex_y": yg * pixel_um,
                    "fov": fov.astype("Int64"),
                }).dropna(subset=["vertex_x", "vertex_y", "fov"])
                table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)
                writer.write_table(table)
        writer.close()
        print(f"[OK] Wrote {out_path}")
        return

    # --- Non-compat mode: estimate break threshold to enhance polylines for plotting (optional) ---
    # --- estimate polyline break threshold (pixels) from a sample ---
    def _estimate_break_px() -> float:
        sample_rows = 2_000_000
        if headerless:
            usecols = [0, 1, 5, 6]  # fov, cellid, x_px, y_px
            names = [f"v{i+1}" for i in range(max(usecols) + 1)]
            df = pd.read_csv(poly_path, header=None, names=names, usecols=usecols, nrows=sample_rows, **read_kwargs)
            fov = pd.to_numeric(df.iloc[:, 0], errors="coerce").astype("Int64").fillna(-1).astype(int)
            cid = df.iloc[:, 1].astype(str)
            xg = pd.to_numeric(df.iloc[:, 2], errors="coerce").to_numpy()
            yg = pd.to_numeric(df.iloc[:, 3], errors="coerce").to_numpy()
        else:
            df = pd.read_csv(poly_path, nrows=sample_rows, **read_kwargs)
            low = {c.lower(): c for c in df.columns}
            def getcol(*ks):
                for k in ks:
                    if k in low:
                        return df[low[k]]
                return None
            fov = pd.to_numeric(getcol("fov"), errors="coerce").astype("Int64").fillna(-1).astype(int)
            cid = getcol("cellid", "cell_id").astype(str)
            xg = pd.to_numeric(getcol("x_global_px", "x_px", "x"), errors="coerce").to_numpy()
            yg = pd.to_numeric(getcol("y_global_px", "y_px", "y"), errors="coerce").to_numpy()
        same = (fov.values[1:] == fov.values[:-1]) & (cid.values[1:] == cid.values[:-1])
        dx = (xg[1:] - xg[:-1])[same]
        dy = (yg[1:] - yg[:-1])[same]
        if dx.size == 0:
            return 32.0
        d = np.hypot(dx, dy)
        med = float(np.median(d[np.isfinite(d)])) if np.isfinite(d).any() else 1.0
        med = max(med, 1.0)
        return float(med * break_mult)

    eff_break_px = float(break_px) if break_px is not None else _estimate_break_px()
    print(f"[INFO] Segmentation polyline break threshold: {eff_break_px:.2f} px (mult={break_mult}, override_px={break_px})")

    # Cross-chunk carry-over state: used only when stitching across chunks
    prev_state: dict[tuple[int, str], tuple[float, float, int, int]] = {}

    def process_chunk(df_in: pd.DataFrame, *, is_headerless: bool) -> pd.DataFrame:
        # Optional vertex subsampling (row-level sampling), applied before later calculations
        if subsample > 1 and len(df_in) > subsample:
            df_in = df_in.iloc[::subsample, :].reset_index(drop=True)

        # Extract columns
        if is_headerless:
            fov = pd.to_numeric(df_in.iloc[:, 0], errors="coerce").astype("Int64")
            cid_raw = df_in.iloc[:, 1].astype(str)
            xg = pd.to_numeric(df_in.iloc[:, 2], errors="coerce")
            yg = pd.to_numeric(df_in.iloc[:, 3], errors="coerce")
        else:
            low = {c.lower(): c for c in df_in.columns}
            def getcol(*ks):
                for k in ks:
                    if k in low:
                        return df_in[low[k]]
                return None
            fov = pd.to_numeric(getcol("fov"), errors="coerce").astype("Int64")
            # Always read the original cell_id (even in fast mode) to align with centroids
            cid_raw = getcol("cellid", "cell_id").astype(str)
            xg = pd.to_numeric(getcol("x_global_px", "x_px", "x"), errors="coerce")
            yg = pd.to_numeric(getcol("y_global_px", "y_px", "y"), errors="coerce")

        # Normalize to a common DataFrame
        # Keep the real cell_id values; fast mode only affects grouping logic, not column contents
        cid = cid_raw

        df = pd.DataFrame({
            "fov": fov,
            "cell_id": cid,
            "x_px": xg,
            "y_px": yg,
        }).dropna(subset=["fov", "x_px", "y_px"]).reset_index(drop=True)

        if df.empty:
            return pd.DataFrame(columns=["cell_id", "label_id", "vertex_x", "vertex_y", "fov", "path_id", "vertex_order"])  # empty

        # -- compute breakpoints --
        fovv = df["fov"].to_numpy()
        cids = df["cell_id"].to_numpy()
        xv = df["x_px"].to_numpy(dtype=float); yv = df["y_px"].to_numpy(dtype=float)

        if mode == "cell":
            same_group = np.zeros(len(df), dtype=bool)
            same_group[1:] = (fovv[1:] == fovv[:-1]) & (cids[1:] == cids[:-1])
        else:
            # Fast mode: only connect when the FOV is unchanged
            same_group = np.zeros(len(df), dtype=bool)
            same_group[1:] = (fovv[1:] == fovv[:-1])

        dx = np.zeros(len(df), dtype=np.float32); dy = np.zeros(len(df), dtype=np.float32)
        dx[1:] = (xv[1:] - xv[:-1]); dy[1:] = (yv[1:] - yv[:-1])
        step = np.hypot(dx, dy)
        break_here = (~same_group) | (step > eff_break_px)
        break_here[0] = True  # Force a new path at the start of each chunk to avoid cross-chunk connections

        # -- compute path_id / vertex_order --
        local_path = np.cumsum(break_here.astype(np.int32)) - 1
        last_break_index = np.maximum.accumulate(np.where(break_here, np.arange(len(df)), -1))
        idx = np.arange(len(df))
        local_order = (idx - last_break_index).astype(np.int32)

        # Optional: stitch across chunks only when mode='cell' and join_chunks is True
        if mode == "cell" and join_chunks and len(df) > 0:
            run_starts = np.where(break_here)[0]
            if 0 not in run_starts:
                run_starts = np.concatenate([[0], run_starts])
            for si in run_starts:
                k = (int(fovv[si]), str(cids[si]))
                if k in prev_state:
                    px_prev, py_prev, pid_prev, ord_prev = prev_state[k]
                    dist = float(np.hypot(xv[si] - px_prev, yv[si] - py_prev))
                    if dist <= eff_break_px:
                        cur_pid = int(local_path[si])
                        local_path[local_path >= cur_pid] += (pid_prev - cur_pid)
                        cur_ord = int(local_order[si])
                        local_order[si:] += ((ord_prev + 1) - cur_ord)
                    else:
                        local_path[si:] += (pid_prev + 1)

            # Update cross-chunk state: keep only the last point of each (fov, cell) in this chunk
            uniq_keys = set(zip(fovv.tolist(), cids.tolist()))
            for k in uniq_keys:
                mask = (fovv == k[0]) & (cids == k[1])
                if not np.any(mask):
                    continue
                last_idx = int(np.flatnonzero(mask)[-1])
                prev_state[k] = (float(xv[last_idx]), float(yv[last_idx]), int(local_path[last_idx]), int(local_order[last_idx]))

        out = pd.DataFrame({
            "cell_id": df["cell_id"].astype(str).to_numpy(),
            "label_id": df["cell_id"].astype(str).to_numpy(),
            "vertex_x": xv * float(pixel_um),
            "vertex_y": yv * float(pixel_um),
            "fov": df["fov"].to_numpy(dtype=np.int64),
            "path_id": local_path.astype(np.int32),
            "vertex_order": local_order.astype(np.int32),
        })
        return out

    # stream and write
    if headerless:
        usecols = [0, 1, 5, 6]
        names = [f"v{i+1}" for i in range(max(usecols) + 1)]
        reader = pd.read_csv(poly_path, header=None, names=names, usecols=usecols, chunksize=chunksize, **read_kwargs)
    else:
        reader = pd.read_csv(poly_path, chunksize=chunksize, **read_kwargs)

    for chunk in reader:
        out_df = process_chunk(chunk, is_headerless=headerless)
        if not out_df.empty:
            table = pa.Table.from_pandas(out_df, schema=schema, preserve_index=False)
            writer.write_table(table)

    writer.close()
    print(f"[OK] Wrote {out_path}")

def convert_segmentation_to_um(path: Path, pixel_um: float) -> None:
    if not path.exists():
        return
    backup = path.with_suffix(".px.parquet")
    if not backup.exists():
        path.replace(backup)
    else:
        # keep a second backup
        path.replace(path.with_suffix(".px2.parquet"))

    pf = pq.ParquetFile(backup)
    writer = None
    for i in range(pf.num_row_groups):
        rg = pf.read_row_group(i)
        cols = [c for c in rg.schema.names]
        # determine XY columns
        x_name = None
        y_name = None
        for cand in ("x", "x_location", "vertex_x"):
            if cand in cols:
                x_name = cand
                break
        for cand in ("y", "y_location", "vertex_y"):
            if cand in cols:
                y_name = cand
                break
        if x_name is None or y_name is None:
            raise ValueError(f"Cannot find X/Y columns in {backup}")
        arr = rg.to_pandas()
        arr[x_name] = pd.to_numeric(arr[x_name], errors="coerce") * pixel_um
        arr[y_name] = pd.to_numeric(arr[y_name], errors="coerce") * pixel_um
        # ensure label_id exists for geneSCOPE segmentation processing
        if "label_id" not in arr.columns:
            if "cell_id" in arr.columns:
                arr["label_id"] = arr["cell_id"].astype(str)
            else:
                arr["label_id"] = ""
        table = pa.Table.from_pandas(arr, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(path, table.schema, compression="zstd")
        writer.write_table(table)
    if writer is not None:
        writer.close()
    print(f"[OK] Converted {backup.name} -> {path.name} (µm)")


def main() -> None:
    ap = argparse.ArgumentParser(description="CosMx flatFiles -> geneSCOPE Parquets (µm)")
    ap.add_argument("--root", required=True, type=Path, help="CosMx run root (contains flatFiles/)")
    ap.add_argument("--dataset-id", default=None, help="Dataset name under flatFiles (auto if omitted)")
    ap.add_argument("--pixel-size-um", type=float, default=0.120280945, help="Micrometers per pixel")
    ap.add_argument("--build-transcripts", action="store_true")
    ap.add_argument("--build-cells", action="store_true")
    ap.add_argument("--convert-segmentation", action="store_true")
    ap.add_argument("--build-segmentation", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--nucleus-distance-unit", choices=["px","um"], default="px",
                    help="Interpret nucleus_distance in tx_file as pixels (convert to µm) or already in µm.")
    # Compatibility toggle: default True matches the original script (only unit scaling); set False to enable enhanced segmentation options
    ap.add_argument("--seg-compat", dest="seg_compat", action="store_true", default=True,
                    help="Compatibility mode with the original script: only coordinate scaling, no path_id/vertex_order (default)")
    ap.add_argument("--no-seg-compat", dest="seg_compat", action="store_false",
                    help="Disable compatibility mode to enable enhanced polyline segmentation options (below)")
    ap.add_argument("--seg-break-mult", type=float, default=10.0,
                    help="Takes effect when compatibility mode is off: break threshold = median step × mult (px)")
    ap.add_argument("--seg-break-px", type=float, default=None,
                    help="Takes effect when compatibility mode is off: absolute pixel threshold, overrides mult")
    ap.add_argument("--seg-mode", choices=["cell","fov"], default="cell",
                    help="Takes effect when compatibility mode is off: 'cell' enforces (fov,cell) grouping; 'fov' restricts within the same FOV")
    ap.add_argument("--no-join-chunks", dest="join_chunks", action="store_false", default=True,
                    help="Takes effect when compatibility mode is off: do not try to stitch paths across chunks (faster)")
    ap.add_argument("--seg-subsample", type=int, default=1,
                    help="Takes effect when compatibility mode is off: vertex subsampling (for plotting speed-up only)")
    args = ap.parse_args()

    root = args.root.resolve()
    if not (root / "flatFiles").exists():
        raise FileNotFoundError(f"{root}/flatFiles not found")
    dataset_dir = find_dataset_dir(root / "flatFiles")

    if args.build_transcripts:
        # communicate nd unit to build_transcripts via a function attribute
        build_transcripts._nd_unit = args.nucleus_distance_unit
        build_transcripts(root, dataset_dir, root / "transcripts.parquet", args.pixel_size_um, args.overwrite)

    if args.build_cells:
        build_cells_from_polygons(dataset_dir, root / "cells.parquet", args.pixel_size_um, args.overwrite)

    if args.build_segmentation:
        build_segmentation_from_polygons(
            dataset_dir,
            root / "segmentation_boundaries.parquet",
            args.pixel_size_um,
            args.overwrite,
            compat=args.seg_compat,
            break_mult=args.seg_break_mult,
            break_px=args.seg_break_px,
            mode=args.seg_mode,
            join_chunks=args.join_chunks,
            subsample=args.seg_subsample,
        )

    if args.convert_segmentation:
        # Avoid double-scaling if we also built segmentation in this run.
        names = ["cell_boundaries.parquet", "nucleus_boundaries.parquet"]
        if not args.build_segmentation:
            names.insert(0, "segmentation_boundaries.parquet")
        for name in names:
            convert_segmentation_to_um(root / name, args.pixel_size_um)

    print("[DONE] Conversion/check complete")


if __name__ == "__main__":
    main()
