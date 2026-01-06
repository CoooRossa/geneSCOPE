# Diagnose sendMaster(SIGPIPE) in createSCOPE Xenium ROI clipping

## Call chain (createSCOPE -> ROI clip)
- `createSCOPE()` -> `R/api_data_construction.R:105-117`
- `.create_scope()` (dispatches Xenium) -> `R/internal_ingest.R:3360-3403`
- `createSCOPE_xenium()` -> `R/internal_ingest.R:4173-4209`
- `build_scope_from_xenium()` -> `R/internal_ingest.R:3985-4055`
- `.clip_points_within_roi()` -> `R/internal_ingest.R:845-871`
- `.clip_points_to_region()` (parallel ROI clip) -> `R/internal_misc.R:1758-1917`
  - `mclapply()` fork backend -> `R/internal_misc.R:1863-1866`
- `.prefetch_roi_molecules()` -> `R/internal_ingest.R:1510-1545`
  - calls `.clip_points_to_region()` for molecule ROI clip -> `R/internal_ingest.R:1531-1537`

## What is parallelized
- The ROI clip splits `dt_pts` (centroids or prefetch molecules) into `chunk_ids` by row index
  (`split(seq_len(nrow(dt_pts)), ceiling(seq_len(nrow(dt_pts)) / chunk_size))`, `R/internal_misc.R:1786-1790`).
- Each worker converts the chunk to an `sf` point layer and runs `sf::st_within()` against the
  ROI `sfc` geometry (`R/internal_misc.R:1834-1844`).
- Fork backend (`mclapply`) is used on non-Windows; PSOCK now broadcasts the ROI geometry and
  ships chunk data in batches (no full `dt_pts` export).

## Likely root causes (code-backed)
1) Fork + `sf`/GEOS instability in workers
   - ROI clip workers call `sf::st_within()` inside `mclapply()` without protection or validation
     beyond `tryCatch` in R (`R/internal_misc.R:1834-1866`). If GEOS/s2 or the C stack aborts,
     the worker dies and `sendMaster()` sees SIGPIPE.
2) Memory pressure on large chunks / many workers
   - The code parallelizes by chunk but each worker allocates an `sf` object for its chunk
     (`st_as_sf` inside `clip_worker`, `R/internal_misc.R:1837-1844`). With large `chunk_size`
     (default 5e5) and many workers (`ncores_safe`), a forked worker can be OOM-killed.
3) Forking after Arrow threads are initialized
   - Arrow thread counts are set in `.load_dataset_dependencies()` (`R/internal_ingest.R:1430-1450`)
     and `arrow::open_dataset()` is called before ROI clipping (`R/internal_ingest.R:1495-1500`).
     Forking a process that already initialized multi-threaded Arrow code can be unsafe and lead
     to worker crashes, surfacing as SIGPIPE in `mclapply`.

## Recommended default behavior change
- Add `parallel_backend = "auto"` for ROI clip stages to choose the safest backend per environment:
  - Linux containers: default to `serial` for large tables, otherwise `psock`
  - Windows: default to `psock`
  - macOS / non-container Linux: default to `fork`
- If the selected parallel backend errors or workers return failures, retry once in serial and
  emit a warning explaining how to override.
- Rationale: preserves correctness and avoids silent worker death, while still allowing users to
  opt into faster parallel backends when stable.

## PSOCK chunk shipping (current)
- PSOCK ROI clip no longer exports full `dt_pts` to each worker. Instead, the master materializes
  per-chunk subsets and ships only those chunks to workers.
- ROI geometry is broadcast once to workers; only chunk data is transmitted per task.
- Batching limits in-flight chunk memory. Default `options(geneSCOPE.psock_batch_chunks = 16)`.
- Large-table heuristic for container defaults uses `options(geneSCOPE.large_table_rows = 5e6)`.

## Debug hooks
- `options(geneSCOPE.debug_parallel = TRUE)` logs chunk sizes, backend, and data size.
- `options(geneSCOPE.debug_parallel_dir = "/tmp")` writes per-worker error logs for ROI clip.

## Recommended settings
- Centroids (cells): `psock` is OK with default `chunk_size`.
- Molecules (Xenium): prefer `serial` in containers; if forcing `psock`, use smaller `chunk_size`
  and conservative `psock_batch_chunks` (e.g. 8-16).

## Patch summary (implemented)
- Added backend selection + serial fallback in `.clip_points_to_region()`.
- PSOCK now ships chunk data with batching; ROI is broadcast to workers (no full `dt_pts` export).
- Plumbed `parallel_backend` through `createSCOPE_xenium()` and `build_scope_from_xenium()` into
  ROI clip stages (`.clip_points_within_roi()` / `.prefetch_roi_molecules()`).
- Added lightweight container detection to support safe `auto` backend selection.
