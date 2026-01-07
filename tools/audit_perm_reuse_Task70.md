# Task70 Permutation Reuse Audit

Date: 2026-01-08
Repo: /Users/haenolabcho/Documents/Code/geneSCOPE-2025-01-05/geneSCOPE

## Candidates found

- C++
  - `src/2.LeeL.cpp` exports `lee_perm`, `lee_perm_block` (permutation counts for Lee's L)
  - `src/7.DeltaLRPerm.cpp` exports `delta_lr_perm*` variants (permutation counts for Delta = Lee's L - Pearson r)
- R wrappers
  - `R/RcppExports.R` exposes `.Call` wrappers for `lee_perm*` and `delta_lr_perm*`
  - `R/native_permutation.R` provides `.delta_perm_pairs()` helper for those functions

## Reuse criteria check

Required for reuse:
- Accept membership label permutation (module size preserved)
- Return permuted memberships or null statistics for member/module scores
- Not tied to Lee's L / DeltaLR-specific computation
- Callable from R without changing existing wrappers

Findings:
- Existing C++ functions take `Xz`, `W`, and `idx_mat` (cell index permutations) and compute Lee's L / DeltaLR exceedance counts.
- They do **not** accept membership labels, do **not** permute module labels, and their outputs are Lee's L / DeltaLR-specific.
- Adapting them would require changing their semantics/inputs, which is not allowed for this task.

## Conclusion

**Not reusable.**

Decision: implement a new C++ permutation routine dedicated to label permutation null for module-quality scoring (member/module metrics), with an R wrapper that is used only by the new Task70 API.
