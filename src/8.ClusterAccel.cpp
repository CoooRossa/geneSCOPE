// Lightweight C++ accelerators for clustering consensus and CMH mapping
// Uses OpenMP when available (Makevars config handles enable/disable)

#include <Rcpp.h>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::export(name="consensus_on_edges_omp")]]
IntegerVector consensus_on_edges_omp(IntegerVector ei, IntegerVector ej,
                                     IntegerMatrix memb,
                                     int n_threads = 1) {
  const R_xlen_t E = ei.size();
  if (E != ej.size()) stop("ei/ej length mismatch");
  const int R = memb.ncol();
  if (R <= 0) return IntegerVector(E); // all zeros

  IntegerVector out(E);

  // Clamp threads if OpenMP present; otherwise ignore
  #ifdef _OPENMP
    if (n_threads < 1) n_threads = 1;
  #else
    n_threads = 1;
  #endif

  #pragma omp parallel for schedule(static) num_threads(n_threads)
  for (R_xlen_t e = 0; e < E; ++e) {
    const int u = ei[e] - 1; // 1-based to 0-based
    const int v = ej[e] - 1;
    int cnt = 0;
    for (int r = 0; r < R; ++r) {
      if (memb(u, r) == memb(v, r)) ++cnt;
    }
    out[e] = cnt;
  }
  return out;
}

static inline uint64_t pair_key32(uint32_t a, uint32_t b) {
  if (a > b) std::swap(a, b);
  return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

// [[Rcpp::export(name="cmh_lookup_rcpp")]]
NumericVector cmh_lookup_rcpp(IntegerVector ei, IntegerVector ej,
                              IntegerVector si, IntegerVector sj, NumericVector sw,
                              double fallback = NA_REAL, int n_threads = 1) {
  const R_xlen_t E = ei.size();
  if (E != ej.size()) stop("ei/ej length mismatch");
  if (si.size() != sj.size() || si.size() != sw.size()) stop("sim edge vectors mismatch");

  // Build map for similarity edges
  std::unordered_map<uint64_t, double> mp;
  mp.reserve(static_cast<size_t>(si.size() * 1.3));
  for (R_xlen_t k = 0; k < si.size(); ++k) {
    const int a = si[k] - 1;
    const int b = sj[k] - 1;
    if (a < 0 || b < 0) continue;
    const uint64_t key = pair_key32(static_cast<uint32_t>(a), static_cast<uint32_t>(b));
    mp[key] = sw[k];
  }

  NumericVector out(E);

  #ifdef _OPENMP
    if (n_threads < 1) n_threads = 1;
  #else
    n_threads = 1;
  #endif

  #pragma omp parallel for schedule(static) num_threads(n_threads)
  for (R_xlen_t e = 0; e < E; ++e) {
    const int a = ei[e] - 1;
    const int b = ej[e] - 1;
    double w = fallback;
    if (a >= 0 && b >= 0) {
      const uint64_t key = pair_key32(static_cast<uint32_t>(a), static_cast<uint32_t>(b));
      auto it = mp.find(key);
      if (it != mp.end()) w = it->second;
    }
    out[e] = w;
  }
  return out;
}

// ---- CI95 edge filtering (moved from 9.CI95Accel.cpp) ----------------------

static inline double interp_linear(const double *x, const double *y, int n, double v) {
  // x is ascending. rule=2: clamp to boundary values outside range.
  if (v <= x[0]) return y[0];
  if (v >= x[n-1]) return y[n-1];
  // binary search
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) >> 1;
    if (x[mid] <= v) lo = mid; else hi = mid;
  }
  const double x0 = x[lo], x1 = x[hi];
  const double y0 = y[lo], y1 = y[hi];
  const double t = (v - x0) / (x1 - x0);
  return y0 + t * (y1 - y0);
}

// rule: 0 = remove_within (drop if L <= hi(r));
//        1 = remove_outside (drop if L < lo(r) || L > hi(r))
// rMat: correlation matrix (g_total x g_total)
// ridx: length G_kept map from kept_genes index -> rMat row index (1-based)
// ei/ej in 1-based index into kept_genes order; L_vals length E.
// [[Rcpp::export(name="ci95_drop_mask_edges_omp")]]
LogicalVector ci95_drop_mask_edges_omp(NumericMatrix rMat,
                                       IntegerVector ridx,
                                       IntegerVector ei,
                                       IntegerVector ej,
                                       NumericVector L_vals,
                                       NumericVector xp,
                                       NumericVector lo95,
                                       NumericVector hi95,
                                       int rule = 0,
                                       int n_threads = 1) {
  const R_xlen_t E = ei.size();
  if (E != ej.size() || E != L_vals.size()) stop("edge vectors length mismatch");
  const int ngrid = xp.size();
  if (ngrid < 2 || lo95.size() != ngrid || hi95.size() != ngrid) stop("invalid CI curves");

  const double* xp_p = REAL(xp);
  const double* lo_p = REAL(lo95);
  const double* hi_p = REAL(hi95);
  const int* ridx_p = INTEGER(ridx);

  LogicalVector drop(E);

  #ifdef _OPENMP
    if (n_threads < 1) n_threads = 1;
  #else
    n_threads = 1;
  #endif

  #pragma omp parallel for schedule(static) num_threads(n_threads)
  for (R_xlen_t k = 0; k < E; ++k) {
    const int gi = ridx_p[ ei[k] - 1 ] - 1; // to 0-based
    const int gj = ridx_p[ ej[k] - 1 ] - 1;
    if (gi < 0 || gj < 0) { drop[k] = false; continue; }
    const double r = rMat(gi, gj);
    const double hlo = interp_linear(xp_p, lo_p, ngrid, r);
    const double hhi = interp_linear(xp_p, hi_p, ngrid, r);
    const double L  = L_vals[k];
    bool d;
    if (rule == 0) {
      d = (L <= hhi);
    } else {
      d = (L < hlo) || (L > hhi);
    }
    drop[k] = d;
  }

  return drop;
}

