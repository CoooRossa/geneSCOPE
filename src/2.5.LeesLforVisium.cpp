// 5.VisiumAccel.cpp (2025-10-26)
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// Sparse × Dense multiply: WX = W %*% X
// For numerical stability, return double; f32 denotes the optional path only.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat spmm_dgc_dense_f64(const arma::mat& X, const arma::sp_mat& W, const int n_threads = 1)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    const uword n = X.n_rows;
    const uword S = X.n_cols;
    arma::mat Y(n, S, fill::zeros);

#pragma omp parallel for schedule(static)
    for (uword j = 0; j < S; ++j) {
        Y.col(j) = W * X.col(j);
    }
    return Y;
}

// [[Rcpp::export]]
arma::mat spmm_dgc_dense_f32(const arma::mat& X, const arma::sp_mat& W, const int n_threads = 1)
{
    // Simplification: compute and return in double to avoid extra copies and casts.
    // For n<=2000 and S<=1e4, memory use is typically acceptable.
    return spmm_dgc_dense_f64(X, W, n_threads);
}

// Convert uint64_t to hex string (with 0x prefix)
static inline std::string to_hex64(uint64_t v) {
    std::ostringstream oss;
    oss << "0x" << std::hex << std::setw(16) << std::setfill('0') << (unsigned long long)v;
    return oss.str();
}

// [[Rcpp::export]]
CharacterMatrix rp_sign_bits(const arma::mat& X,
                             const int bits = 12,
                             const int n_tables = 6,
                             const int seed = 1,
                             const int n_threads = 1)
{
    if (bits <= 0 || bits > 60) stop("bits must be in [1, 60]");
    if (n_tables <= 0) stop("n_tables must be positive");

    const uword n = X.n_rows;  // spots
    const uword S = X.n_cols;  // genes

    CharacterMatrix out(S, n_tables);

    // For each table and each bit, generate a length-n Rademacher sign vector
    // and take signed dot products with each column, packing to 64-bit bucket IDs

#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

#pragma omp parallel for schedule(static)
    for (int t = 0; t < n_tables; ++t) {
        std::mt19937_64 rng((uint64_t)seed + 0x9E3779B97F4A7C15ULL * (t + 1));
        // Pre-generate bits × n sign matrix (bool), true=+1, false=-1
        std::vector< std::vector<uint8_t> > signs(bits, std::vector<uint8_t>(n));
        for (int b = 0; b < bits; ++b) {
            // Simple reproducible Rademacher via RNG low bit
            for (uword i = 0; i < n; ++i) {
                uint64_t r = rng();
                signs[b][i] = (r & 1ULL) ? 1 : 0; // 1=>+1, 0=>-1
            }
        }
        // Compute bucket ID for each gene (column)
        for (uword g = 0; g < S; ++g) {
            uint64_t h = 0ULL;
            for (int b = 0; b < bits; ++b) {
                // Signed projection for the b-th hyperplane
                const uint8_t* sg = signs[b].data();
                const double* xg = X.colptr(g);
                double acc = 0.0;
                for (uword i = 0; i < n; ++i) {
                    acc += sg[i] ? xg[i] : -xg[i];
                }
                const uint64_t bit = (acc >= 0.0) ? 1ULL : 0ULL;
                h |= (bit << b);
            }
            out(g, t) = to_hex64(h);
        }
    }
    return out;
}

// Legacy candidate interface. It cannot compute canonical Lee's L because it
// receives WX but not W (and therefore cannot determine S2). Keep the ABI so
// existing callers fail with a precise error instead of silently receiving a
// bivariate-Moran quantity mislabeled as Lee's L.
// Input: row_ptr (S+1, 1-based), indices (1-based)
// Output: CSR trimmed to Top-K (triplets: row_ptr, indices, values)
// [[Rcpp::export]]
Rcpp::List leeL_topk_candidates(const arma::mat& X,
                                const arma::mat& WX,
                                const Rcpp::IntegerVector& row_ptr,
                                const Rcpp::IntegerVector& indices,
                                const int K_keep = 100,
                                const int n_threads = 1)
{
    (void)X;
    (void)WX;
    (void)row_ptr;
    (void)indices;
    (void)K_keep;
    (void)n_threads;
    stop(
        "leeL_topk_candidates is disabled: its legacy ABI lacks W/S2 and cannot "
        "compute canonical Lee's L. Use computeL() or lee_L()."
    );
    return Rcpp::List::create();
}
