// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// ------------------------------------------------------------
// 1. Single-run Lee's L (corrected)
// ------------------------------------------------------------
// [[Rcpp::export]]
arma::mat leeL_cpp(const arma::mat &Xz, const arma::sp_mat &W)
{
    const unsigned int ncell = Xz.n_rows, ngen = Xz.n_cols;
    arma::mat L(ngen, ngen, arma::fill::zeros);

    double S0 = arma::accu(W); // sum of weights
    if (S0 == 0)
        return L;

    arma::vec den = arma::sum(arma::square(Xz), 0).t(); // squared norm of each gene vector

#pragma omp parallel for schedule(static)
    for (unsigned int f = 0; f < ngen; ++f)
    {
        arma::vec z = Xz.col(f);
        arma::vec Wz = W * z;                     // W * z_f
        arma::vec num = Xz.t() * Wz;              // Xᵀ * (W * z_f) for all genes
        double dz2 = arma::accu(arma::square(z)); // squared norm of z_f
        arma::vec denom = arma::sqrt(dz2 * den);  // sqrt(‖z_f‖² * ‖z_g‖²) for each gene g
        L.col(f) = (static_cast<double>(ncell) / S0) * (num / denom);
    }
    return L;
}

// ------------------------------------------------------------
// 2. Batch permutations: returns g × g × B cube
// ------------------------------------------------------------
// [[Rcpp::export]]
arma::cube leeL_perm_cpp(const arma::mat &Xz, const arma::sp_mat &W,
                         const arma::umat &idx_mat)
{
    const unsigned int ncell = Xz.n_rows,
                       ngen = Xz.n_cols,
                       B = idx_mat.n_cols; // number of permutations
    double S0 = arma::accu(W);
    arma::cube out(ngen, ngen, B, arma::fill::zeros);

#pragma omp parallel for schedule(dynamic)
    for (unsigned int b = 0; b < B; ++b)
    {
        arma::uvec idx = idx_mat.col(b);
        arma::mat Xp = Xz.rows(idx); // permuted Xz
        arma::vec den = arma::sum(arma::square(Xp), 0).t();
        for (unsigned int f = 0; f < ngen; ++f)
        {
            arma::vec z = Xp.col(f);
            arma::vec Wz = W * z;
            arma::vec num = Xp.t() * Wz;
            double dz2 = arma::accu(arma::square(z));
            arma::vec denom = arma::sqrt(dz2 * den);
            out.slice(b).col(f) = (static_cast<double>(ncell) / S0) * (num / denom);
        }
    }
    return out; // return a g × g × B cube
}

// ------------------------------------------------------------
// 3. Bulk βx / βy estimation
// ------------------------------------------------------------
// [[Rcpp::export]]
arma::mat grad_betas_cpp(const arma::mat &Xz,
                         const arma::vec &cx,
                         const arma::vec &cy)
{
    const unsigned int n = Xz.n_rows,
                       g = Xz.n_cols;

    arma::mat X(n, 3, arma::fill::ones); // design matrix [1  cx  cy]
    X.col(1) = cx;
    X.col(2) = cy;

    // Armadillo::solve() automatically selects QR or SVD for least squares
    arma::mat coef = arma::solve(X, Xz); // solves for coefficients, result is 3 × g

    return coef.rows(1, 2).t(); // returns g × 2 matrix with (βx, βy) for each gene
}
