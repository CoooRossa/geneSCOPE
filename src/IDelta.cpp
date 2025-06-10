// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
arma::vec idelta_cpp(const arma::sp_mat &G)
{

    // 1) G: g×n, sparse matrix
    arma::mat denseG = arma::mat(G);
    arma::vec tot = sum(denseG, 1); // g×1

    // 2) Σ x (x-1)
    arma::mat prod = denseG % (denseG - 1.0);
    arma::vec xi_xi1 = sum(prod, 1);

    // 3) Iδ_raw
    const double q = static_cast<double>(G.n_cols);
    arma::vec id_raw = (q * xi_xi1) / (tot % (tot - 1));

    // 4) tot < 2 → NaN
    id_raw.elem(find(tot < 2.0)).fill(datum::nan);

    return id_raw; // g×1
}
