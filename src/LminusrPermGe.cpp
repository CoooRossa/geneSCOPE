// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
arma::mat leeL_minusR_perm_ge_cpp(const arma::mat &Xz,
                                  const arma::sp_mat &W,
                                  const arma::umat &idx_mat,
                                  const arma::mat &Delta_ref)
{
    const uword ncell = Xz.n_rows,
                ngen = Xz.n_cols,
                B = idx_mat.n_cols;

    const double S0 = accu(W); // Σ W_ij
    arma::mat geCnt(ngen, ngen, fill::zeros);

#pragma omp parallel
    {
        arma::mat localCnt(ngen, ngen, fill::zeros);

#pragma omp for schedule(dynamic)
        for (uword b = 0; b < B; ++b)
        {

            arma::mat Xp = Xz.rows(idx_mat.col(b)); // ncell × ngen

            arma::vec den = sum(square(Xp), 0).t(); // ngen × 1

            for (uword f = 0; f < ngen; ++f)
            {

                arma::vec z = Xp.col(f); // ncell × 1
                double dz2 = accu(square(z));

                arma::vec Wz = W * z;         // ncell × 1
                arma::vec numL = Xp.t() * Wz; // ngen × 1
                arma::vec denomL = sqrt(dz2 * den);
                arma::vec Ltmp = (static_cast<double>(ncell) / S0) * (numL / denomL);

                arma::vec numR = Xp.t() * z; // ngen × 1
                arma::vec rtmp = numR / static_cast<double>(ncell);

                arma::vec Dtmp = Ltmp - rtmp;

                localCnt.col(f) += arma::conv_to<arma::vec>::from(
                    abs(Dtmp) >= abs(Delta_ref.col(f)));
            }
        }

#pragma omp critical
        geCnt += localCnt;
    }

    return geCnt; // P = (geCnt + 1)/(B + 1) in R
}