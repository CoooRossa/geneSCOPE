// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
arma::mat leeL_perm_ge_cpp(const arma::mat  &Xz,
                           const arma::sp_mat &W,
                           const arma::umat &idx_mat,
                           const arma::mat  &L_ref)
{
    const uword ncell = Xz.n_rows,
                ngen  = Xz.n_cols,
                B     = idx_mat.n_cols;
    const double S0   = accu(W);

    mat geCnt(ngen, ngen, fill::zeros);

#pragma omp parallel
{
    mat localCnt(ngen, ngen, fill::zeros);

#pragma omp for schedule(dynamic)
    for (uword b = 0; b < B; ++b) {
        mat  Xp  = Xz.rows( idx_mat.col(b) );

        vec den  = sum(square(Xp), 0).t();

        for (uword f = 0; f < ngen; ++f) {
            vec z   = Xp.col(f);
            vec Wz  = W * z;
            vec num = Xp.t() * Wz;
            double dz2  = accu(square(z));
            vec denom   = sqrt(dz2 * den);
            vec Ltmp    = (static_cast<double>(ncell) / S0) * (num / denom);

            // ---- 关键修改：转换成 double 0/1 向量再累加 ----
            localCnt.col(f) += arma::conv_to<arma::vec>::from(
+                                   abs(Ltmp) >= abs(L_ref.col(f)) );
        }
    }

#pragma omp critical
    geCnt += localCnt;
}
    return geCnt;
}