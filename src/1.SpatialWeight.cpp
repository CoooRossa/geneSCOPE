// 1.SpatialWeight.cpp (2025-06-27)
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// Requires OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppEigen.h>
using namespace Rcpp;

/*------------------------------------------------------------------*/
/* 1. Parallel generation of Queen/Rook adjacency list            */
/*------------------------------------------------------------------*/
//' @title Build Queen or Rook Neighbourhood List in Parallel
//'
//' @description
//'   Creates a spatial neighbourhood list for a rectangular grid using either
//'   Queen (8-neighbour) or Rook (4-neighbour) connectivity. OpenMP
//'   parallelisation scales efficiently to millions of grid cells.
//'
//' @param nrow Integer. Number of rows in the grid.
//' @param ncol Integer. Number of columns in the grid.
//' @param queen Logical. TRUE (default) for Queen moves, FALSE for Rook.
//'
//' @return A list of integer vectors of class \code{nb}. Indices are 1-based
//'   and attributes \code{region.id} and \code{queen} are set.
//'
//' @details
//'   Offsets for rook and (optionally) diagonal neighbours are pre-computed and
//'   applied in an OpenMP \code{for} loop. Each thread accumulates a private
//'   neighbour vector which is swapped into the main container to avoid
//'   locking.
//'
//' @examples
//' \dontrun{
//' nb <- grid_nb_omp(1024, 1024, queen = TRUE)
//' head(nb[[1]])
//' }
// [[Rcpp::export]]
List grid_nb_omp(const int nrow,
                 const int ncol,
                 const bool queen = true)
{
    const int n = nrow * ncol;
    std::vector<std::vector<int>> nb(n);

    const int drow[8] = {0, 0, -1, 1, -1, -1, 1, 1};
    const int dcol[8] = {-1, 1, 0, 0, -1, 1, -1, 1};
    const int nOff = queen ? 8 : 4;

#pragma omp parallel for schedule(static)
    for (int idx = 0; idx < n; ++idx)
    {
        int r = idx / ncol, c = idx % ncol;
        std::vector<int> local_nb;
        local_nb.reserve(nOff);
        for (int k = 0; k < nOff; ++k)
        {
            int rr = r + drow[k], cc = c + dcol[k];
            if (rr >= 0 && rr < nrow && cc >= 0 && cc < ncol)
                local_nb.push_back(rr * ncol + cc + 1); // 1-based
        }
        nb[idx].swap(local_nb);
    }

    List out(n);
    for (int i = 0; i < n; ++i)
        out[i] = IntegerVector(nb[i].begin(), nb[i].end());
    out.attr("class") = "nb";
    out.attr("region.id") = Rcpp::seq(1, n);
    out.attr("queen") = queen;
    return out;
}

/*------------------------------------------------------------------*/
/* 2. Build binary sparse matrix (dgCMatrix) in parallel          */
/*------------------------------------------------------------------*/
//' @title Convert Neighbour List to Binary Sparse Matrix (dgCMatrix)
//'
//' @description
//'   Transforms an \code{nb} neighbourhood list into a column-compressed
//'   binary sparse matrix of class \code{dgCMatrix}. The routine streams the
//'   triplet representation in parallel with OpenMP and finally constructs the
//'   sparse matrix with Eigen.
//'
//' @param nb An object of class \code{nb}, typically returned by
//'   \code{grid_nb_omp()} or \pkg{spdep} functions.
//'
//' @return A \code{dgCMatrix} of dimension \eqn{n \times n} where \eqn{n} is
//'   the number of regions. Entries are \code{1} where two regions are
//'   adjacent and \code{0} elsewhere.
//'
//' @details
//'   The neighbour indices are 1-based; they are converted to 0-based before
//'   creation of the triplets for Eigen. Each thread accumulates its own
//'   vector of triplets to avoid synchronisation. After the parallel region
//'   all triplets are concatenated and the resulting sparse matrix is
//'   compressed in column-major order.
//'
//' @examples
//' \dontrun{
//' nb  <- grid_nb_omp(512, 512, queen = FALSE)
//' W   <- listw_B_omp(nb)
//' Matrix::nnzero(W) / length(nb)  # average degree
//' }
// [[Rcpp::export]]
SEXP listw_B_omp(const List nb)
{
    const int n = nb.size();

    /* 2.1 单线程深拷贝到 C++ 容器 —— 彻底隔离 R 内存 */
    std::vector<std::vector<int>> nb_cpp(n);
    for (int i = 0; i < n; ++i)
    {
        IntegerVector v = nb[i];
        nb_cpp[i].assign(v.begin(), v.end());
    }

    /* 2.2 线程私有三元组 */
#ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#else
    const int nThreads = 1;
#endif
    std::vector<std::vector<Eigen::Triplet<double>>> trip_tls(nThreads);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        int tid =
#ifdef _OPENMP
            omp_get_thread_num();
#else
            0;
#endif
        auto &vec = trip_tls[tid];
        vec.reserve(vec.size() + nb_cpp[i].size());
        for (int j : nb_cpp[i])
        {
            if (j == 0)
                continue;
            vec.emplace_back(i, j - 1, 1.0); // 0-based
        }
    }

    /* 2.3 合并三元组 */
    size_t total = 0;
    for (auto &v : trip_tls)
        total += v.size();
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(total);
    for (auto &v : trip_tls)
    {
        std::move(v.begin(), v.end(), std::back_inserter(trip));
        std::vector<Eigen::Triplet<double>>().swap(v); // 释放线程缓存
    }

    /* 2.4 构造压缩稀疏矩阵 */
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> W(n, n);
    W.setFromTriplets(trip.begin(), trip.end());
    W.makeCompressed();

    return Rcpp::wrap(W); // dgCMatrix
}

// 兼容新命名但调用原版函数
// [[Rcpp::export]]
SEXP grid_nb(const int nrow,
             const int ncol,
             const bool queen = true)
{
    return grid_nb_omp(nrow, ncol, queen);
}

// [[Rcpp::export]]
SEXP nb2mat(SEXP nb)
{
    return listw_B_omp(nb);
}