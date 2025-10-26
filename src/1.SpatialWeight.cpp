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
/* 1b. Parallel generation of Hex adjacency (odd-r / even-r)       */
/*------------------------------------------------------------------*/
//' @title Build Hexagonal (6-neighbour) Neighbourhood List in Parallel
//'
//' @description
//'   Creates a spatial neighbourhood list for a hex tiling laid out on a
//'   rectangular index with row-wise offsets (odd-r or even-r). See
//'   https://www.redblobgames.com/grids/hex-grids/ for conventions.
//'
//' @param nrow Integer. Number of rows in the logical lattice.
//' @param ncol Integer. Number of columns in the logical lattice.
//' @param oddr Logical. If TRUE, use odd-r layout (odd rows shifted right);
//'   otherwise use even-r layout (even rows shifted right).
//'
//' @return A list of integer vectors of class \code{nb}. Indices are 1-based
//'   and attributes \code{region.id} and \code{topology} ("hex-oddr" or
//'   "hex-evenr") are set. Attribute \code{queen} is set to FALSE for
//'   compatibility with existing code expecting a boolean.
//'
//' @examples
//' \dontrun{
//' nb_hex <- grid_nb_hex_omp(256, 512, oddr = TRUE)
//' }
// [[Rcpp::export]]
List grid_nb_hex_omp(const int nrow,
                     const int ncol,
                     const bool oddr = true)
{
    const int n = nrow * ncol;
    std::vector<std::vector<int>> nb(n);

#pragma omp parallel for schedule(static)
    for (int idx = 0; idx < n; ++idx)
    {
        int r = idx / ncol, c = idx % ncol; // 0-based
        std::vector<int> local_nb;
        local_nb.reserve(6);

        // Offsets depend on row parity and layout
        const bool is_odd = (r % 2) == 1;

        // Horizontal neighbours (W, E)
        int rr = r, cc = c - 1;
        if (cc >= 0) local_nb.push_back(rr * ncol + cc + 1);
        cc = c + 1;
        if (cc < ncol) local_nb.push_back(rr * ncol + cc + 1);

        // Diagonals depend on odd/even row and chosen layout
        if (oddr) {
            // odd-r: odd rows shifted to the right
            // NE, NW (r-1,*), SE, SW (r+1,*)
            if (r - 1 >= 0) {
                // NE: (r-1, c+1) when odd, (r-1, c) when even
                cc = c + (is_odd ? 1 : 0);
                if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                // NW: (r-1, c) when odd, (r-1, c-1) when even
                cc = c + (is_odd ? 0 : -1);
                if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
            }
            if (r + 1 < nrow) {
                // SE: (r+1, c+1) when odd, (r+1, c) when even
                cc = c + (is_odd ? 1 : 0);
                if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                // SW: (r+1, c) when odd, (r+1, c-1) when even
                cc = c + (is_odd ? 0 : -1);
                if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
            }
        } else {
            // even-r: even rows shifted to the right (mirror of odd-r)
            if (r - 1 >= 0) {
                // NE: (r-1, c+1) when even, (r-1, c) when odd
                cc = c + (is_odd ? 0 : 1);
                if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                // NW: (r-1, c) when even, (r-1, c-1) when odd
                cc = c + (is_odd ? -1 : 0);
                if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
            }
            if (r + 1 < nrow) {
                // SE: (r+1, c+1) when even, (r+1, c) when odd
                cc = c + (is_odd ? 0 : 1);
                if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                // SW: (r+1, c) when even, (r+1, c-1) when odd
                cc = c + (is_odd ? -1 : 0);
                if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
            }
        }

        nb[idx].swap(local_nb);
    }

    List out(n);
    for (int i = 0; i < n; ++i)
        out[i] = IntegerVector(nb[i].begin(), nb[i].end());
    out.attr("class") = "nb";
    out.attr("region.id") = Rcpp::seq(1, n);
    out.attr("queen") = false;
    out.attr("topology") = oddr ? "hex-oddr" : "hex-evenr";
    return out;
}

/*------------------------------------------------------------------*/
/* 1c. Parallel generation of Hex adjacency (odd-q / even-q)       */
/*------------------------------------------------------------------*/
//' @title Build Hexagonal (6-neighbour) Neighbour List with column-offset (odd-q/even-q)
//'
//' @description
//'   Convenience wrapper that reuses the row-offset builder on a transposed
//'   lattice to implement pointy-top hex coordinates with column offsets
//'   (odd-q/even-q). Region indices are 1..(nrow*ncol) in column-major order
//'   consistent with mapping id = (gx-1)*nrow + gy.
//'
//' @param nrow Integer. Number of rows in the logical lattice.
//' @param ncol Integer. Number of columns in the logical lattice.
//' @param oddq Logical. If TRUE, odd columns are shifted; otherwise even.
//'
//' @return An \code{nb} list with attributes \code{region.id}, \code{queen}=FALSE,
//'   and \code{topology} set to \code{"hex-oddq"} or \code{"hex-evenq"}.
// [[Rcpp::export]]
List grid_nb_hexq_omp(const int nrow,
                      const int ncol,
                      const bool oddq = true)
{
    // Build on transposed lattice using row-offset builder
    List nb_t = grid_nb_hex_omp(ncol, nrow, oddq);

    // nb_t already has class and neighbours; region.id remains 1..(ncol*nrow)
    // which matches column-major order for the original (nrow x ncol) lattice.
    nb_t.attr("queen") = false;
    nb_t.attr("topology") = oddq ? "hex-oddq" : "hex-evenq";
    return nb_t;
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

    /* 2.1 Single-thread deep copy into C++ containers — fully isolate from R memory */
    std::vector<std::vector<int>> nb_cpp(n);
    for (int i = 0; i < n; ++i)
    {
        IntegerVector v = nb[i];
        nb_cpp[i].assign(v.begin(), v.end());
    }

    /* 2.2 Thread-private triplet buffers */
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

    /* 2.3 Merge triplets */
    size_t total = 0;
    for (auto &v : trip_tls)
        total += v.size();
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(total);
    for (auto &v : trip_tls)
    {
        std::move(v.begin(), v.end(), std::back_inserter(trip));
        std::vector<Eigen::Triplet<double>>().swap(v); // release thread-local buffer
    }

    /* 2.4 Construct compressed sparse matrix */
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> W(n, n);
    W.setFromTriplets(trip.begin(), trip.end());
    W.makeCompressed();

    return Rcpp::wrap(W); // dgCMatrix
}

// Compatibility alias that calls the original function
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
