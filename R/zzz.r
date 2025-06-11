.onLoad <- function(libname, pkgname) {

  Sys.setenv(OPENBLAS_NUM_THREADS = "1")

  n <- getOption("mypkg.omp_threads", NA_integer_)
  if (is.na(n) || n <= 0) {
    n <- parallel::detectCores(logical = FALSE)   
    if (is.na(n) || n < 1) n <- parallel::detectCores()  
    if (is.na(n) || n < 1) n <- 1               
  }

  Sys.setenv(OMP_NUM_THREADS = as.character(n))

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try({
      RhpcBLASctl::blas_set_num_threads(1)  
      RhpcBLASctl::omp_set_num_threads(n)   
    }, silent = TRUE)
  }

  packageStartupMessage(
    sprintf("[FG2CLI] OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=%d",
        as.integer(Sys.getenv("OMP_NUM_THREADS", 1L)), "You can set options(mypkg.omp_threads = x) to change the number of threads used by OpenMP in this package (x is the number of threads you want to use).")
  )
}