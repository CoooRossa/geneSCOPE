\name{ScoreModuleQuality}
\alias{ScoreModuleQuality}
\title{Score module quality using STRINGdb}
\description{
Computes module-level coherence against STRINGdb by comparing within-module
versus between-module STRING combined scores (and a delta statistic). This is
the only supported module-quality evaluation backend in geneSCOPE.
}
\usage{
ScoreModuleQuality(
  x,
  membership = NULL,
  membership_cols = NULL,
  grid_name = "grid30",
  stats_layer = "LeeStats_Xz",
  graph_slot = "g_consensus",
  edge_weight = "weight",
  species = 9606,
  string_score_threshold = 700,
  n_perm = 200,
  seed = 1,
  writeback = FALSE,
  verbose = getOption("geneSCOPE.verbose", TRUE),
  ...
)
}
\arguments{
\item{x}{A \code{scope_object}, \code{igraph}, or list with \code{graph} and optional \code{membership}.}
\item{membership}{Optional membership vector named by gene.}
\item{membership_cols}{Optional character vector of \code{x@meta.data} column names to score. When provided, \code{membership} is ignored and each column is scored separately (results include \code{membership_col}).}
\item{grid_name}{Grid layer name for \code{scope_object} input.}
\item{stats_layer}{Stats layer name for \code{scope_object} input.}
\item{graph_slot}{Graph slot name for \code{scope_object} input.}
\item{edge_weight}{Edge attribute used as weight.}
\item{species}{NCBI taxonomy ID for STRING.}
\item{string_score_threshold}{STRING combined score threshold for positives.}
\item{n_perm}{Number of label permutations for the within/between delta null.}
\item{seed}{Random seed for permutations.}
\item{writeback}{When TRUE and \code{x} is a \code{scope_object}, store results under \code{scope_obj@stats[[grid_name]][[stats_layer]]$module_quality$STRINGdb}.}
\item{verbose}{Emit progress messages when TRUE.}
\item{...}{Reserved for future extensions.}
}
\value{
A list with \code{summary}, \code{per_module}, and \code{meta}. When
\code{writeback = TRUE} and \code{x} is a \code{scope_object}, the returned list
also includes \code{scope_obj} (updated) and \code{run_id} (storage key).
}
\details{
Per-module results include within/between STRING score summaries and
E_avg-style diagnostics computed from STRING interactions (mean combined score
over all within-module pairs, treating missing pairs as zero).
}
\examples{
\dontrun{
res <- ScoreModuleQuality(scope_obj, grid_name = "grid30", n_perm = 50)
}
}
