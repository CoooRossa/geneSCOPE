\name{GetModuleQualityScores}
\alias{GetModuleQualityScores}
\title{Extract stored STRINGdb module-quality scores}
\description{
Reads module-quality payloads stored by \code{ScoreModuleQuality(writeback = TRUE)}
under \code{scope_obj@stats} and returns a flattened per-module table.
}
\usage{
GetModuleQualityScores(x, grid_name = NULL, stats_layer = NULL)
}
\arguments{
\item{x}{A \code{scope_object}.}
\item{grid_name}{Optional grid filter.}
\item{stats_layer}{Optional stats layer filter.}
}
\value{
A data.frame with one row per module per stored run (empty if none found).
Includes \code{grid_name}, \code{stats_layer}, \code{run_id}, and (when
available) \code{membership_col}.
}
\examples{
\dontrun{
res <- GetModuleQualityScores(scope_obj, grid_name = "grid30")
}
}
