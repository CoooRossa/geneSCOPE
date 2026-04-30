\name{SummarizeModuleQualityScores}
\alias{SummarizeModuleQualityScores}
\title{Summarize stored STRINGdb module-quality scores}
\description{
Summarizes module-quality payloads stored by \code{ScoreModuleQuality(writeback = TRUE)}.
}
\usage{
SummarizeModuleQualityScores(x, grid_name = NULL, stats_layer = NULL)
}
\arguments{
\item{x}{A \code{scope_object} containing stored module-quality scores.}
\item{grid_name}{Optional grid filter.}
\item{stats_layer}{Optional stats layer filter.}
}
\value{
A list with:
\itemize{
  \item \code{per_module}: flattened per-module table (see \code{GetModuleQualityScores()}),
  \item \code{per_run}: per-run summary table,
  \item \code{summary_text}: short description.
}
}
\examples{
\dontrun{
summary <- SummarizeModuleQualityScores(scope_obj)
}
}
