#' @title External Reference Optional Module
#' @description Adapter layer for external reference databases (STRINGdb).
#'   This module requires STRINGdb package for protein-protein interaction validation.
#' @name optional-external-reference
#' @keywords internal
NULL

#' Check if STRINGdb is available
#' @return Logical indicating if STRINGdb is installed
#' @keywords internal
.stringdb_available <- function() {
  requireNamespace("STRINGdb", quietly = TRUE)
}

#' Require STRINGdb or stop with informative error
#' @return Invisible TRUE if available
#' @keywords internal
.require_stringdb_or_stop <- function() {
  if (identical(Sys.getenv("GENESCOPE_DISABLE_STRINGDB"), "1")) {
    stop(
      "STRINGdb is disabled via GENESCOPE_DISABLE_STRINGDB=1.\n",
      "Remove this environment variable to enable STRINGdb functionality."
    )
  }
  if (!.stringdb_available()) {
    stop(
      "External reference evaluation requires STRINGdb package.\n",
      "Install with:\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install('STRINGdb')"
    )
  }
  invisible(TRUE)
}

#' Create STRINGdb instance (adapter)
#' @description Wrapper around STRINGdb::STRINGdb with dependency checking
#' @param species NCBI taxonomy ID (default 9606 for human)
#' @param version STRING version (default "11.5")
#' @param score_threshold Minimum score threshold (default 0)
#' @param input_directory Optional cache directory
#' @return STRINGdb object
#' @keywords internal
.connect_stringdb <- function(species = 9606, version = "11.5", 
                                score_threshold = 0, input_directory = NULL) {
  .require_stringdb_or_stop()
  
  if (is.null(input_directory)) {
    STRINGdb::STRINGdb$new(
      version = version,
      species = species,
      score_threshold = score_threshold
    )
  } else {
    STRINGdb::STRINGdb$new(
      version = version,
      species = species,
      score_threshold = score_threshold,
      input_directory = input_directory
    )
  }
}

#' Map genes to STRING IDs (adapter)
#' @description Map gene symbols to STRING protein IDs
#' @param string_db STRINGdb object
#' @param genes Character vector of gene symbols
#' @return Named vector mapping genes to STRING IDs
#' @keywords internal
.map_genes_to_string <- function(string_db, genes) {
  .require_stringdb_or_stop()
  
  genes <- unique(as.character(genes))
  if (!length(genes)) {
    return(setNames(character(0), character(0)))
  }
  
  df <- data.frame(gene = genes, stringsAsFactors = FALSE)
  mapped <- string_db$map(df, "gene", removeUnmappedRows = FALSE)
  setNames(mapped$STRING_id, mapped$gene)
}

#' Get STRING interactions (adapter)
#' @description Retrieve protein-protein interactions from STRING
#' @param string_db STRINGdb object
#' @param string_ids Character vector of STRING protein IDs
#' @return Data frame with interactions (from, to, combined_score)
#' @keywords internal
.get_string_interactions <- function(string_db, string_ids) {
  .require_stringdb_or_stop()
  
  ids <- unique(na.omit(as.character(string_ids)))
  if (!length(ids)) {
    return(data.frame(from = character(0), to = character(0), 
                      combined_score = numeric(0)))
  }
  
  interactions <- string_db$get_interactions(ids)
  if (!is.data.frame(interactions) || !nrow(interactions)) {
    data.frame(from = character(0), to = character(0), combined_score = numeric(0))
  } else {
    interactions[, c("from", "to", "combined_score"), drop = FALSE]
  }
}
