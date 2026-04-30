#' Validate STRINGdb availability.
#' @keywords internal
.stringdb_require_or_stop <- function() {
    if (identical(Sys.getenv("GENESCOPE_DISABLE_STRINGDB"), "1")) {
        stop(
            "STRINGdb is required for STRING evaluation. Install with:\n",
            "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
            "  BiocManager::install('STRINGdb')"
        )
    }
    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
        stop(
            "STRINGdb is required for STRING evaluation. Install with:\n",
            "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
            "  BiocManager::install('STRINGdb')"
        )
    }
    invisible(TRUE)
}

#' Read update flag from env/option.
#' @keywords internal
.stringdb_update_flag <- function() {
    opt <- getOption("geneSCOPE.stringdb_update", FALSE)
    if (isTRUE(opt)) return(TRUE)
    env <- Sys.getenv("GENESCOPE_STRINGDB_UPDATE", "")
    if (nzchar(env) && !identical(env, "0")) return(TRUE)
    FALSE
}

#' Download STRINGdb file with md5 tracking.
#' @keywords internal
.stringdb_download_with_md5 <- function(url, dest, force_update = FALSE) {
    dest_dir <- dirname(dest)
    if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
    md5_path <- paste0(dest, ".md5")

    if (!force_update && file.exists(dest) && file.info(dest)$size > 0) {
        if (file.exists(md5_path)) return(invisible(TRUE))
        current <- tryCatch(tools::md5sum(dest), error = function(e) NA_character_)
        current <- unname(current)
        if (!is.na(current)) writeLines(current, md5_path)
        return(invisible(TRUE))
    }

    use_tmp <- isTRUE(force_update) && file.exists(dest)
    tmp <- if (use_tmp) paste0(dest, ".tmp") else dest
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    if (!file.exists(tmp) || file.info(tmp)$size == 0) {
        if (use_tmp && file.exists(tmp)) unlink(tmp)
        stop("Failed to download STRINGdb data: ", basename(dest))
    }
    new_md5 <- unname(tryCatch(tools::md5sum(tmp), error = function(e) NA_character_))

    if (use_tmp && file.exists(dest)) {
        old_md5 <- NA_character_
        if (file.exists(md5_path)) {
            old_md5 <- readLines(md5_path, warn = FALSE)
            old_md5 <- if (length(old_md5)) trimws(old_md5[1]) else NA_character_
        } else {
            old_md5 <- unname(tryCatch(tools::md5sum(dest), error = function(e) NA_character_))
        }
        if (!is.na(old_md5) && !is.na(new_md5) && nzchar(old_md5) && identical(old_md5, new_md5)) {
            unlink(tmp)
            if (!file.exists(md5_path) && !is.na(new_md5)) writeLines(new_md5, md5_path)
            return(invisible(TRUE))
        }
        file.rename(tmp, dest)
    }
    if (!is.na(new_md5)) writeLines(new_md5, md5_path)
    invisible(TRUE)
}

#' Prepare STRINGdb download cache with md5 checks.
#' @keywords internal
.stringdb_prepare_cache_files <- function(string_db, cache_dir, force_update = NULL) {
    if (is.null(cache_dir) || !nzchar(cache_dir)) return(invisible(FALSE))
    if (is.null(force_update)) force_update <- .stringdb_update_flag()

    protocol <- string_db$protocol
    species <- string_db$species
    file_version <- string_db$file_version

    aliases_base <- paste0(species, ".protein.aliases.v", file_version, ".txt.gz")
    info_base <- paste0(species, ".protein.info.v", file_version, ".txt.gz")

    aliases_url <- paste0(
        protocol,
        "://stringdb-downloads.org/download/protein.aliases.v",
        file_version,
        "/",
        aliases_base
    )
    info_url <- paste0(
        protocol,
        "://stringdb-downloads.org/download/protein.info.v",
        file_version,
        "/",
        info_base
    )

    network_type_param <- ""
    if (tolower(string_db$network_type) == "physical") {
        network_type_param <- "physical."
    }
    link_data_param <- "links.v"
    link_data <- tolower(string_db$link_data)
    if (link_data == "detailed") {
        link_data_param <- "links.detailed.v"
    } else if (link_data == "full") {
        link_data_param <- "links.full.v"
    }
    links_base <- paste0(
        species,
        ".protein.",
        network_type_param,
        link_data_param,
        file_version,
        ".txt.gz"
    )
    links_url <- paste0(
        protocol,
        "://stringdb-downloads.org/download/protein.",
        network_type_param,
        link_data_param,
        file_version,
        "/",
        links_base
    )

    .stringdb_download_with_md5(aliases_url, file.path(cache_dir, aliases_base), force_update = force_update)
    .stringdb_download_with_md5(info_url, file.path(cache_dir, info_base), force_update = force_update)
    .stringdb_download_with_md5(links_url, file.path(cache_dir, links_base), force_update = force_update)
    invisible(TRUE)
}

.stringdb_subscore_columns <- function() {
    c("nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore")
}

.stringdb_detailed_links_info <- function(string_db, cache_dir = NULL) {
    network_type_param <- ""
    if (tolower(string_db$network_type) == "physical") {
        network_type_param <- "physical."
    }
    file_version <- string_db$file_version
    species <- string_db$species
    file_base <- paste0(species, ".protein.", network_type_param, "links.detailed.v", file_version, ".txt.gz")
    url <- paste0(
        string_db$protocol,
        "://stringdb-downloads.org/download/protein.",
        network_type_param,
        "links.detailed.v",
        file_version,
        "/",
        file_base
    )
    base_dir <- cache_dir
    if (is.null(base_dir) || !nzchar(base_dir)) {
        base_dir <- tryCatch(as.character(string_db$input_directory), error = function(e) "")
    }
    if (is.null(base_dir) || !nzchar(base_dir)) {
        base_dir <- tempdir()
    }
    list(url = url, file_base = file_base, file_path = file.path(base_dir, file_base))
}

.stringdb_read_detailed_subset <- function(file_path, ids, cols) {
    ids <- unique(na.omit(as.character(ids)))
    if (!length(ids)) return(data.table::data.table())

    zcat_cmd <- Sys.which("gzcat")
    if (!nzchar(zcat_cmd)) zcat_cmd <- Sys.which("zcat")
    gzip_cmd <- Sys.which("gzip")
    if (!nzchar(zcat_cmd) && nzchar(gzip_cmd)) zcat_cmd <- gzip_cmd
    awk_cmd <- Sys.which("awk")

    dt <- data.table::data.table()
    if (nzchar(zcat_cmd) && nzchar(awk_cmd)) {
        ids_file <- tempfile(fileext = ".txt")
        on.exit(unlink(ids_file), add = TRUE)
        writeLines(ids, ids_file, useBytes = TRUE)
        decompress <- if (basename(zcat_cmd) == "gzip") {
            paste(zcat_cmd, "-dc", shQuote(file_path))
        } else {
            paste(zcat_cmd, shQuote(file_path))
        }
        cmd <- paste(
            decompress,
            "|",
            awk_cmd,
            shQuote("NR==FNR{a[$1]=1;next} FNR==1{print;next} ($1 in a) && ($2 in a){print}"),
            shQuote(ids_file),
            "-"
        )
        dt <- tryCatch(data.table::fread(cmd = cmd, showProgress = FALSE), error = function(e) data.table::data.table())
    } else {
        # Avoid attempting to fread() the full links.detailed file when streaming tools are missing.
        return(data.table::data.table())
    }
    if (!nrow(dt)) return(data.table::data.table())
    if (!all(cols %in% names(dt))) return(data.table::data.table())
    dt[, ..cols]
}

.stringdb_append_subscores <- function(edge_df, string_db, cache_dir = NULL, force_update = NULL) {
    if (!is.data.frame(edge_df) || !nrow(edge_df)) {
        if (!is.data.frame(edge_df)) edge_df <- as.data.frame(edge_df)
        for (col in .stringdb_subscore_columns()) edge_df[[col]] <- numeric(0)
        return(edge_df)
    }
    if (!all(c("string_from", "string_to") %in% names(edge_df))) {
        for (col in .stringdb_subscore_columns()) edge_df[[col]] <- NA_real_
        return(edge_df)
    }

    ids <- unique(na.omit(c(as.character(edge_df$string_from), as.character(edge_df$string_to))))
    ids <- ids[nzchar(ids)]
    if (!length(ids)) {
        for (col in .stringdb_subscore_columns()) edge_df[[col]] <- NA_real_
        return(edge_df)
    }

    links <- .stringdb_detailed_links_info(string_db, cache_dir = cache_dir)
    dir.create(dirname(links$file_path), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(links$file_path) || file.info(links$file_path)$size == 0) {
        if (is.null(force_update)) force_update <- .stringdb_update_flag()
        .stringdb_download_with_md5(links$url, links$file_path, force_update = isTRUE(force_update))
    }
    if (!file.exists(links$file_path) || file.info(links$file_path)$size == 0) {
        for (col in .stringdb_subscore_columns()) edge_df[[col]] <- NA_real_
        return(edge_df)
    }

    cols <- c(
        "protein1",
        "protein2",
        "neighborhood",
        "fusion",
        "cooccurence",
        "coexpression",
        "experimental",
        "database",
        "textmining",
        "combined_score"
    )
    interactions <- .stringdb_read_detailed_subset(links$file_path, ids, cols)
    if (!nrow(interactions)) {
        for (col in .stringdb_subscore_columns()) edge_df[[col]] <- NA_real_
        return(edge_df)
    }

    interactions <- data.table::as.data.table(interactions)
    interactions <- interactions[protein1 %chin% ids & protein2 %chin% ids]
    if (!nrow(interactions)) {
        for (col in .stringdb_subscore_columns()) edge_df[[col]] <- NA_real_
        return(edge_df)
    }
    interactions[, key := paste(pmin(protein1, protein2), pmax(protein1, protein2), sep = "|")]
    interactions <- unique(interactions, by = "key")
    interactions[, `:=`(
        nscore = as.numeric(neighborhood),
        fscore = as.numeric(fusion),
        pscore = as.numeric(cooccurence),
        ascore = as.numeric(coexpression),
        escore = as.numeric(experimental),
        dscore = as.numeric(database),
        tscore = as.numeric(textmining)
    )]
    sub <- interactions[, c("key", .stringdb_subscore_columns()), with = FALSE]

    dt <- data.table::as.data.table(edge_df)
    dt[, .row := .I]
    dt[, key := paste(pmin(string_from, string_to), pmax(string_from, string_to), sep = "|")]
    dt <- merge(dt, sub, by = "key", all.x = TRUE, sort = FALSE)
    data.table::setorder(dt, .row)
    dt[, c("key", ".row") := NULL]
    as.data.frame(dt)
}

#' Create STRINGdb instance.
#' @keywords internal
.stringdb_connect <- function(species = 9606,
                              score_threshold = 0,
                              version = "11.5",
                              cache_dir = NULL) {
    .stringdb_require_or_stop()
    cache_dir_use <- cache_dir
    if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
        cache_dir_use <- Sys.getenv("GENESCOPE_STRINGDB_DIR", "")
    }
    if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
        cache_dir_use <- getOption("geneSCOPE.stringdb_dir", "")
    }
    if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
        cache_dir_use <- tryCatch(
            tools::R_user_dir("geneSCOPE", "cache"),
            error = function(e) ""
        )
    }
    if (is.null(cache_dir_use) || !nzchar(cache_dir_use)) {
        return(STRINGdb::STRINGdb$new(
            version = version,
            species = species,
            score_threshold = score_threshold
        ))
    }
    if (!dir.exists(cache_dir_use)) {
        ok <- tryCatch({
            dir.create(cache_dir_use, recursive = TRUE)
            TRUE
        }, error = function(e) FALSE)
        if (!isTRUE(ok)) {
            return(STRINGdb::STRINGdb$new(
                version = version,
                species = species,
                score_threshold = score_threshold
            ))
        }
    }
    string_db <- STRINGdb::STRINGdb$new(
        version = version,
        species = species,
        score_threshold = score_threshold,
        input_directory = cache_dir_use
    )
    .stringdb_prepare_cache_files(string_db, cache_dir_use)
    string_db
}

#' Map genes to STRING IDs with optional caching.
#' @keywords internal
.stringdb_map_genes <- function(string_db, genes, cache_dir = NULL) {
    genes <- unique(as.character(genes))
    if (!length(genes)) {
        return(setNames(character(0), character(0)))
    }
    cache_path <- NULL
    if (!is.null(cache_dir) && requireNamespace("digest", quietly = TRUE)) {
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
        key <- digest::digest(list(genes = sort(genes), species = string_db$species))
        cache_path <- file.path(cache_dir, paste0("string_mapping_", key, ".rds"))
        if (file.exists(cache_path)) {
            return(readRDS(cache_path))
        }
    }
    df <- data.frame(gene = genes, stringsAsFactors = FALSE)
    mapped <- string_db$map(df, "gene", removeUnmappedRows = FALSE)
    mapping <- setNames(mapped$STRING_id, mapped$gene)
    if (!is.null(cache_path)) saveRDS(mapping, cache_path)
    mapping
}

#' Fetch STRING interactions for provided STRING IDs with optional caching.
#' @keywords internal
.stringdb_get_interactions <- function(string_db, string_ids, cache_dir = NULL) {
    ids <- unique(na.omit(as.character(string_ids)))
    if (!length(ids)) {
        return(data.frame(from = character(0), to = character(0), combined_score = numeric(0)))
    }
    cache_path <- NULL
    if (!is.null(cache_dir) && requireNamespace("digest", quietly = TRUE)) {
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
        key <- digest::digest(list(ids = sort(ids), species = string_db$species))
        cache_path <- file.path(cache_dir, paste0("string_interactions_", key, ".rds"))
        if (file.exists(cache_path)) {
            return(readRDS(cache_path))
        }
    }
    interactions <- string_db$get_interactions(ids)
    if (!is.data.frame(interactions) || !nrow(interactions)) {
        interactions <- data.frame(from = character(0), to = character(0), combined_score = numeric(0))
    } else {
        interactions <- interactions[, c("from", "to", "combined_score"), drop = FALSE]
    }
    if (!is.null(cache_path)) saveRDS(interactions, cache_path)
    interactions
}

#' Annotate edges with STRING scores and labels.
#' @keywords internal
.stringdb_label_edges <- function(edge_df, gene_to_string, interactions, score_threshold = 700) {
    edge_df$string_from <- gene_to_string[edge_df$from]
    edge_df$string_to <- gene_to_string[edge_df$to]
    comparable <- !is.na(edge_df$string_from) & !is.na(edge_df$string_to)

    if (nrow(interactions)) {
        key_int <- paste(
            pmin(interactions$from, interactions$to),
            pmax(interactions$from, interactions$to),
            sep = "|"
        )
        score_map <- setNames(interactions$combined_score, key_int)
    } else {
        score_map <- setNames(numeric(0), character(0))
    }

    key_edges <- paste(
        pmin(edge_df$string_from, edge_df$string_to),
        pmax(edge_df$string_from, edge_df$string_to),
        sep = "|"
    )
    string_score <- score_map[key_edges]
    string_score[is.na(string_score)] <- 0
    string_score[!comparable] <- NA_real_

    edge_df$string_score <- as.numeric(string_score)
    edge_df$label <- ifelse(!is.na(edge_df$string_score), edge_df$string_score >= score_threshold, NA)
    edge_df
}

#' Compute AUPRC for a score ranking.
#' @keywords internal
.compute_auprc <- function(labels, scores) {
    ok <- !is.na(labels) & !is.na(scores)
    labels <- labels[ok]
    scores <- scores[ok]
    if (!length(labels)) return(NA_real_)
    labels <- as.integer(labels) > 0
    n_pos <- sum(labels)
    n_neg <- sum(!labels)
    if (n_pos == 0 || n_neg == 0) return(NA_real_)
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    tp <- cumsum(labels)
    fp <- cumsum(!labels)
    recall <- tp / n_pos
    precision <- tp / pmax(tp + fp, 1)
    recall <- c(0, recall)
    precision <- c(1, precision)
    sum(diff(recall) * precision[-1])
}

#' Compute AUROC for a score ranking.
#' @keywords internal
.compute_auroc <- function(labels, scores) {
    ok <- !is.na(labels) & !is.na(scores)
    labels <- labels[ok]
    scores <- scores[ok]
    if (!length(labels)) return(NA_real_)
    labels <- as.integer(labels) > 0
    n_pos <- sum(labels)
    n_neg <- sum(!labels)
    if (n_pos == 0 || n_neg == 0) return(NA_real_)
    ranks <- rank(scores, ties.method = "average")
    sum_ranks_pos <- sum(ranks[labels])
    (sum_ranks_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' Precision at K for a score ranking.
#' @keywords internal
.precision_at_k <- function(labels, scores, k = 100L) {
    ok <- !is.na(labels) & !is.na(scores)
    labels <- labels[ok]
    scores <- scores[ok]
    if (!length(labels)) return(NA_real_)
    ord <- order(scores, decreasing = TRUE)
    labels <- labels[ord]
    k <- min(length(labels), as.integer(k))
    if (k <= 0L) return(NA_real_)
    mean(as.integer(labels[seq_len(k)]))
}
