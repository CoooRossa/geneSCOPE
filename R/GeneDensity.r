#' @title Compute Grid-Based Density from coordObj (Custom Output Name)
#' @description
#'   Using the existing grid counts in coordObj, compute molecular density per grid cell
#'   based on a specified gene set or cluster information. The result is saved under
#'   `coordObj@grid[[grid_name]][[density_name]]` with a user‐specified field name.
#'
#' @param coordObj      A coordObj that must contain `coordObj@grid[[grid_name]]$grid_info`
#'                      and `coordObj@grid[[grid_name]]$counts`.
#' @param grid_name     Character (optional). The name of the grid sublayer to use (e.g. `"grid_lenGrid50"`).
#'                      If NULL and there is only one sublayer under `coordObj@grid`, that sublayer is used;
#'                      otherwise, the user must specify it explicitly.
#' @param density_name  Character. The name under which to save the density results in
#'                      `coordObj@grid[[grid_name]][[density_name]]`.
#' @param cluster_col   Character (optional). A column name in `coordObj@meta.data` used to select genes by cluster.
#' @param cluster_num   Integer or character (optional). When `cluster_col` is not NULL, select genes where
#'                      `coordObj@meta.data[[cluster_col]] == cluster_num`.
#' @param genes         Character vector (optional). Directly specify the list of genes for which to compute density.
#'                      If both `genes` and `cluster_col` are provided, `genes` takes priority.
#'
#' @return The modified coordObj, where `coordObj@grid[[grid_name]][[density_name]]` is a data.frame
#'         containing:
#'           - grid_id: Grid cell ID
#'           - gx, gy: Grid indices in the x and y directions
#'           - count: Total count of selected genes in that grid cell
#'           - density: count divided by grid cell area (µm²)
#'           - xmin, xmax, ymin, ymax: Grid cell boundary coordinates
#'
#' @importFrom dplyr filter pull group_by summarise ungroup right_join mutate select
#' @export
computeDensity <- function(coordObj,
                                grid_name         = NULL,
                                density_name      = "densityDF",
                                cluster_col       = NULL,
                                cluster_num       = NULL,
                                layer_name        = NULL, 
                                genes             = NULL,
                                normalize_method  = c("none", "per_grid", "global_gene")) {
normalize_method <- match.arg(normalize_method)

## —— 0. —— 
if (is.null(grid_name)) {
 sub_layers <- names(coordObj@grid)
 if (length(sub_layers) == 1) {
   grid_layer_name <- sub_layers[[1]]
 } else {
   stop("coordObj@grid has multiple sub-layers; please specify grid_name explicitly.")
 }
} else {
 grid_layer_name <- as.character(grid_name)
}
if (!(grid_layer_name %in% names(coordObj@grid))) {
 stop("Not found: coordObj@grid[['", grid_layer_name, "']]")
}

## —— 1. —— 
grid_info <- coordObj@grid[[grid_layer_name]]$grid_info
if (is.null(grid_info)) {
 stop("coordObj@grid[['", grid_layer_name, "']] must contain grid_info.")
}
counts_dt <- coordObj@grid[[grid_layer_name]][[ layer_name ]]
if (is.null(counts_dt)) {
 stop("coordObj@grid[['", grid_layer_name, "']] not contain counts_dt for layer_name: ", layer_name)
}
if (!all(c("grid_id","gene","count") %in% colnames(counts_dt))) {
 stop("Layer '", layer_name, "' must contain columns: grid_id, gene, count.")
}

## —— 2.  —— 
if (!is.null(genes)) {
 sel_genes <- genes
} else if (!is.null(cluster_col) && !is.null(cluster_num)) {
 if (is.null(coordObj@meta.data) ||
     !(cluster_col %in% colnames(coordObj@meta.data))) {
   stop("coordObj@meta.data not found, or does not contain column: ", cluster_col)
 }
 sel_genes <- rownames(coordObj@meta.data)[
   coordObj@meta.data[[cluster_col]] == cluster_num
 ]
 if (length(sel_genes) == 0) {
   stop("In meta.data  col '", cluster_col,
        " == ", cluster_num, " found no genes.")
 }
} else {
 stop("Please specify either `genes` or both `cluster_col` and `cluster_num`.")
}

library(dplyr)

## —— 3. —— 
if (normalize_method == "none") {
 # —— 3.1 —— 
 processed_counts <- counts_dt

} else if (normalize_method == "per_grid") {
 # —— 3.2 —— 
 total_per_grid <- counts_dt %>%
   group_by(grid_id) %>%
   summarise(total_count_in_grid = sum(count), .groups = "drop")

 processed_counts <- counts_dt %>%
   left_join(total_per_grid, by = "grid_id") %>%
   mutate(count = ifelse(total_count_in_grid == 0, 0, count / total_count_in_grid)) %>%
   select(grid_id, gene, count) 

} else if (normalize_method == "global_gene") {
 # —— 3.3 —— 
 global_per_gene <- counts_dt %>%
   group_by(gene) %>%
   summarise(global_count = sum(count), .groups = "drop")

 processed_counts <- counts_dt %>%
   left_join(global_per_gene, by = "gene") %>%
   mutate(count = ifelse(global_count == 0, 0, count / global_count)) %>%
   select(grid_id, gene, count)

}

## —— 4. —— 
sel_counts <- processed_counts %>%
 filter(gene %in% sel_genes) %>%
 group_by(grid_id) %>%
 summarise(count = sum(count), .groups = "drop")

## —— 5. —— 
grid_with_area <- grid_info %>%
 mutate(
   width  = xmax - xmin,
   height = ymax - ymin,
   area   = width * height
 ) %>%
 select(grid_id, area)

## —— 6. —— 
density_tbl <- grid_with_area %>%
 left_join(sel_counts, by = "grid_id") %>%
 mutate(count = ifelse(is.na(count), 0, count),
        density = count / area) %>%
 select(grid_id, density)

## —— 7. —— 
if (is.null(coordObj@grid[[grid_layer_name]]$densityDF)) {
 all_ids <- grid_info$grid_id
 empty_df <- data.frame(row.names = all_ids)
 coordObj@grid[[grid_layer_name]]$densityDF <- empty_df
}
densityDF <- coordObj@grid[[grid_layer_name]]$densityDF

dens_vec <- density_tbl$density
names(dens_vec) <- density_tbl$grid_id

densityDF[[density_name]] <- dens_vec[rownames(densityDF)]
densityDF[[density_name]][is.na(densityDF[[density_name]])] <- 0

coordObj@grid[[grid_layer_name]]$densityDF <- densityDF
return(coordObj)
}
