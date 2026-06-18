#' Plot a psychHypergraph Object
#'
#' Visualize psychometric hypergraphs using the \pkg{HyperG} package.
#' Hyperedges can be selected from the complete hypergraph or from the
#' redundant/synergistic sub-hypergraphs, and may be filtered according
#' to their O-information values.
#'
#' Hyperedge colors are automatically mapped to O-information values using
#' a red-white-blue gradient:
#'
#' \itemize{
#'   \item Red = synergistic hyperedges (negative O-information)
#'   \item White = weak higher-order interaction
#'   \item Blue = redundant hyperedges (positive O-information)
#' }
#'
#' When node groups are provided, node colors can be automatically assigned
#' according to group membership.
#'
#' @param x A \code{psychHypergraph} object produced by
#'   \code{\link{build_hypergraph}}.
#'
#' @param type Character string specifying which hypergraph to visualize:
#'   \describe{
#'     \item{"full"}{All significant hyperedges.}
#'     \item{"redundant"}{Only hyperedges with positive O-information.}
#'     \item{"synergistic"}{Only hyperedges with negative O-information.}
#'   }
#'
#' @param top_n Optional integer specifying the number of strongest
#'   hyperedges (based on absolute O-information) to display.
#'   If \code{NULL}, all available hyperedges are plotted.
#'
#' @param oinfo_threshold Optional minimum absolute O-information value.
#'   Hyperedges with |O| < threshold are excluded.
#'
#' @param groups Optional grouping vector indicating node membership.
#'   Must have length equal to the number of variables in the hypergraph.
#'   Used for automatic node coloring.
#'
#' @param labels Optional character vector of node labels.
#'   Must have length equal to the number of variables.
#'   If \code{NULL}, original variable names are used.
#'
#' @param vertex.color Node colors.
#'
#' Can be:
#' \itemize{
#'   \item a single color applied to all nodes;
#'   \item a vector of length equal to the number of variables;
#'   \item when \code{groups} is supplied, a vector of length equal to the
#'         number of groups.
#' }
#'
#' If \code{NULL}, colors are automatically generated from group membership
#' (if available) or set to white.
#'
#' @param vertex.size Numeric scalar or vector specifying node sizes.
#'
#' @param vertex.label.cex Numeric scaling factor for node labels.
#'
#' @param vertex.label.color Color of node labels.
#'
#' @param vertex.shape Shape of nodes passed to
#'   \code{HyperG::plot.hypergraph()}.
#'
#' @param edge.color Color of graph edges.
#'
#' @param edge.width Width of hyperedge boundaries.
#'
#' @param edge.lty Line type of hyperedge boundaries.
#'
#' @param edge.alpha Transparency level used for hyperedge regions.
#'   Values range from 0 (fully transparent) to 1 (fully opaque).
#'
#' @param layout Optional layout matrix.
#'   If \code{NULL}, a spring layout is generated automatically from the
#'   projected adjacency matrix.
#'
#' @param main Optional plot title.
#'   If \code{NULL}, a title is generated automatically.
#'
#' @param seed Random seed used for layout generation.
#'
#' @param ... Additional arguments passed to
#'   \code{HyperG::plot.hypergraph()}.
#'
#' @details
#'
#' Hyperedges are represented as shaded regions surrounding their member
#' nodes. The color of each hyperedge reflects the sign and magnitude of
#' its O-information value.
#'
#' The plotting procedure consists of:
#'
#' \enumerate{
#'   \item Selecting a hypergraph (\code{full}, \code{redundant},
#'         or \code{synergistic});
#'   \item Filtering hyperedges using \code{top_n} and/or
#'         \code{oinfo_threshold};
#'   \item Constructing a HyperG hypergraph object;
#'   \item Computing a force-directed layout;
#'   \item Rendering nodes and hyperedges.
#' }
#'
#' Large hypergraphs may be difficult to interpret visually.
#' Using \code{top_n} or \code{oinfo_threshold} is recommended when the
#' number of hyperedges is large.
#'
#' @return
#' A list containing:
#' \item{hypergraph}{The HyperG hypergraph object used for plotting.}
#' \item{layout}{The layout matrix used for node positions.}
#' \item{incidence_matrix}{The incidence matrix of the plotted hypergraph.}
#' \item{oinfo}{The O-information values corresponding to the plotted hyperedges.}
#'
#' @seealso
#' \code{\link{build_hypergraph}}
#'
#' @examples
#' \dontrun{
#' library(psych)
#' data(bfi)
#'
#' hg <- build_hypergraph(
#'   data = bfi[,1:25],
#'   k = 3:5,
#'   candidate_method = "all",
#'   cor_method = "pearson",
#'   parallel = TRUE
#' )
#'
#' # Full hypergraph
#' plot(hg)
#'
#' # Only synergistic hyperedges
#' plot(hg, type = "synergistic")
#'
#' # Top 20 strongest hyperedges
#' plot(hg, top_n = 20)
#'
#' # Threshold-based filtering
#' plot(hg, oinfo_threshold = 0.05)
#'
#' # Group coloring
#' groups <- rep(c("A","C","E","N","O"), each = 5)
#' plot(hg, groups = groups)
#'
#' }
#'
#' @export
plot.psychHypergraph <- function(
    x,
    type = c("full", "redundant", "synergistic"),
    top_n = NULL,
    oinfo_threshold = NULL,
    groups = NULL,
    labels = NULL,
    vertex.color = NULL,
    vertex.size = 15,
    vertex.label.cex = 1,
    vertex.label.color = "black",
    vertex.shape = "circle",
    edge.color = NA,
    edge.width = 1,
    edge.lty = 1,
    edge.alpha = 0.2,
    layout = NULL,
    main = NULL,
    seed = 123,
    ...
) {
  set.seed(seed)

  if (!requireNamespace("HyperG", quietly = TRUE)) {
    stop("The 'HyperG' package is required for plotting. Please install it.")
  }

  type <- match.arg(type)

  incidence_mat <- switch(
    type,
    full = x$incidence_matrix,
    redundant = x$redundant_incidence_matrix,
    synergistic = x$synergistic_incidence_matrix
  )

  if (is.null(incidence_mat) || nrow(incidence_mat) == 0) {
    warning("No hyperedges to plot for type = '", type, "'")
    return(invisible(NULL))
  }

  oinfo_df <- x$all_test_results
  oinfo_df$abs_oinfo <- abs(oinfo_df$Oinfo)

  hyperedge_oinfo <- rep(NA_real_, nrow(incidence_mat))

  for (i in seq_len(nrow(incidence_mat))) {

    vars_i <- which(incidence_mat[i, ] == 1)

    match_idx <- which(
      sapply(oinfo_df$vars, function(v) {
        identical(
          sort(as.integer(strsplit(v, ",")[[1]])),
          sort(vars_i)
        )
      })
    )

    if (length(match_idx) > 0) {
      hyperedge_oinfo[i] <- oinfo_df$Oinfo[match_idx[1]]
    }
  }

  selected_indices <- seq_len(nrow(incidence_mat))

  if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
    oinfo_abs <- abs(hyperedge_oinfo)
    # 过滤掉NA值
    valid_idx <- which(!is.na(oinfo_abs))
    if (length(valid_idx) > 0) {
      top_indices <- valid_idx[order(oinfo_abs[valid_idx], decreasing = TRUE)[1:min(top_n, length(valid_idx))]]
      selected_indices <- intersect(selected_indices, top_indices)
      message("Selected top ", length(top_indices), " hyperedges by absolute O-info")
    }
  }

  if (!is.null(oinfo_threshold) && is.numeric(oinfo_threshold)) {
    oinfo_abs <- abs(hyperedge_oinfo)
    threshold_indices <- which(oinfo_abs >= oinfo_threshold & !is.na(oinfo_abs))
    selected_indices <- intersect(selected_indices, threshold_indices)
    message("Selected ", length(threshold_indices), " hyperedges with |O-info| >= ", oinfo_threshold)
  }

  if (length(selected_indices) == 0) {
    warning("No hyperedges remaining after filtering")
    return(invisible(NULL))
  }

  filtered_incidence <- incidence_mat[selected_indices, , drop = FALSE]
  filtered_oinfo <- hyperedge_oinfo[selected_indices]

  if (!is.null(groups)) {
    if (length(groups) != x$n_variables) {
      stop(
        "groups must have length equal to number of variables (",
        x$n_variables, ")."
      )
    }
    groups <- as.character(groups)
    group_levels <- unique(groups)
    n_groups <- length(group_levels)
  } else {
    group_levels <- NULL
    n_groups <- 0
  }

  hg_obj <- HyperG::hypergraph_from_incidence_matrix(
    filtered_incidence
  )

  n_vertices_actual <- ncol(filtered_incidence)

  if (!is.null(colnames(x$normalized_data))) {
    vertex_labels_all <- colnames(x$normalized_data)
  } else {
    vertex_labels_all <- as.character(seq_len(x$n_variables))
  }

  if (is.null(labels)) {
    vertex_labels <- vertex_labels_all
  } else {
    if (length(labels) != x$n_variables) {
      stop("Labels must have length equal to number of variables (",
           x$n_variables, "). Got ", length(labels))
    }
    vertex_labels <- labels
  }

  default_vertex_color <- "white"

  if (is.null(groups)) {
    if (is.null(vertex.color)) {
      vertex_color_vec <- rep(
        default_vertex_color,
        n_vertices_actual
      )
    } else if (length(vertex.color) == 1) {
      vertex_color_vec <- rep(
        vertex.color,
        n_vertices_actual
      )
    } else if (length(vertex.color) == x$n_variables) {
      vertex_color_vec <- vertex.color
    } else if (length(vertex.color) == n_vertices_actual) {
      vertex_color_vec <- vertex.color
    } else {
      warning(
        "vertex.color length not recognized; using white."
      )
      vertex_color_vec <- rep(
        default_vertex_color,
        n_vertices_actual
      )
    }

  } else {
    groups_used <- groups

    if (is.null(vertex.color)) {
      group_palette <- stats::setNames(
        grDevices::hcl.colors(
          n_groups,
          palette = "Set 3"
        ),
        group_levels
      )
      vertex_color_vec <- unname(
        group_palette[groups_used]
      )
    } else if (length(vertex.color) == n_groups) {
      group_palette <- stats::setNames(
        vertex.color,
        group_levels
      )
      vertex_color_vec <- unname(
        group_palette[groups_used]
      )
    } else if (length(vertex.color) == 1) {
      vertex_color_vec <- rep(
        vertex.color,
        n_vertices_actual
      )
    } else if (length(vertex.color) == x$n_variables) {
      vertex_color_vec <- vertex.color
    } else if (length(vertex.color) == n_vertices_actual) {
      vertex_color_vec <- vertex.color
    } else {
      warning(
        paste0(
          "vertex.color length must be 1, ",
          n_groups,
          " (groups), ",
          x$n_variables,
          " (all variables), or ",
          n_vertices_actual,
          " (displayed vertices). Using automatic group colors."
        )
      )

      group_palette <- stats::setNames(
        grDevices::hcl.colors(
          n_groups,
          palette = "Set 3"
        ),
        group_levels
      )

      vertex_color_vec <- unname(
        group_palette[groups_used]
      )
    }
  }

  if (length(vertex.size) == 1) {
    vertex_size_vec <- rep(vertex.size, n_vertices_actual)
  } else if (length(vertex.size) == x$n_variables) {
    vertex_size_vec <- vertex.size
  } else if (length(vertex.size) == n_vertices_actual) {
    vertex_size_vec <- vertex.size
  } else {
    warning("vertex.size length doesn't match number of vertices; using default")
    vertex_size_vec <- rep(15, n_vertices_actual)
  }

  mark_groups_arg <- lapply(
    seq_len(nrow(filtered_incidence)),
    function(i) which(filtered_incidence[i, ] == 1)
  )

  generate_edge_colors <- function(oinfo_vals, alpha = 0.3) {
    max_abs <- max(abs(oinfo_vals), na.rm = TRUE)
    norm_vals <- oinfo_vals / max_abs
    colors <- grDevices::colorRampPalette(c("red", "white", "blue"))(100)
    color_idx <- round(50 * (norm_vals + 1))  # [-1,1] -> [0,100]
    color_idx <- pmin(pmax(color_idx, 1), 100)
    grDevices::adjustcolor(colors[color_idx], alpha.f = alpha)
  }

  edge_colors <- generate_edge_colors(filtered_oinfo, alpha = edge.alpha)
  edge_border_colors <- generate_edge_colors(filtered_oinfo, alpha = 0.8)

  if (is.null(main)) {
    main <- switch(type,
                   full = sprintf("Full Hypergraph (%d hyperedges, %d vertices)",
                                  nrow(filtered_incidence), n_vertices_actual),
                   redundant = sprintf("Redundant Hypergraph (O-info > 0, %d hyperedges, %d vertices)",
                                       nrow(filtered_incidence), n_vertices_actual),
                   synergistic = sprintf("Synergistic Hypergraph (O-info < 0, %d hyperedges, %d vertices)",
                                         nrow(filtered_incidence), n_vertices_actual)
    )

    if (!is.null(top_n)) {
      main <- paste0(main, " [Top ", top_n, "]")
    }
    if (!is.null(oinfo_threshold)) {
      main <- paste0(main, " [|O-info| ≥ ", oinfo_threshold, "]")
    }
  }

  A <- sqrt(t(filtered_incidence) %*% filtered_incidence)
  diag(A) <- 0
  qgraph_ig <- qgraph::qgraph(A,layout = "spring", DoNotPlot = TRUE, repulsion = 0.8)
  if (is.null(layout)) {
    layout <- qgraph_ig$layout
  }


  plot_args <- list(
    x = hg_obj,
    mark.groups = mark_groups_arg,
    mark.col = edge_colors,
    mark.border = edge_border_colors,
    layout = layout,
    main = main,
    vertex.label = vertex_labels,
    vertex.color = vertex_color_vec,
    vertex.size = vertex_size_vec,
    vertex.label.cex = vertex.label.cex,
    vertex.label.color = vertex.label.color,
    vertex.shape = vertex.shape,
    edge.color = edge.color,
    edge.width = edge.width,
    edge.lty = edge.lty,
    ...
  )

  plot_args <- plot_args[!sapply(plot_args, is.null)]

  tryCatch({
    do.call(HyperG::plot.hypergraph, plot_args)
  }, error = function(e) {
    # if plotting fails, print a warning and try again with fewer parameters
    warning("Failed to plot with full parameters: ", e$message)
    warning("Attempting with simplified parameters...")

    simple_args <- list(
      x = hg_obj,
      mark.groups = mark_groups_arg,
      main = main,
      vertex.label = vertex_labels,
      vertex.color = vertex_color_vec,
      vertex.size = vertex_size_vec
    )
    do.call(HyperG::plot.hypergraph, simple_args)
  })

  invisible(list(
    hypergraph = hg_obj,
    layout = layout,
    incidence_matrix = filtered_incidence,
    oinfo = filtered_oinfo
  ))
}
