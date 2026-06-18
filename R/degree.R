#' Calculate Hypergraph Node Degree Centrality
#'
#' Compute node-level degree centrality measures in psychometric
#' hypergraphs.
#'
#' The function calculates both unweighted degree and weighted degree
#' for nodes participating in significant hyperedges.
#'
#' Degree measures are computed separately for:
#'
#' \itemize{
#'   \item Redundant hypergraph (positive O-information)
#'   \item Synergistic hypergraph (negative O-information)
#'   \item full hypergraph (all significant hyperedges)
#' }
#'
#' Unweighted degree corresponds to the number of hyperedges containing
#' a node. Weighted degree corresponds to the sum of absolute
#' O-information values across all hyperedges containing the node.
#'
#' @param object A \code{psychHypergraph} object produced by
#'   \code{\link{build_hypergraph}}.
#'
#' @param scale Scaling method applied separately within each
#'   hypergraph type and degree measure.
#'
#'   \describe{
#'     \item{"raw"}{Raw degree values.}
#'     \item{"zscore"}{Standardized z-scores.}
#'     \item{"relative"}{Min-max scaled values between 0 and 1.}
#'   }
#'
#' @details
#'
#' For a node \eqn{i}, degree is defined as:
#'
#' \deqn{
#' degree(i) = \sum_h I(i \in h)
#' }
#'
#' where \eqn{I(\cdot)} is an indicator function and \eqn{h}
#' denotes a significant hyperedge.
#'
#' Weighted degree is defined as:
#'
#' \deqn{
#' weighted\ degree(i) =
#' \sum_h |O_h| I(i \in h)
#' }
#'
#' where \eqn{O_h} is the O-information value of hyperedge
#' \eqn{h}.
#'
#' Higher values indicate nodes participating in more
#' higher-order interactions.
#'
#' @return
#' A data frame containing:
#'
#' \describe{
#'
#' \item{node_id}{
#' Numeric node index.
#' }
#'
#' \item{node_name}{
#' Node name.
#' }
#'
#' \item{graph_type}{
#' Hypergraph type:
#' \code{"redundant"},
#' \code{"synergistic"},
#' or \code{"full"}.
#' }
#'
#' \item{degree_type}{
#' Degree measure:
#' \code{"degree"} or
#' \code{"weighted_degree"}.
#' }
#'
#' \item{value}{
#' Degree centrality value.
#' }
#'
#' }
#'
#' @seealso
#' \code{\link{plot_node_degree}},
#' \code{\link{build_hypergraph}}
#'
#' @examples
#' \dontrun{
#'
#' deg <- node_degree(hg)
#'
#' head(deg)
#'
#' deg_z <- node_degree(
#'   hg,
#'   scale = "zscore"
#' )
#'
#' deg_relative <- node_degree(
#'   hg,
#'   scale = "relative"
#' )
#'
#' }
#'
#' @export
node_degree <- function(object, scale = c("raw","zscore","relative")){

  if(!inherits(object, "psychHypergraph")){
    stop("object must be a psychHypergraph.")
  }

  results <- object$all_test_results

  if(nrow(results) == 0){
    stop("No hyperedge results found.")
  }

  node_names <- colnames(object$normalized_data)
  n_nodes <- length(node_names)

  # helper
  calc_degree <- function(df){
    deg <- numeric(n_nodes)
    wdeg <- numeric(n_nodes)
    if(nrow(df) == 0){
      return(
        data.frame(
          node_id = seq_len(n_nodes),
          node_name = node_names,
          degree = 0,
          weighted_degree = 0
        )
      )
    }

    for(i in seq_len(nrow(df))){
      vars <- as.integer(
        strsplit(df$vars[i], ",")[[1]]
      )
      deg[vars] <- deg[vars] + 1
      wdeg[vars] <- wdeg[vars] + abs(df$Oinfo[i])
    }

    data.frame(
      node_id = seq_len(n_nodes),
      node_name = node_names,
      degree = deg,
      weighted_degree = wdeg
    )
  }

  # redundant
  redundant_df <- results[results$significant & results$Oinfo > 0, ]
  redundant_res <- calc_degree(redundant_df)

  # synergistic
  synergistic_df <- results[results$significant & results$Oinfo < 0, ]
  synergistic_res <- calc_degree(synergistic_df)

  # full
  full_df <- results[results$significant, ]
  full_res <- calc_degree(full_df)

  # merge
  out <- rbind(
    data.frame(
      node_id = redundant_res$node_id,
      node_name = redundant_res$node_name,
      graph_type = "redundant",
      degree_type = "degree",
      value = redundant_res$degree
    ),
    data.frame(
      node_id = redundant_res$node_id,
      node_name = redundant_res$node_name,
      graph_type = "redundant",
      degree_type = "weighted_degree",
      value = redundant_res$weighted_degree
    ),
    data.frame(
      node_id = synergistic_res$node_id,
      node_name = synergistic_res$node_name,
      graph_type = "synergistic",
      degree_type = "degree",
      value = synergistic_res$degree
    ),
    data.frame(
      node_id = synergistic_res$node_id,
      node_name = synergistic_res$node_name,
      graph_type = "synergistic",
      degree_type = "weighted_degree",
      value = synergistic_res$weighted_degree
    ),
    data.frame(
      node_id = full_res$node_id,
      node_name = full_res$node_name,
      graph_type = "full",
      degree_type = "degree",
      value = full_res$degree
    ),
    data.frame(
      node_id = full_res$node_id,
      node_name = full_res$node_name,
      graph_type = "full",
      degree_type = "weighted_degree",
      value = full_res$weighted_degree
    )
  )


  # scaling
  scale <- match.arg(scale)
  if(scale != "raw"){
    out <- do.call(
      rbind,
      lapply(
        split(out, list(out$graph_type,
                        out$degree_type),
              drop = TRUE),
        function(df){
          x <- df$value
          if(scale == "zscore"){
            if(stats::sd(x) == 0){
              df$value <- 0
            }else{
              df$value <- as.numeric(scale(x))
            }
          }
          if(scale == "relative"){
            rng <- max(x) - min(x)
            if(rng == 0){
              df$value <- 0
            }else{
              df$value <- (x - min(x)) / rng
            }
          }
          df
        }
      )
    )

    rownames(out) <- NULL
  }

  out
}

#' Plot Hypergraph Node Degree Centrality
#'
#' Visualize node degree centrality measures in psychometric
#' hypergraphs.
#'
#' The function displays degree and weighted degree centrality
#' for redundant, synergistic, and full hypergraphs using
#' line plots, point plots, or a combination of both.
#'
#' @param object A \code{psychHypergraph} object produced by
#'   \code{\link{build_hypergraph}}.
#'
#' @param scale Scaling method passed to
#'   \code{\link{node_degree}}.
#'
#' @param degree_type Degree measure(s) to display.
#' One or both of:
#'
#' \itemize{
#'   \item \code{"degree"}
#'   \item \code{"weighted_degree"}
#' }
#'
#' @param graph_type Hypergraph type(s) to display.
#' Any combination of:
#'
#' \itemize{
#'   \item \code{"redundant"}
#'   \item \code{"synergistic"}
#'   \item \code{"full"}
#' }
#'
#' @param style Plot style.
#'
#' \describe{
#'   \item{"line"}{Line plot.}
#'   \item{"point"}{Point plot.}
#'   \item{"both"}{Line and point plot.}
#' }
#'
#' @param orientation Axis orientation.
#'
#' \describe{
#'   \item{"vertical"}{Nodes on the x-axis.}
#'   \item{"horizontal"}{Nodes on the y-axis.}
#' }
#'
#' @param rotate_labels Logical; rotate node labels to improve
#' readability.
#'
#' @param sort_nodes Logical; if \code{TRUE}, nodes are ordered
#' according to centrality values.
#'
#' @details
#'
#' The plot compares node centrality across different hypergraph
#' types and degree measures.
#'
#' Two centrality measures are available:
#'
#' \itemize{
#'   \item Degree: number of significant hyperedges involving a node.
#'   \item Weighted degree: sum of absolute O-information values
#'   across significant hyperedges involving a node.
#' }
#'
#' When multiple graph types are selected, they are shown together
#' using different colors.
#'
#' Separate panels are created for each degree measure.
#'
#' The resulting plot can be used to identify variables that
#' participate most strongly in higher-order redundant or
#' synergistic interactions.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' The returned object can be further customized using
#' standard \pkg{ggplot2} syntax.
#'
#' @seealso
#' \code{\link{node_degree}},
#' \code{\link{build_hypergraph}}
#'
#' @examples
#' \dontrun{
#'
#' # Default plot
#' plot_node_degree(hg)
#'
#' # Weighted degree only
#' plot_node_degree(
#'   hg,
#'   degree_type = "weighted_degree"
#' )
#'
#' # full hypergraph only
#' plot_node_degree(
#'   hg,
#'   graph_type = "full"
#' )
#'
#' # Horizontal layout
#' plot_node_degree(
#'   hg,
#'   orientation = "horizontal"
#' )
#'
#' # Sorted nodes
#' plot_node_degree(
#'   hg,
#'   sort_nodes = TRUE
#' )
#'
#' }
#'
#' @export
plot_node_degree <- function(
    object,
    scale = c("raw", "zscore", "relative"),
    degree_type = c("degree", "weighted_degree"),
    graph_type = c("redundant", "synergistic", "full"),
    style = c("line", "point", "both"),
    orientation = c("vertical", "horizontal"),
    rotate_labels = TRUE,
    sort_nodes = FALSE
){

  if(!inherits(object, "psychHypergraph")){
    stop("object must be a psychHypergraph.")
  }

  scale <- match.arg(scale)
  style <- match.arg(style)
  orientation <- match.arg(orientation)

  # degree data
  deg_df <- node_degree(
    object = object,
    scale = scale
  )
  deg_df <- deg_df[
    deg_df$degree_type %in% degree_type &
      deg_df$graph_type %in% graph_type,
  ]
  if(nrow(deg_df) == 0){
    stop("No data remaining after filtering.")
  }

  # pretty labels
  deg_df$degree_label <- dplyr::recode(
    deg_df$degree_type,
    degree = "Degree",
    weighted_degree = "Weighted Degree"
  )

  deg_df$graph_type <- factor(
    deg_df$graph_type,
    levels = c(
      "redundant",
      "synergistic",
      "full"
    )
  )

  # node order
  if(sort_nodes){
    tmp <- deg_df[deg_df$degree_type == degree_type[1] & deg_df$graph_type == graph_type[1],]
    node_order <- tmp$node_name[order(tmp$value, decreasing = TRUE)]
  }else{
    node_order <- colnames(object$normalized_data)
  }
  deg_df$node_name <- factor(deg_df$node_name, levels = node_order)

  if(orientation == "horizontal"){
    deg_df <- deg_df[order(deg_df$node_name), ]
  }

  # labels rotation
  if(rotate_labels && orientation == "vertical"){
    x_text <- ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  }else{
    x_text <- ggplot2::element_text()
  }

  # axis labels
  #ylab_text <- switch(
  #  scale,
  #  raw = "Degree Value",
  #  zscore = "Z-score Degree Value",
  #  relative = "Relative Degree Value"
  #)
  ylab_text <- "Degree Value"

  if(orientation == "vertical"){
    p <- ggplot2::ggplot(
      deg_df,
      ggplot2::aes(
        x = node_name,
        y = value,
        color = graph_type,
        group = graph_type
      )
    ) +
      ggplot2::facet_wrap(
        ~ degree_label,
        ncol = 1,
        scales = "free_y"
      ) +
      ggplot2::labs(
        x = "Node",
        y = ylab_text,
        color = "Hyperedge Type"
      )
  }else{
    p <- ggplot2::ggplot(
      deg_df,
      ggplot2::aes(
        x = value,
        y = node_name,
        color = graph_type,
        group = graph_type
      )
    ) +
      ggplot2::facet_grid(
        rows = NULL,
        cols = ggplot2::vars(degree_label),
        scales = "free_x"
      ) +
      ggplot2::labs(
        x = ylab_text,
        y = "Node",
        color = "Graph Type"
      )
  }

  if(style == "line" || style == "both"){
    if(orientation == "vertical"){
      # 垂直布局使用 geom_line
      p <- p + ggplot2::geom_line(linewidth = 0.8)
    } else {
      p <- p + ggplot2::geom_path(
        data = deg_df[order(deg_df$node_name), ],
        linewidth = 0.8,
        mapping = ggplot2::aes(
          x = value,
          y = node_name,
          color = graph_type,
          group = graph_type
        )
      )
    }
  }

  if(style == "point" || style == "both"){
    p <- p + ggplot2::geom_point(size = 2.5)
  }

  # final theme
  p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background =
        ggplot2::element_rect(
          fill = "grey95"
        ),
      strip.text =
        ggplot2::element_text(
          face = "bold"
        ),
      panel.grid.minor =
        ggplot2::element_blank(),
      axis.text.x = x_text
    )

  # additional adjustment for horizontal orientation
  if(orientation == "horizontal"){
    p <- p +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(
          size = ifelse(rotate_labels, 8, 10)
        )
      )
  }

  return(p)

}
