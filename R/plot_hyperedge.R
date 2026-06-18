#' Plot Hyperedge O-information Estimates
#'
#' Visualize O-information estimates for hyperedges identified during
#' hypergraph construction. Each hyperedge is represented by its
#' O-information value and, optionally, its bootstrap confidence interval.
#'
#' Positive O-information indicates redundant higher-order interactions,
#' whereas negative O-information indicates synergistic interactions.
#'
#' The resulting figure resembles a forest plot, allowing users to inspect
#' the magnitude, direction, and uncertainty of hyperedge-level
#' O-information estimates.
#'
#' @param hg A \code{psychHypergraph} object produced by
#'   \code{\link{build_hypergraph}}.
#'
#' @param significant_only Logical; if \code{TRUE}, only hyperedges passing
#'   the significance testing pipeline are displayed.
#'
#' @param sort Logical; if \code{TRUE}, hyperedges are sorted by
#'   O-information values from largest to smallest.
#'
#' @param show_ci Logical; if \code{TRUE}, bootstrap confidence intervals
#'   are displayed as horizontal error bars.
#'
#' @param color_type Logical; if \code{TRUE}, points are colored according
#'   to hyperedge type:
#'   \itemize{
#'     \item Redundant (positive O-information)
#'     \item Synergistic (negative O-information)
#'   }
#'
#' @param show_edge_label Logical; if \code{TRUE}, hyperedge labels are
#'   displayed on the y-axis.
#'
#' @param point_size Numeric value controlling point size.
#'
#' @param redundant_color Color for redundant hyperedges (positive O-information).
#'
#' @param synergistic_color Color for synergistic hyperedges (negative O-information).
#'
#' @details
#'
#' Hyperedges are plotted on the y-axis and O-information values on the
#' x-axis. A vertical dashed line at zero separates redundant
#' (\eqn{O > 0}) and synergistic (\eqn{O < 0}) interactions.
#'
#' When bootstrap confidence intervals are available, they are displayed
#' using the BCa intervals stored in the psychHypergraph object.
#'
#' Hyperedge labels correspond to the variable indices defining each
#' hyperedge. For example:
#'
#' \preformatted{
#' 1-3-5
#' }
#'
#' represents the hyperedge formed by variables 1, 3, and 5.
#'
#' For large hypergraphs, displaying all hyperedges may result in crowded
#' figures. In such cases, setting \code{significant_only = TRUE} or
#' subsetting the hypergraph beforehand is recommended.
#'
#' @return
#' A \code{ggplot2} object.
#'
#' The returned object can be further customized using standard
#' \pkg{ggplot2} syntax.
#'
#' @seealso
#' \code{\link{build_hypergraph}},
#' \code{\link{plot.psychHypergraph}}
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
#' # Plot all hyperedges
#' plot_hyperedge(hg)
#'
#' }
#'
#' @export
plot_hyperedge <- function(
    hg,
    significant_only = FALSE,
    sort = TRUE,
    show_ci = TRUE,
    color_type = TRUE,
    show_edge_label = TRUE,
    point_size = 1.8,
    redundant_color = "#0072B2",
    synergistic_color = "#D55E00"
){

  if(!inherits(hg, "psychHypergraph")){
    stop("hg must be a psychHypergraph object.")
  }

  df <- hg$all_test_results

  if(significant_only){
    df <- df[df$significant, ]
  }

  if(nrow(df) == 0){
    stop("No hyperedges available.")
  }

  # Hyperedge labels
  df$edge_label <- gsub(
    pattern = ",",
    replacement = "-",
    x = df$vars
  )

  # Type
  df$type <- ifelse(
    df$Oinfo > 0,
    "Redundant",
    "Synergistic"
  )

  # Sort
  if(sort){

    ord <- order(
      df$Oinfo,
      decreasing = TRUE
    )

    df <- df[ord, ]
  }

  df$edge_label <- factor(
    df$edge_label,
    levels = rev(df$edge_label)
  )

  # Base plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Oinfo,
      y = edge_label
    )
  )

  # Confidence intervals
  if(show_ci){

    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          xmin = boot_ci_lower,
          xmax = boot_ci_upper
        ),
        width = 0.15,
        linewidth = 0.4,
        alpha = 0.7
      )
  }

  # Points
  if(color_type){

    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(
          color = type
        ),
        size = point_size
      )

  } else {

    p <- p +
      ggplot2::geom_point(
        size = point_size
      )
  }

  # Zero line
  p <- p +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    )

  # Colors
  if(color_type){

    p <- p +
      ggplot2::scale_color_manual(
        values = c(
          Redundant = redundant_color,
          Synergistic = synergistic_color
        )
      )
  }

  # Axis labels
  p <- p +
    ggplot2::scale_y_discrete(
      labels = if(show_edge_label) levels(df$edge_label) else rep("", nrow(df))
    )

  # Theme
  p <- p +
    ggplot2::labs(
      x = "O-information",
      y = "Hyperedge",
      color = "Type"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right"
    )

  return(p)
}
