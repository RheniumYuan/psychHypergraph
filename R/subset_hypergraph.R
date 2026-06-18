#' Subset psychHypergraph based on various criteria
#'
#' @param x psychHypergraph object
#' @param type "full" (default), "redundant", or "synergistic"
#' @param oinfo_threshold Minimum absolute O-info value to retain
#' @param oinfo_quantile Quantile threshold for O-info (e.g., 0.9 to keep top 10%)
#' @param top_n Keep only top N hyperedges by absolute O-info
#' @param k Filter hyperedges by order (e.g., 3 or 4)
#' @param nodes Filter hyperedges that include specific nodes (variable indices)
#' @param node_mode "any" (default) to keep hyperedges that include any of the specified nodes, or "all" to keep only those that include all specified nodes
#' @param significant_only Logical. If TRUE (default), only keep hyperedges that passed all significance tests
#' @return A new psychHypergraph object containing only the subset of hyperedges that meet the specified criteria
#' @export
subset_hypergraph <- function(
    x,
    type = c("full", "redundant", "synergistic"),
    oinfo_threshold = NULL,
    oinfo_quantile = NULL,
    top_n = NULL,
    k = NULL,
    nodes = NULL,
    node_mode = c("any", "all"),
    significant_only = TRUE
){

  if(!inherits(x, "psychHypergraph")){
    stop("x must be a psychHypergraph object.")
  }

  type <- match.arg(type)
  node_mode <- match.arg(node_mode)

  df <- x$all_test_results

  if(type == "redundant"){
    df <- df[df$Oinfo > 0, ]
  }

  if(type == "synergistic"){
    df <- df[df$Oinfo < 0, ]
  }

  # siginificant_only
  if(significant_only){
    df <- df[df$significant, ]
  }

  # abs(Oinfo)
  df$abs_Oinfo <- abs(df$Oinfo)

  if(!is.null(oinfo_threshold)){
    df <- df[df$abs_Oinfo >= oinfo_threshold, ]
  }

  # quantile
  if(!is.null(oinfo_quantile)){
    cutoff <- stats::quantile(
      df$abs_Oinfo,
      probs = oinfo_quantile,
      na.rm = TRUE
    )

    df <- df[df$abs_Oinfo >= cutoff, ]
  }

  if(!is.null(top_n)){
    ord <- order(
      df$abs_Oinfo,
      decreasing = TRUE
    )

    df <- df[
      ord[seq_len(min(top_n,length(ord)))],
    ]
  }

  if(!is.null(k)){
    edge_order <- sapply(
      strsplit(df$vars,","),
      length
    )

    df <- df[edge_order %in% k, ]
  }

  if(!is.null(nodes)){
    keep <- sapply(
      strsplit(df$vars,","),
      function(v){

        v <- as.integer(v)

        if(node_mode == "any"){
          any(nodes %in% v)
        }else{
          all(nodes %in% v)
        }
      }
    )

    df <- df[keep, ]
  }

  hyperedges <- lapply(
    strsplit(df$vars,","),
    function(x) as.integer(x)
  )

  names(hyperedges) <- paste0("HE", df$hyperedge_id)

  redundant_hypergraph <- hyperedges[df$Oinfo > 0]

  synergistic_hypergraph <- hyperedges[df$Oinfo < 0]

  out <- x

  out$hypergraph <- hyperedges
  out$redundant_hypergraph <- redundant_hypergraph
  out$synergistic_hypergraph <- synergistic_hypergraph
  out$hyperedge_details <- df
  out$n_hyperedges <- length(hyperedges)
  out$n_redundant <- length(redundant_hypergraph)
  out$n_synergistic <- length(synergistic_hypergraph)

  out
}
