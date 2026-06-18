#' Convert psychHypergraph to HyperG hypergraph object
#'
#' Converts a psychHypergraph object to the hypergraph format from the
#' HyperG package, supporting full, redundant, and synergistic sub-hypergraphs.
#'
#' @param object A psychHypergraph object returned by build_hypergraph()
#' @param which Character: which hypergraph to convert. Options are:
#'   "full" (default), "redundant", "synergistic"
#'
#' @return A HyperG hypergraph object, or a named list of hypergraph objects.
#' @export
as_HyperG <- function(object,
                      which = c("full", "redundant", "synergistic")) {

  if (!inherits(object, "psychHypergraph")) {
    stop("Input must be a psychHypergraph object.")
  }

  which <- match.arg(which)

  if (!requireNamespace("HyperG", quietly = TRUE)) {
    stop("Package 'HyperG' is required. Please install it first.")
  }

  convert_one <- function(edge_list, name) {
    if (length(edge_list) == 0) {
      warning("The '", name, "' hypergraph is empty. Returning an empty hypergraph.")
      return(HyperG::hypergraph_from_edgelist(list()))
    }
    HyperG::hypergraph_from_edgelist(edge_list)
  }

  switch(
    which,
    full = {
      convert_one(object$hypergraph, "full")
    },
    redundant = {
      convert_one(object$redundant_hypergraph, "redundant")
    },
    synergistic = {
      convert_one(object$synergistic_hypergraph, "synergistic")
    }
  )
}
