#' Generate All Candidate Variable Multiplets
#'
#' Generates all possible variable subsets of a given order k (k ≥ 3).
#'
#' This function implements the brute-force candidate generation
#' strategy described by Marinazzo et al. (2024), where every
#' possible combination of variables of size \code{k} is evaluated.
#' Although exhaustive, the number of candidate multiplets grows
#' combinatorially with the number of variables \eqn{\binom{p}{k}}
#'
#' where \eqn{p} is the number of variables.
#'
#' For large psychometric datasets this approach may become
#' computationally expensive, motivating more selective
#' candidate generation procedures such as those proposed by
#' Possenti et al. (2026).
#'
#' @param n_vars Total number of variables.
#' @param k Order of candidate multiplets (\eqn{k \ge 3}).
#'
#' @return A list of integer vectors, each representing a variable
#' subset of size \code{k}.
#'
#' @references
#' Marinazzo, D., Van Roozendaal, J., Rosas, F. E., Stella, M.,
#' Comolatti, R., Colenbier, N., Stramaglia, S., & Rosseel, Y. (2024).
#' An information-theoretic approach to build hypergraphs in psychometrics.
#' *Behavior Research Methods*, *56*(7), 8057-8079.
#' https://doi.org/10.3758/s13428-024-02471-8
#'
#' @examples
#' generate_all_multiplets(
#'   n_vars = 5,
#'   k = 3
#' )
#'
#' @export
generate_all_multiplets <- function(n_vars, k = 3) {
  if (k < 3) stop("k must be ≥ 3 for hypergraph construction")
  if (k > n_vars) stop("k cannot exceed total number of variables")

  # Generate all combinations
  combs <- utils::combn(1:n_vars, k, simplify = FALSE)
  return(combs)
}

#' Compute Community Affinity Matrix
#'
#' @noRd
compute_affinity_matrix <- function(g, membership){

  comms <- sort(unique(membership))
  n_comm <- length(comms)

  W <- matrix(0, n_comm, n_comm)

  Edf <- igraph::as_data_frame(g, what = "edges")

  for(i in seq_len(nrow(Edf))){

    v1 <- Edf$from[i]
    v2 <- Edf$to[i]

    c1 <- membership[as.character(v1)]
    c2 <- membership[as.character(v2)]

    w <- abs(Edf$weight[i])

    W[c1,c2] <- W[c1,c2] + w
    W[c2,c1] <- W[c2,c1] + w
  }

  counts <- matrix(1, n_comm, n_comm)

  for(i in seq_along(comms)){
    for(j in seq_along(comms)){

      ni <- sum(membership==comms[i])
      nj <- sum(membership==comms[j])

      counts[i,j] <- ni*nj
    }
  }

  W <- W / counts

  return(W)
}

#' Score a Clique Based on Community Affinity
#'
#' @noRd
score_clique <- function(
    clique,
    membership,
    W
){

  k <- length(clique)

  score <- 0

  pairs <- utils::combn(clique,2)

  for(i in seq_len(ncol(pairs))){

    v1 <- pairs[1,i]
    v2 <- pairs[2,i]

    c1 <- membership[as.character(v1)]
    c2 <- membership[as.character(v2)]

    score <- score + W[c1,c2]
  }

  score <- score / factorial(k)

  score
}

#' Greedy Expansion of a Clique
#'
#' @noRd
greedy_expand <- function(
    seed,
    g,
    membership,
    W,
    improvement = 0.05,
    max_size = 5
){

  current <- seed

  current_score <-
    score_clique(
      current,
      membership,
      W
    )

  repeat{

    if(length(current) >= max_size)
      break

    neighbors <- unique(unlist(
      lapply(current,function(v){

        as.integer(
          igraph::neighbors(
            g,
            v
          )
        )
      })
    ))

    neighbors <- setdiff(
      neighbors,
      current
    )

    if(length(neighbors)==0)
      break

    best_gain <- -Inf
    best_node <- NULL

    for(v in neighbors){

      new_set <- sort(c(current,v))

      new_score <-
        score_clique(
          new_set,
          membership,
          W
        )

      gain <- new_score-current_score

      if(gain > best_gain){

        best_gain <- gain
        best_node <- v
      }
    }

    if(best_gain < improvement)
      break

    current <- sort(c(current,best_node))

    current_score <- score_clique(
      current,
      membership,
      W
    )
  }

  current
}

#' Generate Candidate Hyperedges Using the Possenti et al. (2026) Pipeline
#'
#' Implements the candidate-selection strategy proposed by
#' Possenti et al. (2026) for identifying promising higher-order
#' interactions prior to O-information estimation.
#'
#' The method combines sparse pairwise network estimation, community
#' detection, clique scoring, and greedy expansion to reduce
#' the combinatorial burden associated with exhaustive
#' hyperedge search.
#'
#' The procedure consists of:
#'
#' \enumerate{
#' \item Estimating an EBICglasso network;
#' \item Detecting communities using the spin-glass algorithm;
#' \item Constructing a community affinity matrix;
#' \item Extracting maximal cliques;
#' \item Ranking cliques according to community affinity;
#' \item Selecting top-scoring seed cliques;
#' \item Greedily expanding seeds;
#' \item Generating all candidate subsets within the requested
#' order range.
#' }
#'
#' Candidate selection is motivated by the assumption that
#' meaningful higher-order dependencies are more likely to
#' occur among variables located in densely connected and
#' topologically coherent regions of the psychometric network.
#'
#' @param data Numeric data matrix or data frame.
#' @param k_range Vector of hyperedge orders to generate.
#' Typical values are \code{3:5}.
#' @param ebic.gamma EBIC tuning parameter for
#' \code{\link[qgraph]{EBICglasso}}.
#' @param top_seed_prop Proportion of top-scoring maximal
#' cliques retained as expansion seeds.
#' @param expansion_gain Minimum score improvement required
#' for adding a node during greedy expansion.
#'
#' @return A list containing:
#' \describe{
#'   \item{candidates}{Generated candidate hyperedges.}
#'   \item{adjacency_matrix}{EBICglasso adjacency matrix.}
#'   \item{graph}{Estimated igraph network.}
#'   \item{communities}{Community assignment vector.}
#'   \item{affinity_matrix}{Community affinity matrix.}
#'   \item{n_candidates}{Number of generated candidates.}
#' }
#'
#' @details
#' This candidate-generation procedure follows the conceptual
#' framework introduced by Possenti et al. (2026) to address
#' the exponential growth of candidate hyperedges in
#' psychometric hypergraph analysis.
#'
#' Rather than evaluating all
#' \eqn{\binom{p}{k}}
#' possible variable subsets, the method restricts the search
#' to topologically plausible regions of an EBICglasso network,
#' substantially reducing computational complexity while
#' preserving interpretable higher-order structures.
#'
#' @references
#' Possenti, F., Girelli, L., Tieri, P., & Petti, M. (2026). Multiplex
#' Hypergraph Modeling of Higher Order Structures in Psychometric Networks
#' (Version 1). arXiv. https://doi.org/10.48550/ARXIV.2604.22744
#'
#' @export
generate_candidate_multiplets_possenti <- function(
    data,
    k_range = 3:5,
    ebic.gamma = 0.5,
    top_seed_prop = 0.2,
    expansion_gain = 0.05
){
  S <- stats::cor(data)

  adj <-
    qgraph::EBICglasso(
      S,
      nrow(data),
      gamma = ebic.gamma
    )

  g <-
    igraph::graph_from_adjacency_matrix(
      adj,
      mode="undirected",
      weighted=TRUE,
      diag=FALSE
    )

  sg <-
    igraph::cluster_spinglass(
      g,
      weights = abs(igraph::E(g)$weight)
    )

  membership <- igraph::membership(sg)

  W <- compute_affinity_matrix(
    g,
    membership
  )

  cliques <-
    igraph::max_cliques(
      g,
      min=min(k_range),
      max=max(k_range)
    )

  cliques <- lapply(
    cliques,
    as.integer
  )

  scores <- sapply(
    cliques,
    score_clique,
    membership = membership,
    W = W
  )

  n_keep <-
    ceiling(
      length(cliques) *
        top_seed_prop
    )

  ord <- order(scores,decreasing=TRUE)

  seeds <- cliques[
    ord[1:n_keep]
  ]

  expanded <- lapply(
    seeds,
    greedy_expand,
    g = g,
    membership = membership,
    W = W,
    improvement = expansion_gain,
    max_size=max(k_range)
  )

  candidates <- list()

  for(mp in expanded){

    for(k in k_range){

      if(length(mp) >= k){

        candidates <-
          c(
            candidates,
            utils::combn(
              mp,
              k,
              simplify=FALSE
            )
          )
      }
    }
  }

  candidates <-
    unique(
      lapply(
        candidates,
        sort
      )
    )

  list(
    candidates = candidates,
    adjacency_matrix = adj,
    graph = g,
    communities = membership,
    affinity_matrix = W,
    n_candidates = length(candidates)
  )

}

#' Generate Candidate Hyperedges from a Pairwise Network
#'
#' Generates candidate variable multiplets using a sparse
#' Gaussian graphical model estimated with EBICglasso,
#' followed by maximal clique extraction and local clique expansion.
#'
#' This procedure is inspired by the candidate-selection strategy
#' proposed by Possenti et al. (2026) for reducing the combinatorial
#' search space of higher-order interactions.
#'
#' The algorithm consists of:
#'
#' \enumerate{
#' \item Estimating a pairwise network using EBICglasso;
#' \item Extracting maximal cliques from the network;
#' \item Generating candidate multiplets from cliques;
#' \item Expanding small cliques using neighboring nodes;
#' \item Removing duplicated candidates.
#' }
#'
#' Compared with exhaustive enumeration, this approach focuses
#' the search on densely connected regions of the psychometric
#' network and substantially reduces computational cost.
#'
#' @param data Numeric data matrix or data frame.
#' Rows correspond to observations and columns to variables.
#' @param k Target hyperedge order.
#' @param ebic.gamma EBIC tuning parameter used in
#' \code{\link[qgraph]{EBICglasso}}.
#' @param min_clique_size Minimum maximal clique size used as
#' candidate seeds.
#'
#' @return A list containing:
#' \describe{
#'   \item{candidates}{Candidate variable subsets.}
#'   \item{adjacency_matrix}{Estimated EBICglasso adjacency matrix.}
#'   \item{n_candidates}{Number of generated candidates.}
#' }
#'
#' @details
#' This implementation provides a lightweight approximation
#' of the candidate-generation framework described by
#' Possenti et al. (2026). For a closer implementation of the
#' original procedure, see
#' \code{\link{generate_candidate_multiplets_possenti}}.
#'
#' @references
#' Possenti, F., Girelli, L., Tieri, P., & Petti, M. (2026). Multiplex
#' Hypergraph Modeling of Higher Order Structures in Psychometric Networks
#' (Version 1). arXiv. https://doi.org/10.48550/ARXIV.2604.22744
#'
#' @export
generate_candidate_multiplets <- function(
    data,
    k = 3,
    ebic.gamma = 0.5,
    min_clique_size = 2
) {
  # Step 1: Estimate network with EBICglasso
  S <- stats::cor(data)
  n <- nrow(data)
  ebic_result <- qgraph::EBICglasso(
    S = S, n = n,
    gamma = ebic.gamma
  )
  adj_mat <- ebic_result
  g <- igraph::graph_from_adjacency_matrix(
    adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE
  )

  # Step 2: Extract maximal cliques
  cliques <- igraph::max_cliques(g, min = min_clique_size)
  clique_list <- lapply(cliques, as.integer)

  # Step 3: Expand cliques to generate k-order candidate multiplets
  candidates <- list()
  for (clq in clique_list) {
    if (length(clq) >= k) {
      subcombs <- utils::combn(clq, k, simplify = FALSE)
      candidates <- c(candidates, subcombs)
    } else if (length(clq) >= 2) {
      # Greedy expansion
      neighbors <- unique(unlist(lapply(clq, function(v) {
        igraph::neighbors(g, v)$name
      })))
      pool <- unique(c(clq, as.integer(neighbors)))
      if (length(pool) >= k) {
        addcombs <- utils::combn(pool, k, simplify = FALSE)
        candidates <- c(candidates, addcombs)
      }
    }
  }

  # Deduplicate
  candidates <- unique(candidates)
  if (length(candidates) == 0) {
    warning("No candidates generated; falling back to full combination.")
    candidates <- generate_all_multiplets(ncol(data), k = k)
  }

  return(list(
    candidates = candidates,
    adjacency_matrix = adj_mat,
    n_candidates = length(candidates)
  ))
}
