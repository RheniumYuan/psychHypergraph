#' Hypergraph Spectral Clustering
#'
#' Detect node communities in psychometric hypergraphs using
#' hypergraph spectral clustering.
#'
#' The implementation follows the normalized hypergraph Laplacian
#' framework proposed by Zhou, Huang, and Schölkopf (2006),
#' where clustering is performed in a low-dimensional spectral
#' embedding obtained from the hypergraph Laplacian. This function is
#' similar to hypergraph community detection module \code{HySC} in python
#' library \code{hypergraphx} (HGX).
#'
#' Nodes are first embedded using the eigenvectors associated with
#' the smallest non-trivial eigenvalues of the Laplacian matrix,
#' and then partitioned into communities using K-means clustering.
#'
#' @param psych_hg A \code{psychHypergraph} object produced by
#'   \code{\link{build_hypergraph}}.
#'
#' @param K Integer specifying the number of clusters.
#'
#' @param weighted_l Logical; if \code{TRUE}, construct a weighted
#'   incidence matrix using O-information values. If \code{FALSE},
#'   use a binary incidence matrix.
#'
#' @param norm_laplacian Logical; if \code{TRUE}, use the normalized
#'   hypergraph Laplacian proposed by Zhou et al. (2006).
#'
#' @param seed Random seed used for K-means initialization.
#'
#' @param n_init Number of random K-means starts.
#'
#' @details
#'
#' Let \eqn{H} denote the hypergraph incidence matrix,
#' \eqn{D_v} the node degree matrix,
#' and \eqn{D_e} the hyperedge degree matrix.
#'
#' The normalized hypergraph Laplacian is defined as:
#'
#' \deqn{
#' L =
#' I -
#' D_v^{-1/2}
#' H
#' D_e^{-1}
#' H^T
#' D_v^{-1/2}
#' }
#'
#' Spectral clustering proceeds as follows:
#'
#' \enumerate{
#'   \item Construct the incidence matrix.
#'   \item Compute the hypergraph Laplacian.
#'   \item Extract the smallest eigenvectors.
#'   \item Normalize rows of the embedding matrix.
#'   \item Apply K-means clustering.
#' }
#'
#' Nodes that do not participate in any hyperedge are treated as
#' isolated nodes and assigned cluster membership 0.
#'
#' @return
#' An object of class \code{"psychHypergraphClustering"} containing:
#'
#' \describe{
#'
#' \item{membership}{Cluster assignment vector.}
#' \item{membership_matrix}{Binary node-by-cluster membership matrix.}
#' \item{eigenvalues}{Smallest Laplacian eigenvalues.}
#' \item{eigenvectors}{Corresponding eigenvectors.}
#' \item{laplacian}{Hypergraph Laplacian matrix.}
#' \item{embed_vectors}{Spectral embedding used for clustering.}
#' \item{clusters}{List of nodes belonging to each cluster.}
#' \item{summary}{Cluster size summary.}
#' \item{parameters}{Analysis settings and clustering information.}
#' }
#'
#' @references
#' Zhou, D., Huang, J., & Schölkopf, B. (2007). Learning with Hypergraphs:
#' Clustering, Classification, and Embedding. In B. Schölkopf, J. Platt, & T.
#' Hofmann (Eds), *Advances in Neural Information Processing Systems 19 (pp.
#' 1601–1608)*. The MIT Press. https://doi.org/10.7551/mitpress/7503.003.0205
#'
#' @seealso
#' \code{\link{suggest_k}},
#' \code{\link{plot.psychHypergraphClustering}}
#'
#' @examples
#' \dontrun{
#'
#' cl <- spectral_clustering(
#'   psych_hg = hg,
#'   K = 3
#' )
#'
#' cl
#'
#' plot(cl)
#'
#' }
#'
#' @export
spectral_clustering <- function(psych_hg,
                                K = 2,
                                weighted_l = FALSE,
                                norm_laplacian = TRUE,
                                seed = 42,
                                n_init = 10) {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is needed for sparse matrix operations")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("Package 'stats' is needed for kmeans")
  }

  # 1. hypergraph representation
  hyperedges <- psych_hg$hypergraph

  N <- psych_hg$n_variables

  if (weighted_l) {
    H <- build_weighted_incidence_matrix(hyperedges, psych_hg$hyperedge_details, N)
  } else {
    H <- build_binary_incidence_matrix(hyperedges, N)
  }

  # 2. calculate node degrees and hyperedge sizes
  if (weighted_l) {
    D_v <- Matrix::rowSums(H)
    D_e <- Matrix::colSums(H)
  } else {
    D_v <- Matrix::rowSums(H)
    D_e <- Matrix::colSums(H)
  }

  isolates <- which(D_v == 0)
  non_isolates <- which(D_v > 0)

  if (length(isolates) > 0) {
    warning(sprintf("Found %d isolated nodes (degree = 0). They will be assigned to cluster 0.",
                    length(isolates)))
  }

  # 3. Laplacian construction
  # L = I - D_v^{-1/2} * H * D_e^{-1} * H^T * D_v^{-1/2}

  inv_sqrt_Dv <- 1 / sqrt(D_v)
  inv_sqrt_Dv[is.infinite(inv_sqrt_Dv) | D_v == 0] <- 0

  inv_De <- 1 / D_e
  inv_De[is.infinite(inv_De) | D_e == 0] <- 0

  Dv_sqrt_inv <- Matrix::Diagonal(x = inv_sqrt_Dv)
  De_inv <- Matrix::Diagonal(x = inv_De)

  # H %*% De_inv %*% t(H)
  H_sparse <- Matrix::Matrix(H, sparse = TRUE)

  H_De_inv_Ht <- H_sparse %*% De_inv %*% Matrix::t(H_sparse)

  L <- Matrix::Diagonal(N) - Dv_sqrt_inv %*% H_De_inv_Ht %*% Dv_sqrt_inv

  L <- as.matrix(L)

  # 4. handle isolated nodes by subsetting the Laplacian to non-isolated nodes
  if (length(non_isolates) > 0) {
    L_sub <- L[non_isolates, non_isolates, drop = FALSE]
  } else {
    stop("All nodes are isolated! Cannot perform clustering.")
  }

  # 5. eigen decomposition
  eig <- eigen(L_sub, symmetric = TRUE)

  order_idx <- order(eig$values)
  eigenvalues <- eig$values[order_idx]
  eigenvectors <- eig$vectors[, order_idx, drop = FALSE]

  # the first eigenvector corresponds to the trivial solution (all ones) and
  #is not informative for clustering
  if (norm_laplacian) {
    if (K > ncol(eigenvectors)) {
      warning(sprintf("K=%d is larger than number of available eigenvectors (%d)",
                      K, ncol(eigenvectors)))
      K <- ncol(eigenvectors)
    }
    if (K >= 1) {
      embed_vecs <- eigenvectors[, 2:min(K+1, ncol(eigenvectors)), drop = FALSE]
    } else {
      embed_vecs <- matrix(0, nrow = nrow(eigenvectors), ncol = 1)
    }
  } else {
    embed_vecs <- eigenvectors[, 1:min(K, ncol(eigenvectors)), drop = FALSE]
  }

  # 6. normalize rows of the embedding matrix
  row_norms <- sqrt(rowSums(embed_vecs^2))
  row_norms[row_norms == 0] <- 1
  embed_vecs_norm <- embed_vecs / row_norms

  # 7. K-means clustering
  set.seed(seed)
  km_result <- stats::kmeans(embed_vecs_norm,
                      centers = K,
                      nstart = n_init,
                      iter.max = 100)

  # 8. construct membership vector and matrix
  membership <- rep(0, N)  # 0 for isolated nodes
  if (length(non_isolates) > 0) {
    membership[non_isolates] <- km_result$cluster
  }

  membership_matrix <- matrix(0, nrow = N, ncol = K)
  for (i in 1:N) {
    if (membership[i] > 0) {
      membership_matrix[i, membership[i]] <- 1
    }
  }

  clusters <- list()
  for (k in 1:K) {
    clusters[[k]] <- which(membership == k)
  }
  if (length(isolates) > 0) {
    clusters[[K+1]] <- isolates
    names(clusters)[K+1] <- "isolated"
  }

  # 9. summary statistics
  summary_stats <- data.frame(
    cluster = 1:K,
    size = sapply(clusters[1:K], length),
    proportion = sapply(clusters[1:K], length) / length(non_isolates)
  )

  # 10. final
  result <- list(
    membership = membership,
    membership_matrix = membership_matrix,
    eigenvalues = eigenvalues[1:min(K+1, length(eigenvalues))],
    eigenvectors = eigenvectors,
    laplacian = L,
    embed_vectors = embed_vecs_norm,
    clusters = clusters,
    summary = summary_stats,
    parameters = list(
      K = K,
      weighted = weighted_l,
      norm_laplacian = norm_laplacian,
      n_nodes = N,
      n_hyperedges = length(hyperedges),
      n_non_isolates = length(non_isolates),
      n_isolates = length(isolates)
    ),
    kmeans_object = km_result
  )

  class(result) <- "psychHypergraphClustering"

  return(result)
}

#' Build binary incidence matrix from hypergraph
#'
#' @param hyperedges List of hyperedges, each a vector of node indices
#' @param N Number of nodes
#' @return Binary incidence matrix (N x E)
#' @noRd
build_binary_incidence_matrix <- function(hyperedges, N) {
  E <- length(hyperedges)
  H <- Matrix::Matrix(0, nrow = N, ncol = E, sparse = TRUE)

  for (e_idx in 1:E) {
    edge_nodes <- hyperedges[[e_idx]]
    if (length(edge_nodes) > 0) {
      H[edge_nodes, e_idx] <- 1
    }
  }

  return(H)
}

#' Build weighted incidence matrix from hypergraph
#'
#' @param hyperedges List of hyperedges
#' @param hyperedge_details Data frame with hyperedge details including Oinfo
#' @param N Number of nodes
#' @return Weighted incidence matrix (N x E)
#' @noRd
build_weighted_incidence_matrix <- function(hyperedges, hyperedge_details, N) {
  E <- length(hyperedges)
  H <- Matrix::Matrix(0, nrow = N, ncol = E, sparse = TRUE)

  edge_weights <- rep(1, E)

  if (!is.null(hyperedge_details) && nrow(hyperedge_details) == E) {
    oinfo <- hyperedge_details$Oinfo
    if (!all(is.na(oinfo))) {
      edge_weights <- oinfo - min(oinfo)
      if (max(edge_weights) > 0) {
        edge_weights <- edge_weights / max(edge_weights)
      } else {
        edge_weights <- rep(1, E)
      }
    }
  }

  for (e_idx in 1:E) {
    edge_nodes <- hyperedges[[e_idx]]
    if (length(edge_nodes) > 0) {
      H[edge_nodes, e_idx] <- edge_weights[e_idx]
    }
  }

  return(H)
}

#' Print method for psychHypergraphClustering
#'
#' @param x psychHypergraphClustering object
#' @param ... Additional arguments
#' @export
print.psychHypergraphClustering <- function(x, ...) {
  cat("Hypergraph Spectral Clustering Results\n")
  cat("=====================================\n")
  cat(sprintf("Number of nodes: %d\n", x$parameters$n_nodes))
  cat(sprintf("Number of hyperedges: %d\n", x$parameters$n_hyperedges))
  cat(sprintf("Number of clusters (K): %d\n", x$parameters$K))
  cat(sprintf("Clustering method: %s Laplacian\n",
              ifelse(x$parameters$weighted, "Weighted", "Binary")))
  cat(sprintf("Non-isolated nodes: %d\n", x$parameters$n_non_isolates))
  cat(sprintf("Isolated nodes: %d\n", x$parameters$n_isolates))
  cat("\nCluster sizes:\n")
  print(x$summary)
  cat("\nFirst few eigenvalues:\n")
  print(round(x$eigenvalues[1:min(5, length(x$eigenvalues))], 4))
}

#' Plot Hypergraph Clustering Results
#'
#' Visualize the results of hypergraph spectral clustering.
#'
#' The function can display:
#'
#' \itemize{
#'   \item Laplacian eigenvalues
#'   \item Spectral embedding colored by cluster membership
#'   \item Both visualizations
#' }
#'
#' @param x A \code{psychHypergraphClustering} object.
#'
#' @param type Type of visualization.
#'
#' \describe{
#'   \item{"eigenvalues"}{
#'   Plot Laplacian eigenvalues.
#'   }
#'
#'   \item{"clusters"}{
#'   Plot spectral embedding colored by cluster.
#'   }
#'
#'   \item{"both"}{
#'   Produce both plots.
#'   }
#' }
#'
#' @param ... Additional arguments.
#'
#' @details
#'
#' The eigenvalue plot can be used to inspect the
#' spectral structure of the hypergraph.
#'
#' The embedding plot shows node positions in the
#' low-dimensional spectral space obtained from the
#' first non-trivial eigenvectors of the Laplacian.
#'
#' Nodes assigned to the same cluster should appear
#' close together in the embedding space.
#'
#' @return
#'
#' If \code{type = "eigenvalues"},
#' returns a ggplot object.
#'
#' If \code{type = "clusters"},
#' returns a ggplot object.
#'
#' If \code{type = "both"},
#' returns a list containing both plots.
#'
#' @seealso
#' \code{\link{spectral_clustering}}
#'
#' @examples
#' \dontrun{
#'
#' cl <- spectral_clustering(hg, K = 3)
#'
#' plot(cl)
#'
#' plot(cl, type = "eigenvalues")
#'
#' plot(cl, type = "clusters")
#'
#' }
#'
#' @export
plot.psychHypergraphClustering <- function(x, type = "both", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is needed for plotting")
    return(invisible(NULL))
  }

  if (type == "eigenvalues" || type == "both") {

    eig_df <- data.frame(
      index = 1:length(x$eigenvalues),
      value = x$eigenvalues
    )
    p1 <- ggplot2::ggplot(eig_df, ggplot2::aes(x = index, y = value)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_line() +
      ggplot2::labs(title = "Laplacian Eigenvalues",
                    x = "Index", y = "Eigenvalue") +
      ggplot2::theme_minimal()

    if (type == "eigenvalues") {
      print(p1)
      return(invisible(p1))
    }
  }

  if (type == "clusters" || type == "both") {

    if (ncol(x$embed_vectors) >= 2) {
      embed_df <- data.frame(
        dim1 = x$embed_vectors[, 1],
        dim2 = x$embed_vectors[, 2],
        cluster = factor(x$membership[x$membership > 0])
      )

      p2 <- ggplot2::ggplot(embed_df, ggplot2::aes(x = dim1, y = dim2, color = cluster)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(title = "Spectral Embedding (First Two Eigenvectors)",
                      x = "First eigenvector", y = "Second eigenvector") +
        ggplot2::theme_minimal()

      print(p2)

      if (type == "both") {
        return(invisible(list(eigen_plot = p1, embed_plot = p2)))
      } else {
        return(invisible(p2))
      }
    } else {
      warning("Insufficient dimensions for embedding plot")
    }
  }
}

#' Suggest the Number of Clusters Using the Eigengap Heuristic
#'
#' Estimate an appropriate number of clusters for hypergraph
#' spectral clustering using the eigengap heuristic.
#'
#' The method computes Laplacian eigenvalues and identifies
#' the largest gap between consecutive eigenvalues.
#'
#' A large eigengap indicates a natural partition of the
#' hypergraph into communities.
#'
#' @param psych_hg A \code{psychHypergraph} object.
#'
#' @param max_K Maximum number of clusters to evaluate.
#'
#' @param weighted_l Logical; whether to use a weighted
#' incidence matrix.
#'
#' @param ... Additional arguments.
#'
#' @details
#'
#' Let
#' \eqn{\lambda_1 \le \lambda_2 \le \cdots}
#' denote Laplacian eigenvalues.
#'
#' The eigengap for cluster number \eqn{k} is:
#'
#' \deqn{
#' gap_k = \lambda_{k+1} - \lambda_k
#' }
#'
#' The suggested number of clusters corresponds to the
#' largest observed eigengap.
#'
#' @return
#' A list containing:
#'
#' \describe{
#' \item{suggested_K}{Recommended number of clusters.}
#' \item{eigenvalues}{Evaluated eigenvalues.}
#' \item{eigengaps}{Consecutive eigenvalue differences.}
#' \item{eigenvalues_all}{Full eigenvalue spectrum.}
#' }
#'
#' @seealso
#' \code{\link{spectral_clustering}}
#'
#' @examples
#' \dontrun{
#'
#' suggest_k(hg)
#'
#' k_info <- suggest_k(
#'   hg,
#'   max_K = 8
#' )
#'
#' k_info$suggested_K
#'
#' }
#'
#' @export
suggest_k <- function(psych_hg, max_K = 10, weighted_l = FALSE, ...) {
  N <- psych_hg$n_variables
  hyperedges <- psych_hg$hypergraph

  H <- build_binary_incidence_matrix(hyperedges, N)
  D_v <- Matrix::rowSums(H)
  D_e <- Matrix::colSums(H)

  D_v[D_v == 0] <- 1

  inv_sqrt_Dv <- 1 / sqrt(D_v)
  inv_De <- 1 / D_e
  inv_De[is.infinite(inv_De) | D_e == 0] <- 0

  Dv_sqrt_inv <- Matrix::Diagonal(x = inv_sqrt_Dv)
  De_inv <- Matrix::Diagonal(x = inv_De)

  H_sparse <- Matrix::Matrix(H, sparse = TRUE)
  H_De_inv_Ht <- H_sparse %*% De_inv %*% Matrix::t(H_sparse)
  L <- Matrix::Diagonal(N) - Dv_sqrt_inv %*% H_De_inv_Ht %*% Dv_sqrt_inv

  L <- as.matrix(L)
  non_isolates <- which(D_v > 0)

  if (length(non_isolates) > 0) {
    L_sub <- L[non_isolates, non_isolates]
    eig <- eigen(L_sub, symmetric = TRUE)
    eigenvalues <- eig$values

    n_eig <- min(max_K + 1, length(eigenvalues))
    eigengaps <- diff(eigenvalues[1:n_eig])

    suggested_K <- which.max(eigengaps)

    result <- list(
      suggested_K = suggested_K,
      eigenvalues = eigenvalues[1:n_eig],
      eigengaps = eigengaps,
      eigenvalues_all = eigenvalues
    )

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      gap_df <- data.frame(
        K = 1:length(eigengaps),
        gap = eigengaps
      )
      p <- ggplot2::ggplot(gap_df, ggplot2::aes(x = K, y = gap)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept = suggested_K, linetype = "dashed", color = "red") +
        ggplot2::labs(title = "Eigengap Heuristic for K Selection",
                      x = "Number of clusters (K)",
                      y = "Eigenvalue gap") +
        ggplot2::theme_minimal()
      print(p)
    }

    return(result)
  } else {
    warning("All nodes are isolated")
    return(list(suggested_K = 1, eigenvalues = numeric(0), eigengaps = numeric(0)))
  }
}
