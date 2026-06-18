#' Build Information-Theoretic Psychometric Hypergraph
#'
#' Constructs psychometric hypergraphs based on higher-order statistical
#' dependencies quantified by O-information. The function implements a complete
#' analysis pipeline inspired by Marinazzo et al. (2024) and extended by
#' Possenti et al. (2026), including:
#'
#' \enumerate{
#'   \item Data preprocessing and correlation estimation
#'   \item Candidate hyperedge generation
#'   \item O-information computation
#'   \item Permutation-based significance testing
#'   \item Bootstrap-based robustness assessment
#'   \item Hierarchical validation of hyperedges
#'   \item Multiple-testing correction
#'   \item Construction of redundant and synergistic hypergraphs
#' }
#'
#' Hyperedges represent groups of variables exhibiting statistically significant
#' higher-order interactions. Positive O-information indicates dominant
#' redundancy, whereas negative O-information indicates dominant synergy.
#'
#' @param data A data frame or matrix containing psychometric observations.
#' Rows correspond to subjects and columns correspond to variables/items.
#'
#' @param cor_method Correlation estimation method.
#' \describe{
#'   \item{"pearson"}{Pearson correlation matrix.}
#'   \item{"polychoric"}{Polychoric correlation matrix for ordinal variables.}
#' }
#'
#' @param k Integer vector specifying hyperedge orders to evaluate.
#' For example, `k = 3:5` evaluates all candidate triplets, quadruplets,
#' and quintuplets.
#'
#' @param candidate_method Method used to generate candidate multiplets.
#' \describe{
#'   \item{"full"}{Exhaustive search over all possible combinations
#'   (Marinazzo et al., 2024). See \code{\link{generate_all_multiplets}}.}
#'
#'   \item{"clique"}{Generate candidates from maximal cliques identified
#'   in an EBICglasso network. See \code{\link{generate_candidate_multiplets}}.}
#'
#'   \item{"possenti"}{Topology-guided candidate generation strategy
#'   inspired by Possenti et al. (2026). See \code{\link{generate_candidate_multiplets_possenti}}.}
#' }
#'
#' @param user_multiplets Optional list of user-defined multiplets.
#' Each element should be an integer vector containing node indices.
#' User-defined multiplets are added to automatically generated candidates.
#'
#' @param ebic.gamma EBIC tuning parameter used when estimating
#' pairwise networks using EBICglasso.
#'
#' @param min_clique_size Minimum clique size retained during
#' clique-based candidate generation.
#'
#' @param n_perm Number of permutations used in significance testing.
#'
#' @param alpha_perm Significance threshold for permutation testing.
#'
#' @param Oinfo_threshold Minimum absolute O-information value required
#' before a hyperedge can be considered significant.
#'
#' @param n_boot Number of bootstrap samples used for robustness assessment.
#'
#' @param alpha_boot Significance threshold for bootstrap confidence intervals.
#'
#' @param correct_method Multiple-testing correction method.
#' One of `"fdr"` or `"holm"`.
#'
#' @param handle_na Missing-data handling strategy.
#' \describe{
#'   \item{"impute"}{Impute missing values using `impute_fun`.}
#'   \item{"omit"}{Remove observations containing missing values.}
#' }
#'
#' @param impute_fun Function used for imputation when
#' `handle_na = "impute"`.
#'
#' @param seed Random seed for reproducibility.
#'
#' @param parallel Logical. Should computations be performed in parallel?
#'
#' @param n_cores Number of CPU cores used when `parallel = TRUE`.
#'
#' @details
#'
#' The analysis proceeds in three major stages:
#'
#' \strong{Stage 1: Data preprocessing}
#'
#' Variables are optionally transformed and a correlation matrix is estimated.
#'
#' \strong{Stage 2: Candidate generation}
#'
#' Candidate hyperedges are generated using either exhaustive enumeration,
#' network-based clique expansion, or topology-guided search.
#'
#' \strong{Stage 3: Hyperedge construction}
#'
#' Hyperedge construction follows a rigorous three-step testing pipeline proposed
#' by Possenti et al. (2026).
#'
#' O-information is computed for each candidate multiplet and evaluated using
#' permutation tests, bootstrap stability assessment, and hierarchical validation.
#' Significant multiplets are retained as hyperedges. Detailed three-step testing
#' pipelines are implemented in \code{\link{sig_multiplet_testing_pipeline}}.
#'
#' Hyperedges with positive O-information are assigned to the redundant
#' hypergraph, whereas hyperedges with negative O-information are assigned
#' to the synergistic hypergraph.
#'
#' @return
#' An object of class `"psychHypergraph"` containing:
#'
#' \describe{
#'
#' \item{hypergraph}{
#' List of all significant hyperedges.
#' }
#'
#' \item{redundant_hypergraph}{
#' Hyperedges with positive O-information.
#' }
#'
#' \item{synergistic_hypergraph}{
#' Hyperedges with negative O-information.
#' }
#'
#' \item{incidence_matrix}{
#' Hypergraph incidence matrix.
#' Rows correspond to hyperedges and columns correspond to nodes.
#' }
#'
#' \item{redundant_incidence_matrix}{
#' Incidence matrix for redundant hyperedges.
#' }
#'
#' \item{synergistic_incidence_matrix}{
#' Incidence matrix for synergistic hyperedges.
#' }
#'
#' \item{hyperedge_details}{
#' Data frame containing significant hyperedges and associated statistics.
#' }
#'
#' \item{all_test_results}{
#' Complete testing results for all evaluated candidate multiplets.
#' }
#'
#' \item{cor_matrix}{
#' Correlation matrix used in analysis.
#' }
#'
#' \item{ebicglasso_network}{
#' Estimated EBICglasso adjacency matrix (if applicable).
#' }
#'
#' \item{candidate_multiplets}{
#' All evaluated candidate multiplets.
#' }
#'
#' \item{parameters}{
#' Analysis settings used to generate the hypergraph.
#' }
#' }
#'
#' @references
#'
#' Marinazzo, D., Van Roozendaal, J., Rosas, F. E., Stella, M.,
#' Comolatti, R., Colenbier, N., Stramaglia, S., & Rosseel, Y. (2024).
#' An information-theoretic approach to build hypergraphs in psychometrics.
#' *Behavior Research Methods*, *56*(7), 8057-8079.
#'
#' Possenti, F., Girelli, L., Tieri, P., & Petti, M. (2026). Multiplex
#' Hypergraph Modeling of Higher Order Structures in Psychometric Networks
#' (Version 1). arXiv. https://doi.org/10.48550/ARXIV.2604.22744
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
#' hg
#'
#' plot(hg)
#'
#' }
#'
#' @export
build_hypergraph <- function(
    data,
    cor_method = c("pearson", "polychoric"),
    k = 3,
    candidate_method = c("full", "clique", "possenti"),
    user_multiplets = NULL,
    ebic.gamma = 0.5,
    min_clique_size = 2,
    n_perm = 1000,
    alpha_perm = 0.05,
    Oinfo_threshold = 0.01,
    n_boot = 1000,
    alpha_boot = 0.05,
    correct_method = c("fdr", "holm"),
    handle_na = c("impute", "omit"),
    impute_fun = stats::median,
    seed = 123,
    parallel = TRUE,
    n_cores = 1
) {

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data frame or matrix.")
  }
  if (any(k < 3)) stop("All k values must be ≥ 3.")
  if (!is.null(user_multiplets) && !is.list(user_multiplets)) {
    stop("user_multiplets must be a list of integer vectors.")
  }

  cor_method <- match.arg(cor_method)
  candidate_method <- match.arg(candidate_method)
  correct_method <- match.arg(correct_method)
  handle_na <- match.arg(handle_na)

  if(parallel){
    future::plan(
      future::multisession,
      workers = n_cores
    )
  }else{
    future::plan(
      future::sequential
    )
  }

  set.seed(seed)

  params <- as.list(environment())

  # Step 1: data preprocessing
  message("Step 0/3: Preprocessing data...")
  prep <- data_preprocess_pipeline(
    data = data,
    cor_method = cor_method,
    handle_na = handle_na,
    impute_fun = impute_fun
  )
  norm_data <- prep$normalized_data
  cor_mat <- prep$cor_matrix
  n_vars <- ncol(norm_data)

  # Step 2: candidate multiplet generation
  message("Generating candidate multiplets using method: ", candidate_method)
  candidate_mp <- list()
  ebic_net <- NULL

  if (candidate_method == "clique") {
    for (order in k) {
      if (order > n_vars) next
      mp_res <- generate_candidate_multiplets(
        data = norm_data,
        k = order,
        ebic.gamma = ebic.gamma,
        min_clique_size = min_clique_size
      )
      candidate_mp <- c(candidate_mp, mp_res$candidates)
      if (is.null(ebic_net)) ebic_net <- mp_res$adjacency_matrix
    }
  } else if (candidate_method == "full") {
    for (order in k) {
      if (order > n_vars) next
      candidate_mp <- c(candidate_mp, generate_all_multiplets(n_vars, k = order))
    }
  } else if (candidate_method == "possenti") {
    mp_res <- generate_candidate_multiplets_possenti(
      data = norm_data,
      k_range = k,
      ebic.gamma = ebic.gamma,
      top_seed_prop = 0.2,
      expansion_gain = 0.05
    )
    candidate_mp <- mp_res$candidates
    ebic_net <- mp_res$adjacency_matrix
  } else {
    # for future methods...
    stop("Unknown candidate_method: ", candidate_method,
         ". Available methods: 'full', 'clique'")
  }

  if (!is.null(user_multiplets)) {
    valid_user_mp <- Filter(function(x) length(x) %in% k, user_multiplets)
    if (length(valid_user_mp) > 0) {
      candidate_mp <- unique(c(candidate_mp, valid_user_mp))
      message("Added ", length(valid_user_mp), " user-specified multiplets.")
    }
  }

  candidate_mp <- unique(lapply(candidate_mp, sort))
  if (length(candidate_mp) == 0) {
    stop("No valid candidate multiplets generated.")
  }
  message("Total candidate multiplets: ", length(candidate_mp))

  # Step 3: three-step testing pipeline
  test_results <- sig_multiplet_testing_pipeline(
    data = norm_data,
    multiplets = candidate_mp,
    n_perm = n_perm,
    n_boot = n_boot,
    alpha_perm = alpha_perm,
    alpha_boot = alpha_boot,
    Oinfo_threshold = Oinfo_threshold,
    correct_method = correct_method,
    cor_method = cor_method,
    parallel = parallel,
    n_cores = n_cores
  )

  future::plan(
    future::sequential
  )

  # Step 4: hypergraph construction
  message("Building hypergraphs...")
  sig_hyperedges <- test_results[test_results$significant, ]

  full_hypergraph <- lapply(strsplit(sig_hyperedges$vars, ","), as.integer)
  names(full_hypergraph) <- paste0("HE", sig_hyperedges$hyperedge_id)

  redundant_idx <- which(sig_hyperedges$Oinfo > 0)
  redundant_hypergraph <- full_hypergraph[redundant_idx]

  synergistic_idx <- which(sig_hyperedges$Oinfo < 0)
  synergistic_hypergraph <- full_hypergraph[synergistic_idx]

  build_incidence_matrix <- function(
    hyperedges,
    weights = NULL,
    n_nodes
  ){

    H <- matrix(
      0,
      nrow = length(hyperedges),
      ncol = n_nodes
    )

    rownames(H) <- names(hyperedges)

    for(i in seq_along(hyperedges)){

      value <- if(is.null(weights)) 1 else weights[i]

      H[i, hyperedges[[i]]] <- value
    }

    H
  }

  # Incidence matrices
  incidence_matrix <- build_incidence_matrix(
    full_hypergraph,
    n_nodes = n_vars
  )

  redundant_incidence_matrix <- build_incidence_matrix(
    redundant_hypergraph,
    n_nodes = n_vars
  )

  synergistic_incidence_matrix <- build_incidence_matrix(
    synergistic_hypergraph,
    n_nodes = n_vars
  )

  # Step 5: hypergraph object
  result <- list(
    hypergraph = full_hypergraph,
    redundant_hypergraph = redundant_hypergraph,
    synergistic_hypergraph = synergistic_hypergraph,

    incidence_matrix = incidence_matrix,
    redundant_incidence_matrix = redundant_incidence_matrix,
    synergistic_incidence_matrix = synergistic_incidence_matrix,

    hyperedge_details = sig_hyperedges,
    all_test_results = test_results,

    raw_data = data,
    cor_matrix = cor_mat,
    normalized_data = norm_data,
    ebicglasso_network = ebic_net,
    candidate_multiplets = candidate_mp,
    candidate_method_used = candidate_method,

    parameters = params,
    n_subjects = nrow(norm_data),
    n_variables = n_vars,
    n_hyperedges = length(full_hypergraph),
    n_redundant = length(redundant_hypergraph),
    n_synergistic = length(synergistic_hypergraph)
  )

  class(result) <- "psychHypergraph"

  return(result)
}


#' Print method for psychHypergraph objects
#'
#' @param x A psychHypergraph object
#' @param ... Additional arguments
#'
#' @export
print.psychHypergraph <- function(x, ...) {
  cat("psychHypergraph Object Summary\n")
  cat("==============================\n\n")
  cat(sprintf("Subjects: %d\n", x$n_subjects))
  cat(sprintf("Variables: %d\n", x$n_variables))
  cat(sprintf("Candidate Multiplets Evaluated: %d\n", nrow(x$all_test_results)))
  cat(sprintf("Total Significant Hyperedges: %d\n", x$n_hyperedges))
  cat(sprintf("  - Redundant (O-info > 0): %d\n", x$n_redundant))
  cat(sprintf("  - Synergistic (O-info < 0): %d\n", x$n_synergistic))
  cat("\nO-Info Statistics:\n")
  cat(sprintf("  Range: [%.3f, %.3f]\n",
              min(x$hyperedge_details$Oinfo),
              max(x$hyperedge_details$Oinfo)))
  cat(sprintf("  Mean (abs): %.3f\n", mean(abs(x$hyperedge_details$Oinfo))))
  cat("\nParameters:\n")
  cat(sprintf("  cor_method: %s\n", x$parameters$cor_method))
  cat(sprintf("  k: %s\n", paste(x$parameters$k, collapse = ", ")))
  cat(sprintf("  Candidate selection method: %s\n", x$parameters$candidate_method))
  cat("\nUse $hypergraph to access full hyperedge list\n")
  cat("Use $hyperedge_details to view complete test results\n")
  cat("\nTo plot:\n")
  cat("  plot(x)                    # Full hypergraph\n")
  cat("  plot(x, type='redundant')  # Only redundant\n")
  cat("  plot(x, type='synergistic')# Only synergistic\n")
  cat("  plot(x, top_n=10)          # Top 10 by |O-info|\n")
}

#' Summary method for psychHypergraph objects
#'
#' @param object A psychHypergraph object
#' @param ... Additional arguments
#' @export
summary.psychHypergraph <- function(object, ...) {
  print(object)
}
