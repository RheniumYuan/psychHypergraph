#' Rank-based Normalization
#'
#' Converts psychological ordinal/non-normal data to standard normal using inverse Gaussian CDF
#' @param data Data frame of psychological scale data
#' @param cols Columns to normalize (default: all numeric columns)
#' @param ties.method Method for handling ties
#' @param eps Numerical stability bound
#' @return Normalized data frame
#' @noRd
rank_normalize <- function(data, cols = NULL, ties.method = "average", eps = 1e-10) {
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }

  # Auto-select numeric columns if not specified
  if (is.null(cols)) {
    cols <- names(data)[apply(data, 2, is.numeric)]
  }

  # Check columns exist
  if (!all(cols %in% names(data))) {
    stop("Some specified columns are not present in the data.")
  }

  df_norm <- data

  # Normalize each column
  for (col in cols) {
    x <- data[[col]]
    if (!is.numeric(x)) next

    # Compute ranks
    r <- rank(x, ties.method = ties.method, na.last = "keep")
    n <- sum(!is.na(x))

    # Avoid division by zero
    if (n == 0) next

    # Rank to uniform probability
    p <- (r - 0.5) / n

    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)

    df_norm[[col]] <- stats::qnorm(p)
  }

  return(df_norm)
}

#' Compute Correlation Matrix for Psychological Data
#'
#' Supports Pearson and polychoric correlations
#' @param data Numeric data frame
#' @param method Correlation method: pearson, polychoric
#' @param cols Columns to use (default: all numeric)
#' @return Correlation matrix
#' @noRd
compute_correlation <- function(
    data,
    method = c("pearson", "polychoric"),
    cols = NULL
) {
  method <- match.arg(method)

  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }

  # Select numeric columns
  if (is.null(cols)) {
    cols <- names(data)[apply(data, 2, is.numeric)]
  }
  df <- data[, cols, drop = FALSE]

  # Compute correlation based on method
  if (method == "pearson") {
    cor_mat <- stats::cor(df, use = "complete.obs")
  } else if (method == "polychoric") {
    cor_list <- psych::polychoric(df)
    cor_mat <- cor_list$rho
  }

  # Ensure positive definiteness
  eig <- eigen(cor_mat,
               symmetric = TRUE)$values

  if(min(eig) < -1e-8){
    warning(
      "Correlation matrix is not positive definite."
    )
  }

  cor_mat <- as.matrix(cor_mat)

  return(cor_mat)
}

#' Preprocessing Pipeline for Psychometric Hypergraph Analysis
#'
#' A comprehensive pipeline that handles missing data, applies rank normalization, and computes the correlation matrix for psychological scale data. Designed to prepare data for hypergraph analysis.
#' @param data Raw data frame of psychometric data
#' @param cor_method Correlation method: pearson, polychoric
#' @param handle_na "impute" or "omit"
#' @param impute_fun Imputation function
#' @return List with normalized data and correlation matrix
#'
#' @export
data_preprocess_pipeline <- function(
    data,
    cor_method = c("pearson", "polychoric"),
    handle_na = c("impute", "omit"),
    impute_fun = stats::median
) {
  cor_method <- match.arg(cor_method)
  handle_na <- match.arg(handle_na)

  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }

  # Keep only numeric columns
  df <- data[, apply(data, 2, is.numeric)]

  # Missing value handling
  if (handle_na == "omit") {
    df <- stats::na.omit(df)
  } else if (handle_na == "impute") {
    df <- as.data.frame(lapply(df, function(x) {
      x[is.na(x)] <- impute_fun(x, na.rm = TRUE)
      x
    }))
  }

  # Step 1: Rank normalization
  if (cor_method == "pearson") {
    df_proc <- rank_normalize(df)
  } else {
    # polychoric should use raw ordinal data
    df_proc <- df
  }

  # Step 2: Compute correlation matrix
  cor_mat <- compute_correlation(df_proc, method = cor_method)

  # Return pipeline result
  return(list(
    normalized_data = df_proc,
    cor_matrix = cor_mat,
    cor_method = cor_method,
    n_obs = nrow(df)
  ))
}
