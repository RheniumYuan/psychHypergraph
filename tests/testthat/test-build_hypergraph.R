test_that("build_hypergraph works", {
  library(psych)
  data(bfi)
  data_bfi <- bfi[, 1:10]
  hg_bfi <- build_hypergraph(
    data = data_bfi,
    cor_method = "pearson",
    k = 3:5,
    candidate_method = "clique",
    n_perm = 20,
    n_boot = 20,
    seed = 123,
    parallel = FALSE,
    n_cores = 1
  )

  expect_s3_class(
    hg_bfi,
    "psychHypergraph"
  )
})
