# Section 4.3: reproducible base-R checks, without private fitted_rates.csv.
# Run at repository root:
# Rscript 04_results/4_3_multi_population/selftest_4_3.R
# These tests validate the code, NOT the full GBD fit or paper Figure 9.

args <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(args)) stop("Run selftest_4_3.R using Rscript, not source().")
here <- dirname(normalizePath(sub("^--file=", "", args[1]), mustWork = TRUE))
root <- normalizePath(file.path(here, "../.."), mustWork = TRUE)
source(file.path(root, "04_results/4_1_single_population/model_functions.R"),
       encoding = "UTF-8")
source(file.path(here, "multi_model_functions.R"), encoding = "UTF-8")
source(file.path(here, "author_appendix_inventory.R"), encoding = "UTF-8")

# §4.2 now requires genuine seven-CVI review and no longer stores any
# historical Silhouette-only cluster CSV. The REAL 126-row grouping contract
# is checked by run_section43(source42=...) when reviewed results exist.
cat("PASS: §4.3 synthetic test does not import obsolete §4.2 groupings\\n")

# Synthetic population: J=2 different causes, six age points, 36 years.
# Explicitly constructed WITHOUT reading or inventing any real GBD death counts.
A <- length(AGES)
TT <- length(YEARS)
J <- 2L
tt <- seq_len(TT)
v <- (mean(AGES) - AGES) / 10

synthetic_case <- function(model) {
  is_lc1 <- endsWith(model, "LC1")
  is_lc2 <- endsWith(model, "LC2")
  is_cbd <- model == "CBD"
  is_plat <- model %in% c("Plat2", "Plat3")
  is_binom <- startsWith(model, "Binomial") || is_cbd
  nf <- if (model == "Plat3") 3L else if (is_lc2 || is_cbd || model == "Plat2") 2L else 1L
  D <- E <- array(0, c(A, TT, J))
  singles <- vector("list", J)
  names(singles) <- c("synthetic_cause_A", "synthetic_cause_B")
  common1 <- 0.27 * sin(tt / 7) - 0.008 * tt
  common2 <- 0.16 * cos(tt / 9) + 0.006 * tt
  for (j in seq_len(J)) {
    if (is_lc1) {
      B <- matrix(0.85 + 0.08 * seq_len(A) + 0.015 * j, A, 1L)
      K <- matrix(common1 + 0.015 * j * sin(tt / 5), TT, 1L)
    } else if (is_lc2) {
      B <- cbind(0.8 + 0.07 * seq_len(A),
                 (seq_len(A) - 3.5) / 3 + 0.04 * j)
      K <- cbind(0.24 * sin(tt / (5 + j)), common2)
    } else if (is_cbd) {
      B <- cbind(rep(1, A), -v)
      K <- cbind(-6.2 + common1, 0.10 * cos(tt / (7 + j)))
    } else if (is_plat) {
      B <- cbind(rep(1, A), v)
      if (nf == 3L) B <- cbind(B, pmax(v, 0))
      K <- cbind(common1, 0.09 * cos(tt / (8 + j)))
      if (nf == 3L) K <- cbind(K, 0.07 * sin(tt / (6 + j)))
    } else stop("Unsupported synthetic model")
    a <- if (is_cbd) rep(0, A) else -7.0 + 0.24 * seq_len(A) + 0.13 * j
    singles[[j]] <- list(a = a, B = B, K = K)
    eta <- a + B %*% t(K)
    E[,,j] <- outer(seq_len(A), tt, function(x, y)
      80000 + 1300 * x + 200 * y + j * 10000)
    D[,,j] <- pmax(1, round(E[,,j] *
      if (is_binom) plogis(eta) else exp(eta)))
  }
  list(D = D, E = E, singles = singles)
}

results <- list()
for (model in MODELS) {
  z <- synthetic_case(model)
  p <- make_multi_problem(z$D, z$E, model, z$singles)
  stopifnot(p$J == 2L, length(p$start) > 0L,
            all(is.finite(p$start)))
  c0 <- p$calc(p$start, return_details = TRUE)
  stopifnot(is.finite(c0$value), all(is.finite(c0$gradient)),
            identical(dim(c0$mu), dim(z$D)), all(c0$mu > 0))
  if (endsWith(model,"LC1") || endsWith(model,"LC2")) {
    # SVD column sign has no scientific meaning: flipping B and K together
    # must not change the initial joint eta, even when sign alignment runs.
    zflip <- z
    h <- if (endsWith(model,"LC2")) 2L else 1L
    zflip$singles[[2]]$B[,h] <- -zflip$singles[[2]]$B[,h]
    zflip$singles[[2]]$K[,h] <- -zflip$singles[[2]]$K[,h]
    pf <- make_multi_problem(z$D,z$E,model,zflip$singles)
    ef <- pf$calc(pf$start, return_details=TRUE)$eta
    stopifnot(max(abs(ef-c0$eta)) < 1e-8)
    if (endsWith(model,"LC2")) {
      # Here each synthetic cause has exactly the SAME second index:
      # the constrained rotation and pooling must reproduce both original
      # single-population fitted surfaces up to numerical rounding.
      for (j in seq_len(J)) {
        sg <- z$singles[[j]]
        independent_eta <- sg$a + sg$B %*% t(sg$K)
        stopifnot(max(abs(independent_eta-c0$eta[,,j])) < 1e-8)
      }
    }
  }
  if (model != "CBD") {
    stopifnot(abs(mean(c0$u$ks)) < 1e-8)
    if (!is.null(c0$u$Kc)) {
      for (h in seq_len(dim(c0$u$Kc)[3]))
        stopifnot(max(abs(colMeans(c0$u$Kc[,,h]))) < 1e-8)
    }
  }
  if (endsWith(model, "LC1") || endsWith(model, "LC2")) {
    stopifnot(max(abs(colSums(c0$u$B1) - 1)) < 1e-10)
    if (endsWith(model, "LC2"))
      stopifnot(max(abs(colSums(c0$u$B2) - 1)) < 1e-10)
  }
  if (!p$binom) {
    # Stress the Poisson overflow branch: score must remain finite and NONZERO.
    huge <- p$start
    huge[1] <- huge[1] + 55
    stress <- p$calc(huge)
    stopifnot(is.finite(stress$value),
              all(is.finite(stress$gradient)),
              max(abs(stress$gradient)) > 0)
  }
  check <- check_multi_gradient(p, eps = 2e-6, ncheck = 60L)
  stopifnot(isTRUE(check$gradient_check_ok))
  # Also test a random direction away from the initialization.
  theta <- p$start + 0.0001 * sin(seq_along(p$start))
  direction <- cos(seq_along(theta))
  direction <- direction / sqrt(sum(direction^2))
  step <- 1e-4
  fd <- (p$calc(theta + step * direction)$value -
         p$calc(theta - step * direction)$value) / (2 * step)
  exact <- sum(p$calc(theta)$gradient * direction)
  dir_error <- abs(fd - exact) / max(1, abs(fd), abs(exact))
  stopifnot(is.finite(dir_error), dir_error < 1e-4)
  results[[model]] <- c(checked = check$parameters_checked,
                        score_error = check$maximum_relative_gradient_error,
                        directional_error = dir_error)
  cat("PASS:", model, "| numeric score error = ",
      format(check$maximum_relative_gradient_error, digits = 4), "\n")
}

# Run the actual optimizer on a small synthetic, correctly specified example.
for (model in c("Poisson LC1", "Plat3")) {
  z <- synthetic_case(model)
  p <- make_multi_problem(z$D, z$E, model, z$singles)
  before <- p$calc(p$start)$value
  f <- fit_multi_group(p, z$D, z$E, maxit = 180L,
                       restarts = 1L, tol = 5e-3)
  stopifnot(is.finite(f$metrics$objective),
            f$metrics$objective <= before + 1e-6,
            identical(dim(f$rate), dim(z$D)),
            all(is.finite(f$rate)), all(f$rate > 0))
  cat("PASS: end-to-end synthetic optimizer for", model, "\n")
}

# Appendix source filenames must match the 52 known published image names.
tmp <- tempfile(fileext = ".csv")
items <- author_appendix_inventory(tmp)
stopifnot(nrow(items) == 52L, sum(items$image_type == "joint_parameters") == 26L,
          sum(items$image_type == "indices_comparison") == 26L,
          length(unique(items$source_url)) == 52L)
unlink(tmp)
cat("PASS: Appendix 9 inventory = 52 unique reference images\n")
cat("ALL SECTION 4.3 SYNTHETIC TESTS PASSED; REAL GBD FIT STILL REQUIRES INPUT CSV.\n")