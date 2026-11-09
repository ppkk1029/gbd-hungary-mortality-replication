# 03 | Varga (2025) §4.1 / Figure 3: one-factor Poisson Lee-Carter
# GBD 2023 Hungary: METHOD replication, not original HCSO numerical replication
#
# Run in RStudio: open this file -> Source
# Or from Console: source(file.choose())  [first choose this run.R]
# Then: choose Population CSV, then choose Deaths CSV, in that order.
# No fixed input paths. Preview in Console + Plots BEFORE choosing output folder.
#
# This stage: 9 GBD Level-2 causes selected from 1980-2015 deaths;
# fit 1980-2015, hold out 2016-2019, set aside pandemic 2020-2023.
# Original paper Figure 3 fitted on a longer FULL historical sample.
# Here GBD deaths are modeled real-valued estimates, not raw count data:
# Poisson is a MEAN/QUASI-LIKELIHOOD structure, NOT a literal count likelihood.
# Requires BASE R only (no extra packages).

cat("\n[1/2] 请选择 GBD Population CSV ...\n")
population_path <- file.choose()
cat("[2/2] 请选择 GBD Deaths CSV (Number + Rate) ...\n")
deaths_path <- file.choose()
if (normalizePath(population_path, mustWork = TRUE) ==
    normalizePath(deaths_path, mustWork = TRUE)) {
  stop("Population 与 Deaths 不能是同一个文件。")
}
pop_raw <- read.csv(population_path, check.names = FALSE, stringsAsFactors = FALSE)
deaths_raw <- read.csv(deaths_path, check.names = FALSE, stringsAsFactors = FALSE)

years <- 1980:2023
fit_years <- 1980:2015
sexes <- c("Male", "Female")
ages <- c("60-64 years", "65-69 years", "70-74 years",
          "75-79 years", "80-84 years", "85-89 years")
age_mid <- c(62, 67, 72, 77, 82, 87)
A <- length(ages)
Tn <- length(fit_years)

# ---- 1. Verify Population and Deaths scope ----
pf <- c("year", "location_name", "sex_name", "age_name",
        "measure_name", "metric_name", "val")
df <- c("year", "location_name", "sex_name", "age_name",
        "cause_name", "measure_name", "metric_name", "val")
if (!all(pf %in% names(pop_raw))) stop("Population CSV 列缺失。")
if (!all(df %in% names(deaths_raw))) stop("Deaths CSV 列缺失。")
if (anyNA(pop_raw[pf]) || nrow(pop_raw) != 44L * 2L * 6L ||
    !all(pop_raw$location_name == "Hungary") ||
    !all(pop_raw$measure_name == "Population") ||
    !all(pop_raw$metric_name == "Number") ||
    !setequal(unique(pop_raw$year), years) ||
    !setequal(unique(pop_raw$sex_name), sexes) ||
    !setequal(unique(pop_raw$age_name), ages) ||
    anyDuplicated(pop_raw[c("year", "sex_name", "age_name")]) ||
    any(!is.finite(pop_raw$val) | pop_raw$val <= 0)) {
  stop("Population 数据范围、主键或数值未通过检查。")
}
if ("population_group_name" %in% names(pop_raw) &&
    !all(pop_raw$population_group_name == "All Population")) {
  stop("Population 的 population_group_name 必须为 All Population。")
}
if (anyNA(deaths_raw[df]) ||
    !all(deaths_raw$location_name == "Hungary") ||
    !all(deaths_raw$measure_name == "Deaths") ||
    !setequal(unique(deaths_raw$year), years) ||
    !setequal(unique(deaths_raw$sex_name), sexes) ||
    !setequal(unique(deaths_raw$age_name), ages) ||
    !setequal(unique(deaths_raw$metric_name), c("Number", "Rate")) ||
    length(unique(deaths_raw$cause_name)) != 19L ||
    nrow(deaths_raw) != 44L * 2L * 6L * 19L * 2L ||
    anyDuplicated(deaths_raw[c("year", "sex_name", "age_name",
                               "cause_name", "metric_name")]) ||
    any(!is.finite(deaths_raw$val) | deaths_raw$val < 0)) {
  stop("Deaths 数据范围、主键或数值未通过检查。")
}

pop <- pop_raw[c("year", "sex_name", "age_name", "val")]
names(pop)[4] <- "exposure_proxy"
keys <- c("year", "sex_name", "age_name", "cause_name")
num <- deaths_raw[deaths_raw$metric_name == "Number", c(keys, "val")]
rate <- deaths_raw[deaths_raw$metric_name == "Rate", c(keys, "val")]
names(num)[5] <- "deaths"
names(rate)[5] <- "rate_per_100k"
data0 <- merge(num, rate, by = keys, all = TRUE, sort = FALSE)
data0 <- merge(data0, pop, by = c("year", "sex_name", "age_name"),
               all.x = TRUE, sort = FALSE)
if (nrow(data0) != 44L * 2L * 6L * 19L ||
    anyNA(data0[c("deaths", "rate_per_100k", "exposure_proxy")])) {
  stop("Deaths / Population 没有正确配对。")
}
data0$rate_from_number <- data0$deaths / data0$exposure_proxy * 100000
data0$relative_error <- ifelse(
  data0$rate_per_100k > 0,
  abs(data0$rate_from_number - data0$rate_per_100k) /
    data0$rate_per_100k,
  ifelse(data0$rate_from_number == 0, 0, Inf)
)
if (!all(is.finite(data0$relative_error))) {
  stop("Number/Rate 出现不一致的零单元，不能继续。")
}
qc_med <- median(data0$relative_error)
qc_p95 <- as.numeric(quantile(data0$relative_error, 0.95))
qc_max <- max(data0$relative_error)
cat("\n人口与死亡数据读取成功：Population ", nrow(pop_raw),
    " 行；Deaths ", nrow(deaths_raw), " 行。\n", sep = "")
cat("Number/Rate/Population QC（相对误差）：median=", qc_med,
    "  p95=", qc_p95, "  max=", qc_max, "\n")
if (qc_med > 0.02 || qc_p95 > 0.10) {
  print(head(data0[order(-data0$relative_error), ], 10))
  stop("死亡人数、Rate 与人口分母的一致性筛查未通过。")
}

# Same 9 GBD Level-2 causes as Figure 2: training-period deaths.
r <- aggregate(deaths ~ cause_name,
               data = subset(data0, year %in% fit_years), FUN = sum)
r <- r[order(-r$deaths, r$cause_name), ]
chosen <- head(r$cause_name, 9L)
cat("\n九类 GBD Level-2 死因：\n")
print(data.frame(rank = 1:9, cause_name = chosen), row.names = FALSE)

# ---- 2. One-factor Poisson Lee-Carter by numerical optimization ----
# D_xt = E_xt * exp(alpha_x + beta_x * kappa_t)
# Identification: max_x abs(beta_x)=1; mean_t kappa_t=0 (stable alternative)
# We minimize HALF Poisson deviance (ignoring count factorial terms) using
# continuous GBD modeled deaths. The optimization is pseudo likelihood.
#
# Numerical parameterization: fix an SVD-selected reference beta to 1, center kappa.
# After convergence normalize beta to max(abs(beta))=1 while preserving fit.

fit_lc <- function(D, E) {
  stopifnot(is.matrix(D), is.matrix(E), all(dim(D) == c(A, Tn)),
            all(dim(E) == c(A, Tn)), all(is.finite(D)),
            all(is.finite(E)), all(D >= 0), all(E > 0), sum(D) > 0)
  # A loadings vector with near-zero sum makes sum(beta)=1 ill-conditioned.
  # Instead fix beta at the largest-|SVD-loading| reference age to 1
  # during optimization; center kappa; AFTER optimization rescale so
  # max_x abs(beta_x)=1. This DOES NOT change beta*kappa or fitted rates.
  # It is an alternative identification, NOT a change in mean model.
  log_r <- log((D + 0.5) / E)
  alpha0 <- rowMeans(log_r)
  centered <- sweep(log_r, 1L, alpha0, "-")
  sv <- svd(centered, nu = 1L, nv = 1L)
  if (!is.finite(sv$d[1]) || sv$d[1] < 1e-7) {
    stop("训练数据中几乎没有可识别的时期变化。")
  }
  u <- sv$u[, 1L]
  anchor <- which.max(abs(u))
  beta0 <- u / u[anchor]  # beta[anchor] = 1, ensures stable initial scale
  kappa0 <- as.numeric(sv$d[1L] * sv$v[, 1L] * u[anchor])
  kbar <- mean(kappa0)
  alpha0 <- alpha0 + beta0 * kbar
  kappa0 <- kappa0 - kbar

  initial_par <- c(alpha0, beta0[-anchor], kappa0[-Tn])
  nA <- A - 1L
  nT <- Tn - 1L
  logE <- log(E)
  totD <- sum(D)
  nz <- D > 0

  unpack <- function(theta) {
    alpha <- theta[seq_len(A)]
    beta <- numeric(A)
    beta[anchor] <- 1
    beta[-anchor] <- theta[A + seq_len(nA)]
    kappa <- c(theta[A + nA + seq_len(nT)],
               -sum(theta[A + nA + seq_len(nT)]))
    list(alpha = alpha, beta = beta, kappa = kappa)
  }
  calc <- function(theta) {
    q <- unpack(theta)
    eta <- sweep(outer(q$beta, q$kappa), 1, q$alpha, "+")
    lp <- logE + eta
    if (any(!is.finite(lp)) || max(lp) > 50) return(NULL)
    lambda <- exp(lp)
    if (any(!is.finite(lambda))) return(NULL)
    residual <- lambda - D
    dev <- 2 * (sum(residual) +
                sum(D[nz] * (log(D[nz]) - lp[nz])))
    # Scale objective & analytic gradient by total GBD estimated deaths
    ga <- rowSums(residual)
    gb <- as.numeric(residual %*% q$kappa)
    gk <- as.numeric(crossprod(q$beta, residual))
    grad <- c(ga, gb[-anchor], gk[-Tn] - gk[Tn]) / totD
    list(value = dev / (2 * totD), deviance = dev, grad = grad,
         eta = eta, lambda = lambda, alpha = q$alpha,
         beta = q$beta, kappa = q$kappa)
  }
  fn <- function(theta) {
    z <- calc(theta)
    if (is.null(z)) return(1e100)
    z$value
  }
  gr <- function(theta) {
    z <- calc(theta)
    if (is.null(z)) return(rep(0, length(theta)))
    z$grad
  }

  set.seed(20261109)
  starts <- list(initial_par)
  for (i in 1:2) {
    v <- initial_par
    # Independent, small jitter to identify accidental local minima.
    v[seq_len(A)] <- v[seq_len(A)] + rnorm(A, sd = 0.015)
    v[A + seq_len(nA)] <- v[A + seq_len(nA)] + rnorm(nA, sd = 0.04)
    v[A + nA + seq_len(nT)] <- v[A + nA + seq_len(nT)] +
      rnorm(nT, sd = max(sd(kappa0), 0.01) * 0.03)
    starts[[i+1L]] <- v
  }
  attempts <- lapply(starts, function(theta0) {
    fit <- tryCatch(
      optim(theta0, fn = fn, gr = gr, method = "BFGS",
            control = list(maxit = 6000L, reltol = 1e-11)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    fitted <- calc(fit$par)
    if (is.null(fitted)) return(NULL)
    list(fitted = fitted, code = fit$convergence,
         max_grad = max(abs(fitted$grad)), value = fitted$value)
  })
  valid <- Filter(function(v)
    !is.null(v) && v$code == 0L &&
    is.finite(v$value) && v$max_grad < 0.02, attempts)
  if (length(valid) < 1L) {
    stop("Poisson LC 未找到达到收敛/梯度要求的解：建议检查数据或改用更稳健算法。")
  }
  best <- valid[[which.min(vapply(valid, function(v) v$value, numeric(1)))]]
  q <- best$fitted
  # Stable identification: max(abs(beta)) = 1.
  beta_scale <- max(abs(q$beta))
  beta <- q$beta / beta_scale
  kappa <- q$kappa * beta_scale
  # Consistent sign across causes: prefer positive sum of beta.
  if (sum(beta) < 0) {
    beta <- -beta
    kappa <- -kappa
  }
  if (max(abs(beta)) > 1.00000001 ||
      abs(mean(kappa)) > 1e-7) stop("归一化/中心化失败。")
  # Explicitly verify reconstructed mortality surface is unchanged.
  eta_new <- sweep(outer(beta, kappa), 1, q$alpha, "+")
  reconstruction_error <- max(abs(eta_new - q$eta))
  if (reconstruction_error > 1e-8) {
    stop("参数尺度重识别改变了拟合死亡率，不允许继续。")
  }
  ratio <- abs(sum(beta)) / sum(abs(beta))
  # An almost-cancelled beta sum indicates that the original
  # sum(beta)=1 identification would explode; flag for later DTW.
  flag <- if (ratio < 0.08) "REVIEW_near_zero_beta_sum" else "OK"
  list(alpha = q$alpha, beta = beta, kappa = kappa,
       rate_fitted = exp(q$eta), count_fitted = q$lambda,
       deviance = q$deviance,
       objective_final = q$value,
       convergence = best$code,
       max_abs_gradient = best$max_grad,
       n_converged_starts = length(valid),
       beta_cancellation_ratio = ratio,
       beta_flag = flag,
       kappa_sd = sd(kappa),
       beta_scale_raw = beta_scale,
       reconstruction_error = reconstruction_error)
}

# We keep all data in memory until user has seen the plots and approved saving.
results <- list()
diagnostics <- data.frame()
alpha_table <- beta_table <- kappa_table <- fitted_table <- data.frame()
for (sx in sexes) {
  for (cause in chosen) {
    dat <- data0[data0$sex_name == sx &
                   data0$cause_name == cause &
                   data0$year %in% fit_years, ]
    dat <- dat[order(match(dat$age_name, ages),
                     match(dat$year, fit_years)), ]
    if (nrow(dat) != A * Tn ||
        !identical(rep(ages, each = Tn), dat$age_name) ||
        !identical(rep(fit_years, times = A), dat$year)) {
      stop("拟合年龄×年份网格不完整：", sx, ", ", cause)
    }
    D <- matrix(dat$deaths, nrow = A, ncol = Tn, byrow = TRUE)
    E <- matrix(dat$exposure_proxy, nrow = A, ncol = Tn, byrow = TRUE)
    cat("\n正在拟合：", sx, " | ", cause, "\n", sep = "")
    m <- fit_lc(D, E)
    key <- paste(sx, cause, sep = "__")
    results[[key]] <- m
    diagnostics <- rbind(diagnostics, data.frame(
      sex_name = sx, cause_name = cause,
      convergence_code = m$convergence,
      half_deviance_per_death = m$objective_final,
      poisson_deviance = m$deviance,
      max_abs_gradient = m$max_abs_gradient,
      n_converged_starts = m$n_converged_starts,
      beta_cancellation_ratio = m$beta_cancellation_ratio,
      beta_flag = m$beta_flag,
      kappa_sd = m$kappa_sd,
      reconstruction_error = m$reconstruction_error,
      n_years = Tn, n_ages = A))
    alpha_table <- rbind(alpha_table, data.frame(
      sex_name = sx, cause_name = cause, age_name = ages,
      age_midpoint = age_mid, alpha = m$alpha))
    beta_table <- rbind(beta_table, data.frame(
      sex_name = sx, cause_name = cause, age_name = ages,
      age_midpoint = age_mid, beta = m$beta))
    kappa_table <- rbind(kappa_table, data.frame(
      sex_name = sx, cause_name = cause, year = fit_years,
      kappa = m$kappa))
    fitted_table <- rbind(fitted_table, data.frame(
      sex_name = sx, cause_name = cause,
      age_name = rep(ages, each = Tn),
      year = rep(fit_years, times = A),
      deaths_gbd = as.vector(t(D)),
      exposure_proxy = as.vector(t(E)),
      fitted_count = as.vector(t(m$count_fitted)),
      observed_rate_per_100k = as.vector(t(D / E * 100000)),
      fitted_rate_per_100k = as.vector(t(m$rate_fitted * 100000))))
  }
}
if (any(diagnostics$beta_flag != "OK")) {
  cat("\n警告：以下年龄敏感度正负抵消明显，不宜在未核查前解释或用于 DTW。\n")
  print(diagnostics[diagnostics$beta_flag != "OK",
                    c("sex_name", "cause_name", "beta_cancellation_ratio",
                      "kappa_sd", "beta_flag")], row.names = FALSE)
}
if (any(diagnostics$convergence_code != 0L)) {
  print(diagnostics[diagnostics$convergence_code != 0L, ])
  stop("至少一个模型未达到优化器收敛判据；不保存正式结果。")
}

# ---- 3. Preview parameter table and Figure 3 in RStudio ----
cat("\n========================================================\n")
cat("PREVIEW ONLY — §4.1 SINGLE-POPULATION. NO OUTPUT FILES SAVED.\n")
cat("========================================================\n")
cat("\n模型拟合汇总（检查 convergence_code 都为 0）：\n")
print(diagnostics, row.names = FALSE)
cat("\n示例：男性心血管疾病（如包含）的一维时间因子前 10 年：\n")
print(head(kappa_table[kappa_table$sex_name == "Male" &
                          kappa_table$cause_name == chosen[1], ], 10),
      row.names = FALSE)
cat("\n稳健识别：每个死因/性别 max(abs(beta))=1，mean(kappa)=0。\n")
for (sx in sexes) {
  for (cause in chosen) {
    m <- results[[paste(sx, cause, sep = "__")]]
    if (abs(max(abs(m$beta)) - 1) > 1e-8 ||
        abs(mean(m$kappa)) > 1e-8) {
      stop("模型识别约束错误。")
    }
  }
}
if (interactive()) {
  View(diagnostics, title = "Poisson LC fit diagnostics")
  View(kappa_table, title = "Poisson LC kappa period factors")
  View(beta_table, title = "Poisson LC beta age sensitivities")
}

draw_figure3 <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    layout(1)
    par(oldpar)
  })
  # Six parameter panels and a dedicated, full-width legend at the bottom.
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7, 7),
                nrow = 3, byrow = TRUE),
         heights = c(1, 1, 0.36))
  par(mar = c(3.5, 4.0, 2.9, 0.8), oma = c(0.2, 0.2, 2.2, 0.2))
  cols <- grDevices::hcl.colors(9, "Dark 3")
  for (sx in sexes) {
    for (parameter in c("alpha", "beta", "kappa")) {
      curves <- lapply(chosen, function(cause)
        results[[paste(sx, cause, sep = "__")]][[parameter]])
      X <- if (parameter == "kappa") fit_years else age_mid
      Y <- range(unlist(curves), finite = TRUE)
      if (diff(Y) == 0) Y <- Y + c(-1, 1) * 1e-5
      plot(X, curves[[1]], type = "n", ylim = Y, xlab =
             if (parameter == "kappa") "Year" else "Age midpoint",
           ylab = parameter, main = paste(sx, parameter), cex.main = 1.0)
      for (i in seq_along(chosen))
        lines(X, curves[[i]], col = cols[i], lwd = 1.55)
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = chosen, col = cols, lty = 1, lwd = 1.8,
         cex = 0.72, ncol = 3, bty = "n")
  mtext("GBD Hungary | Poisson Lee-Carter (Level 2 causes)",
        side = 3, outer = TRUE, line = 0.5, cex = 1.15, font = 2)
}
cat("\n请查看 RStudio Plots 窗格 §4.1 的六面板参数图：alpha、beta、kappa × 男女。\n")
draw_figure3()

answer <- tolower(trimws(readline(
  "\n检查预览结果后，是否保存 Figure 3 结果？输入 y 保存，其他键取消："
)))
if (!identical(answer, "y")) {
  cat("未保存文件。对象仍在 R 会话内存（diagnostics、kappa_table 等）。\n")
} else {
  if (.Platform$OS.type == "windows") {
    output_parent <- utils::choose.dir(default = getwd(),
                   caption = "请选择 Figure 3 结果的保存位置")
  } else if (requireNamespace("tcltk", quietly = TRUE)) {
    output_parent <- tcltk::tk_choose.dir(default = getwd(),
                   caption = "请选择 Figure 3 结果的保存位置")
  } else {
    stop("该系统缺少文件夹选择器，请在桌面 RStudio 中操作。")
  }
  if (is.na(output_parent) || !nzchar(output_parent)) {
    cat("取消选择保存位置，没有保存文件。\n")
  } else {
    outdir <- file.path(output_parent, "Figure3_Poisson_LC_results")
    ok <- TRUE
    if (dir.exists(outdir)) {
      ok <- identical(tolower(trimws(readline(
        "目标文件夹已存在，输入 y 允许覆盖同名文件；其他键取消："
      ))), "y")
    }
    if (!ok) {
      cat("用户取消覆盖，没有保存文件。\n")
    } else {
      if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
      write.csv(alpha_table, file.path(outdir, "alpha.csv"), row.names = FALSE)
      write.csv(beta_table, file.path(outdir, "beta.csv"), row.names = FALSE)
      write.csv(kappa_table, file.path(outdir, "kappa.csv"), row.names = FALSE)
      write.csv(diagnostics, file.path(outdir, "fit_diagnostics.csv"),
                row.names = FALSE)
      write.csv(fitted_table, file.path(outdir, "fitted_rates.csv"),
                row.names = FALSE)
      write.csv(data.frame(rank = 1:9, cause_name = chosen),
                file.path(outdir, "selected_causes.csv"), row.names = FALSE)
      grDevices::png(file.path(outdir, "figure3_poisson_lc.png"),
                     width = 2600, height = 1550, res = 160)
      draw_figure3()
      grDevices::dev.off()
      writeLines(c("Figure 3 GBD Poisson Lee-Carter method replication",
                   "Note: GBD modeled deaths; quasi-likelihood, not actual integer count MLE.",
                   "Fit: 1980-2015. Test reserved: 2016-2019. 2020-2023 excluded.",
                   "Age: 60-89 grouped in six five-year bins.",
                   "GBD Level-2 nine causes selected using 1980-2015 deaths.",
                   "Sexes fitted independently. Population serves as exposure proxy.",
                   "Identification max(abs(beta))=1, mean(kappa)=0 (alternative stable scaling).",
                   "beta_flag=REVIEW warns against uncritical clustering for near-zero beta sum.",
                   paste("Run date:", Sys.time()),
                   paste("Population:", basename(population_path)),
                   paste("Deaths:", basename(deaths_path)),
                   paste("QC median relative error:", qc_med),
                   paste("QC p95 relative error:", qc_p95)),
                 file.path(outdir, "run_summary.txt"))
      cat("\nRESULTS SAVED TO:\n",
          normalizePath(outdir, winslash = "/", mustWork = TRUE), "\n")
      print(list.files(outdir))
    }
  }
}
