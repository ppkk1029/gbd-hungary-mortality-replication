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
# Identification: sum_x beta_x = 1, sum_t kappa_t = 0
# We minimize HALF Poisson deviance (ignoring count factorial terms) using
# continuous GBD modeled deaths. The optimization is pseudo likelihood.
#
# Unknowns: alpha[A], beta_free[A-1], kappa_free[T-1].
# beta[A] = 1 - sum(beta_free), kappa[T] = -sum(kappa_free).
# The following analytic gradient is with respect to FREE parameters.

fit_lc <- function(D, E) {
  stopifnot(is.matrix(D), is.matrix(E), all(dim(D) == c(A, Tn)),
            all(dim(E) == c(A, Tn)), all(is.finite(D)),
            all(is.finite(E)), all(D >= 0), all(E > 0),
            sum(D) > 0)

  log_r <- log((D + 0.5) / E) # Only for INITIALIZATION (not the objective).
  alpha0 <- rowMeans(log_r)
  sv <- svd(sweep(log_r, 1, alpha0, "-"), nu = 1L, nv = 1L)
  u <- sv$u[, 1]
  if (abs(sum(u)) >= 0.15) {
    beta0 <- u / sum(u)
    kappa0 <- as.numeric(sv$d[1] * sv$v[, 1] * sum(u))
  } else {
    beta0 <- rep(1 / A, A)
    kappa0 <- A * colMeans(sweep(log_r, 1, alpha0, "-"))
  }
  kbar <- mean(kappa0)
  alpha0 <- alpha0 + beta0 * kbar
  kappa0 <- kappa0 - kbar
  theta0 <- c(alpha0, beta0[1:(A - 1L)], kappa0[1:(Tn - 1L)])

  unpack <- function(th) {
    a <- th[seq_len(A)]
    bb <- th[A + seq_len(A - 1L)]
    kk <- th[2L * A - 1L + seq_len(Tn - 1L)]
    list(alpha = a, beta = c(bb, 1 - sum(bb)),
         kappa = c(kk, -sum(kk)))
  }
  total_d <- sum(D)
  positives <- D > 0
  log_E <- log(E)
  calc <- function(th) {
    p <- unpack(th)
    eta <- sweep(outer(p$beta, p$kappa), 1L, p$alpha, "+")
    lp <- log_E + eta
    if (any(!is.finite(lp)) || max(lp) > 40) return(NULL)
    L <- exp(lp)
    residual <- L - D
    obj <- (sum(residual) +
              sum(D[positives] * (log(D[positives]) - lp[positives]))) /
           total_d
    ga <- rowSums(residual)
    gb <- as.numeric(residual %*% p$kappa)
    gk <- as.numeric(crossprod(p$beta, residual))
    grad <- c(ga,
              gb[1:(A - 1L)] - gb[A],
              gk[1:(Tn - 1L)] - gk[Tn]) / total_d
    list(value = obj, grad = grad, eta = eta, lambda = L,
         alpha = p$alpha, beta = p$beta, kappa = p$kappa)
  }
  fn <- function(th) {
    ans <- calc(th)
    if (is.null(ans)) return(1e100)
    ans$value
  }
  gr <- function(th) {
    ans <- calc(th)
    if (is.null(ans)) return(rep(0, length(th)))
    ans$grad
  }
  initial <- fn(theta0)
  opt <- optim(theta0, fn = fn, gr = gr, method = "BFGS",
               control = list(maxit = 4000L, reltol = 1e-10))
  fitted <- calc(opt$par)
  if (is.null(fitted)) stop("优化产生了不可用的参数。")
  if (!is.finite(fitted$value) || fitted$value > initial + 1e-6) {
    stop("拟合目标函数没有相对于 SVD 初始化改善。")
  }
  if (max(abs(fitted$grad)) > 0.02) {
    warning("梯度仍然偏大：需要进一步检查该死因的收敛。")
  }
  stopifnot(abs(sum(fitted$beta) - 1) < 1e-8,
            abs(sum(fitted$kappa)) < 1e-8)
  list(alpha = fitted$alpha, beta = fitted$beta,
       kappa = fitted$kappa, rate_fitted = exp(fitted$eta),
       count_fitted = fitted$lambda,
       deviance = 2 * total_d * fitted$value,
       objective_initial = initial,
       objective_final = fitted$value,
       convergence = opt$convergence,
       max_abs_gradient = max(abs(fitted$grad)))
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
if (any(diagnostics$convergence_code != 0L)) {
  print(diagnostics[diagnostics$convergence_code != 0L, ])
  stop("至少一个模型未达到优化器收敛判据；不保存正式结果。")
}

# ---- 3. Preview parameter table and Figure 3 in RStudio ----
cat("\n========================================================\n")
cat("PREVIEW ONLY — FIGURE 3. NO OUTPUT FILES SAVED.\n")
cat("========================================================\n")
cat("\n模型拟合汇总（检查 convergence_code 都为 0）：\n")
print(diagnostics, row.names = FALSE)
cat("\n示例：男性心血管疾病（如包含）的一维时间因子前 10 年：\n")
print(head(kappa_table[kappa_table$sex_name == "Male" &
                          kappa_table$cause_name == chosen[1], ], 10),
      row.names = FALSE)
cat("\n识别约束核验：每个死因/性别 sum(beta)=1，mean(kappa)=0。\n")
for (sx in sexes) {
  for (cause in chosen) {
    m <- results[[paste(sx, cause, sep = "__")]]
    if (abs(sum(m$beta) - 1) > 1e-8 ||
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
  on.exit(par(oldpar))
  par(mfrow = c(2, 3), mar = c(3.5, 3.9, 3.0, 0.8),
      oma = c(6.7, 0.3, 2.6, 0.3))
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
  par(xpd = NA)
  legend("bottom", inset = c(0, -0.13), legend = chosen,
         col = cols, lty = 1, lwd = 1.8, cex = 0.61,
         ncol = 3, bty = "n")
  mtext("GBD Hungary | Poisson Lee-Carter (Level 2 causes)",
        side = 3, outer = TRUE, line = 1.1, cex = 1.15, font = 2)
  mtext("Fit: 1980-2015 | Ages 60-89 | Male/Female separate | GBD quasi-likelihood",
        side = 1, outer = TRUE, line = 1.0, cex = 0.75)
}
cat("\n请查看 RStudio Plots 窗格的 Figure 3 六面板图：alpha、beta、kappa × 男女。\n")
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
                   "Identification sum(beta)=1, sum(kappa)=0.",
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
