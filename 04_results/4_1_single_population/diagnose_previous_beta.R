# §4.1 | Audit the PREVIOUS beta.csv / kappa.csv, no Deaths CSV required.
# This script does NOT refit any model or silently overwrite any files.
# In RStudio: Source this script. Select old beta.csv, then old kappa.csv.
cat("请选择旧版 Figure 3 的 beta.csv：\n")
b <- read.csv(file.choose(), stringsAsFactors = FALSE)
cat("请选择旧版 Figure 3 的 kappa.csv：\n")
k <- read.csv(file.choose(), stringsAsFactors = FALSE)
need_b <- c("sex_name", "cause_name", "age_name", "beta")
need_k <- c("sex_name", "cause_name", "year", "kappa")
stopifnot(all(need_b %in% names(b)), all(need_k %in% names(k)))
key_b <- c("sex_name", "cause_name", "age_name")
key_k <- c("sex_name", "cause_name", "year")
stopifnot(!anyNA(b[need_b]), !anyNA(k[need_k]),
          !anyDuplicated(b[key_b]), !anyDuplicated(k[key_k]),
          all(is.finite(b$beta)), all(is.finite(k$kappa)))
groups <- unique(b[c("sex_name", "cause_name")])
report <- data.frame()
normalized_beta <- data.frame()
normalized_kappa <- data.frame()
for (i in seq_len(nrow(groups))) {
  sx <- groups$sex_name[i]
  cs <- groups$cause_name[i]
  bb <- b[b$sex_name == sx & b$cause_name == cs, ]
  kk <- k[k$sex_name == sx & k$cause_name == cs, ]
  if (nrow(bb) != 6L || nrow(kk) != 36L ||
      !identical(sort(kk$year), 1980:2015)) stop("旧参数表不完整：", cs, sx)
  scale <- max(abs(bb$beta))
  if (scale == 0) stop("beta 全零，无法识别时期因子：", cs, sx)
  bn <- bb$beta / scale
  kn <- kk$kappa * scale
  if (sum(bn) < 0) {
    bn <- -bn
    kn <- -kn
  }
  # Numerical identity: beta*kappa exactly preserved at each age/year.
  original_surface <- outer(bb$beta, kk$kappa)
  rescaled_surface <- outer(bn, kn)
  max_error <- max(abs(original_surface - rescaled_surface))
  stopifnot(max_error < 1e-8)
  cancellation <- abs(sum(bb$beta)) / sum(abs(bb$beta))
  report <- rbind(report, data.frame(
    sex_name = sx, cause_name = cs,
    sum_beta_old = sum(bb$beta), max_abs_beta_old = scale,
    cancellation_ratio = cancellation,
    sd_kappa_old = sd(kk$kappa),
    sd_kappa_normalized = sd(kn),
    identical_fit_max_abs_error = max_error,
    flag = if (cancellation < 0.08) "REVIEW" else "OK"
  ))
  bb$beta <- bn
  kk$kappa <- kn
  normalized_beta <- rbind(normalized_beta, bb)
  normalized_kappa <- rbind(normalized_kappa, kk)
}
report <- report[order(report$flag, -report$max_abs_beta_old), ]
cat("\n旧模型参数异常检查：\n")
print(report, row.names = FALSE, digits = 5)
cat("\n结论：仅仅归一化不会改变模型拟合值，也不能证明模型稳定。\n")
cat("若出现 REVIEW，需检查 Deaths/Population + 重新拟合，并对 DTW 进行敏感性分析。\n")
if (interactive()) {
  View(report, "§4.1 旧版 beta/kappa 缩放异常检查")
}
cat("\n本程序默认不保存任何文件。可在 RStudio 环境中查看 report、normalized_beta、normalized_kappa。\n")
