# Figure 2 | Varga (2025), s3.1: age-standardized cause-specific mortality
# GBD 2023 Hungary METHOD replication; not the original HCSO numerical result.
#
# Run in RStudio from the GitHub REPOSITORY ROOT:
#   source("02_figure2/run.R")
#
# Input CSV files (keep the original exported column names):
#   data/raw/IHME-GBD_2023_DATA-pop.csv
#   data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv  [Deaths Number + Rate]
# If the deaths file was renamed on your PC, change deaths_path below.
# Do not commit raw GBD files to GitHub.

population_path <- "data/raw/IHME-GBD_2023_DATA-pop.csv"
deaths_path     <- "data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv"
outdir          <- "outputs/02_figure2"

years <- 1980:2023
standard_year <- 1980L
sexes <- c("Male", "Female")
ages <- c("60-64 years", "65-69 years", "70-74 years",
          "75-79 years", "80-84 years", "85-89 years")
n_population <- length(years) * length(sexes) * length(ages)
n_causes <- 19L
n_rate_rows <- n_population * n_causes

# We use a common, BOTH-SEX 1980 age distribution as a transparent
# comparability convention. The original used 1970 standard population,
# included ages 90+, and may have used different sex weights.
# Thus this is NOT a direct replication of the original ASMR values.

# ---- 1. Independently downloaded GBD Population, all 44 years ----
if (!file.exists(population_path)) {
  stop("Missing Population CSV: ", population_path,
       "\nPut the user-uploaded original file in data/raw/.")
}
pop <- read.csv(population_path, check.names = FALSE,
                stringsAsFactors = FALSE)
pop_fields <- c("year", "location_name", "sex_name", "age_name",
                "measure_name", "metric_name", "val", "lower", "upper")
stopifnot(all(pop_fields %in% names(pop)))

if (anyNA(pop[, pop_fields]) ||
    !setequal(unique(pop$year), years) ||
    !setequal(unique(pop$sex_name), sexes) ||
    !setequal(unique(pop$age_name), ages) ||
    !all(pop$location_name == "Hungary") ||
    !all(pop$measure_name == "Population") ||
    !all(pop$metric_name == "Number")) {
  stop("Population scope/field check failed: verify country, years, sexes, ages and metric.")
}
if ("population_group_name" %in% names(pop) &&
    !all(pop$population_group_name == "All Population")) {
  stop("Population group must be 'All Population'.")
}
if (nrow(pop) != n_population ||
    anyDuplicated(pop[c("year", "sex_name", "age_name")]) > 0L ||
    !all(is.finite(pop$val)) || !all(pop$val > 0) ||
    !all(is.finite(pop$lower)) || !all(is.finite(pop$upper)) ||
    !all(pop$lower <= pop$val & pop$val <= pop$upper)) {
  stop("Population records have missing/duplicate/invalid population or intervals.")
}
pop <- pop[c("year", "sex_name", "age_name", "val", "lower", "upper")]
names(pop)[names(pop) == "val"] <- "population"
pop <- pop[order(pop$year, match(pop$age_name, ages),
                 match(pop$sex_name, sexes)), ]

# ---- 2. Standardization weights from BOTH sexes in 1980 ----
base <- pop[pop$year == standard_year, ]
stopifnot(nrow(base) == length(sexes) * length(ages))
standard <- aggregate(population ~ age_name, data = base, FUN = sum)
standard <- standard[match(ages, standard$age_name), ]
stopifnot(nrow(standard) == length(ages), !anyNA(standard$population))
names(standard)[names(standard) == "population"] <- "population_1980_both_sexes"
standard$weight <- standard$population_1980_both_sexes /
                   sum(standard$population_1980_both_sexes)
stopifnot(all(is.finite(standard$weight)), all(standard$weight > 0),
          abs(sum(standard$weight) - 1) < 1e-12)

# Known data check (rounded total only): independently inspected uploaded CSV.
# This is a reproducibility assertion, not a replacement for actual data.
stopifnot(abs(sum(standard$population_1980_both_sexes) -
              1789572.6757876906) < 1)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(pop, file.path(outdir, "population_1980_2023_clean.csv"),
          row.names = FALSE)
write.csv(standard, file.path(outdir, "standard_age_weights_1980.csv"),
          row.names = FALSE)

cat("Population file PASSED: ", nrow(pop),
    " rows; 1980 pooled ages 60-89 population = ",
    format(sum(standard$population_1980_both_sexes), digits = 12),
    "; weight sum = ", sum(standard$weight), "\n", sep = "")

# ---- 3. Check Deaths Number and age-specific Rate data ----
if (!file.exists(deaths_path)) {
  stop("Population checks passed and weights were saved, but Deaths CSV is absent: ",
       deaths_path,
       "\nChange deaths_path to your actual local Deaths filename.")
}
gbd <- read.csv(deaths_path, check.names = FALSE,
                stringsAsFactors = FALSE)
death_fields <- c("year", "location_name", "measure_name",
                  "sex_name", "age_name", "cause_name",
                  "metric_name", "val")
stopifnot(all(death_fields %in% names(gbd)))
if (anyNA(gbd[, death_fields]) ||
    !all(gbd$location_name == "Hungary") ||
    !all(gbd$measure_name == "Deaths") ||
    !setequal(unique(gbd$year), years) ||
    !setequal(unique(gbd$sex_name), sexes) ||
    !setequal(unique(gbd$age_name), ages) ||
    !setequal(unique(gbd$metric_name), c("Number", "Rate")) ||
    length(unique(gbd$cause_name)) != n_causes ||
    nrow(gbd) != 2L * n_rate_rows ||
    anyDuplicated(gbd[c("year", "sex_name", "age_name",
                         "cause_name", "metric_name")]) > 0L ||
    !all(is.finite(gbd$val)) || !all(gbd$val >= 0)) {
  stop("Deaths CSV is incomplete, duplicated or from the wrong GBD scope.")
}
if (all(c("lower", "upper") %in% names(gbd))) {
  if (anyNA(gbd[c("lower", "upper")]) ||
      !all(is.finite(gbd$lower)) || !all(is.finite(gbd$upper)) ||
      !all(gbd$lower <= gbd$val & gbd$val <= gbd$upper)) {
    stop("Deaths CSV uncertainty interval validation failed.")
  }
}

keys <- c("year", "sex_name", "age_name", "cause_name")
rate <- gbd[gbd$metric_name == "Rate", c(keys, "val")]
num  <- gbd[gbd$metric_name == "Number", c(keys, "val")]
names(rate)[names(rate) == "val"] <- "rate_per_100k"
names(num)[names(num) == "val"] <- "death_number"
paired <- merge(rate, num, by = keys, all = FALSE, sort = FALSE)
stopifnot(nrow(rate) == n_rate_rows, nrow(num) == n_rate_rows,
          nrow(paired) == n_rate_rows)

# ---- 4. Cross-check Deaths Number/Rate against INDEPENDENT Population ----
# This uses the WHO/IHME convention for age-specific GBD mortality:
# rate is deaths per 100,000 of the corresponding population.
# GBD posterior means need not obey the identity EXACTLY because
# E[D/P] need not equal E[D]/E[P]; deviations are audited, not hidden.
qc <- merge(paired, pop[c("year", "sex_name", "age_name", "population")],
            by = c("year", "sex_name", "age_name"),
            all.x = TRUE, sort = FALSE)
stopifnot(nrow(qc) == n_rate_rows,
          !anyNA(qc$population),
          all(qc$population > 0))
qc$rate_recomputed_per_100k <- 100000 * qc$death_number / qc$population
qc$relative_error <- ifelse(
  qc$rate_per_100k > 0,
  abs(qc$rate_recomputed_per_100k - qc$rate_per_100k) /
    qc$rate_per_100k,
  ifelse(qc$rate_recomputed_per_100k == 0, 0, Inf)
)
qc <- qc[order(qc$year, qc$sex_name, qc$age_name,
               qc$cause_name), ]
write.csv(qc, file.path(outdir, "number_rate_population_check.csv"),
          row.names = FALSE)

if (!all(is.finite(qc$relative_error))) {
  stop("Inconsistent zero Rate/Number observations; inspect QC file.")
}
med_error <- median(qc$relative_error)
p95_error <- as.numeric(quantile(qc$relative_error, 0.95, names = FALSE))
max_error <- max(qc$relative_error)
cat("Rate consistency: median relative difference=", med_error,
    "; p95=", p95_error, "; max=", max_error, "\n")
# QC alert thresholds (not an inferential test of GBD uncertainty).
# Large systematic discrepancies could mean mismatched population rounds,
# non-matching mortality metrics, or incorrectly labelled population files.
if (med_error > 0.02 || p95_error > 0.10) {
  stop("Rates and independent population fail QC screen (median >2% or p95 >10%).",
       "\nOpen outputs/02_figure2/number_rate_population_check.csv.")
}
if (max_error > 0.10) {
  warning("Some individual cells have >10% Rate vs Number/Population difference.",
          " Review number_rate_population_check.csv.")
}

# ---- 5. Compute ASMR for each sex x year x cause ----
# Never divide the age-specific GBD Rate by population AGAIN.
# Convert GBD per-100,000 Rate to per-person-year by dividing by 100000.
rates_w <- merge(rate, standard[c("age_name", "weight")],
                 by = "age_name", all.x = TRUE, sort = FALSE)
stopifnot(nrow(rates_w) == n_rate_rows, !anyNA(rates_w$weight))
rates_w$per_person_year <- rates_w$rate_per_100k / 100000
rates_w$weighted <- rates_w$per_person_year * rates_w$weight
asmr <- aggregate(weighted ~ year + sex_name + cause_name,
                  data = rates_w, FUN = sum)
names(asmr)[names(asmr) == "weighted"] <- "asmr_per_person_year"
asmr$asmr_per_100k <- 100000 * asmr$asmr_per_person_year
asmr$standard_year <- standard_year
asmr <- asmr[order(asmr$cause_name, asmr$year, asmr$sex_name), ]
stopifnot(nrow(asmr) == length(years) * length(sexes) * n_causes,
          all(is.finite(asmr$asmr_per_person_year)),
          all(asmr$asmr_per_person_year >= 0))

# A weighted arithmetic mean must lie between the age-specific min/max.
lo <- aggregate(per_person_year ~ year + sex_name + cause_name,
                data = rates_w, FUN = min)
hi <- aggregate(per_person_year ~ year + sex_name + cause_name,
                data = rates_w, FUN = max)
names(lo)[names(lo) == "per_person_year"] <- "min_age_rate"
names(hi)[names(hi) == "per_person_year"] <- "max_age_rate"
bounds <- merge(lo, hi, by = c("year", "sex_name", "cause_name"))
check <- merge(asmr, bounds, by = c("year", "sex_name", "cause_name"))
stopifnot(nrow(check) == nrow(asmr),
          all(check$asmr_per_person_year >= check$min_age_rate - 1e-12),
          all(check$asmr_per_person_year <= check$max_age_rate + 1e-12))
write.csv(asmr, file.path(outdir, "asmr_19_causes_1980_2023.csv"),
          row.names = FALSE)

# ---- 6. Choose the nine main GBD Level-2 causes, using TRAINING years ----
# These are NOT the same nine ICD chapters in the original Varga Figure 2.
rankdata <- num[num$year <= 2015, c("cause_name", "death_number")]
ranking <- aggregate(death_number ~ cause_name, data = rankdata, FUN = sum)
ranking <- ranking[order(-ranking$death_number, ranking$cause_name), ]
top9 <- head(ranking$cause_name, 9)
write.csv(data.frame(rank = seq_along(top9), cause_name = top9),
          file.path(outdir, "figure2_selected_9_causes.csv"),
          row.names = FALSE)

# ---- 7. Figure 2-style plot: 3 x 3 panels, Male solid, Female dashed ----
png(file.path(outdir, "figure2_gbd_asmr_9_panels.png"),
    width = 2250, height = 1950, res = 180)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3.7, 4.7, 4.7, 1.4),
    oma = c(2.7, 1, 3.2, 1))
cols <- grDevices::hcl.colors(9L, palette = "Dark 3")
for (i in seq_along(top9)) {
  z <- asmr[asmr$cause_name == top9[i], ]
  male <- z[z$sex_name == "Male", ]
  female <- z[z$sex_name == "Female", ]
  male <- male[order(male$year), ]
  female <- female[order(female$year), ]
  stopifnot(nrow(male) == length(years),
            nrow(female) == length(years),
            identical(male$year, female$year))
  ylim <- range(c(male$asmr_per_person_year,
                  female$asmr_per_person_year))
  if (diff(ylim) == 0) ylim <- ylim + c(-1, 1) * 1e-10
  ttl <- paste(strwrap(top9[i], width = 26), collapse = "\n")
  plot(male$year, male$asmr_per_person_year, type = "l",
       col = cols[i], lwd = 1.7, lty = 1,
       xlab = "Year", ylab = "ASMR", main = ttl,
       ylim = ylim, xaxt = "n", cex.main = 0.82)
  axis(1, at = c(1980, 2000, 2020))
  lines(female$year, female$asmr_per_person_year,
        col = cols[i], lwd = 1.7, lty = 2)
  legend("topright", legend = c("Males", "Females"),
         col = cols[i], lty = c(1, 2), lwd = 1.7, cex = 0.65, bty = "n")
}
mtext("GBD Hungary, ages 60-89, 1980-2023", outer = TRUE,
      side = 3, line = 1.2, font = 2, cex = 1.05)
mtext("1980 pooled age weights (independent GBD Population); not HCSO Figure 2",
      outer = TRUE, side = 1, line = 1.0, cex = 0.72)
par(oldpar)
dev.off()

cat("Figure 2 calculation finished.\n")
cat("1980 pooled standard weights: ",
    paste(round(standard$weight, 6), collapse = ", "), "\n", sep = "")
cat("ASMR series: ", nrow(asmr), " rows; nine panels: ",
    paste(top9, collapse = "; "), "\n", sep = "")
cat("Outputs in ", outdir, "\n", sep = "")
