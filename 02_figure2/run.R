# Figure 2: Age-standardized cause-specific mortality rates (ASMR)
# Reviewed 2026-11-09: require an INDEPENDENT standard-population file.
# Method replication, not HCSO 1970-2021 numerical replication.
#
# Execute from the project root in RStudio:
# source("02_figure2/run.R")
#
# Inputs:
# 1) Existing GBD 2023 Hungary cause-specific deaths CSV, Rate AND Number.
# 2) NEW: data/raw/gbd_population_1980.csv from GBD 2023 Population.
#    Population is required ONLY for baseline-year fixed age weights.
#
# GBD sample: Hungary; 1980-2023; Male/Female; 60-89 in six 5-year bins.
# ASMR = sum_{age} (age-specific rate per 100,000 / 100,000) * weight_1980.
# The age-specific GBD Rate is NOT a precomputed age-standardized rate.

deaths_path <- "data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv"
population_path <- "data/raw/gbd_population_1980.csv"
ages <- c("60-64 years", "65-69 years", "70-74 years",
          "75-79 years", "80-84 years", "85-89 years")
sexes <- c("Male", "Female")
years <- 1980:2023
outdir <- "outputs/02_figure2"

if (!file.exists(deaths_path)) {
  stop("Deaths CSV not found: ", deaths_path,
       "\nUpdate deaths_path if your downloaded filename differs.")
}
if (!file.exists(population_path)) {
  stop("Population CSV not found: ", population_path,
       "\nDownload GBD 2023 Population, Hungary, year 1980,",
       " Male/Female, ages 60-64...85-89, metric Number.",
       "\nSave as data/raw/gbd_population_1980.csv.",
       "\nThe previous Number/Rate population proxy is NOT used in this revised script.")
}

# --- Part A: validate deaths and age-specific rates ---
gbd <- read.csv(deaths_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("year", "sex_name", "age_name", "cause_name",
              "metric_name", "val", "measure_name", "location_name")
stopifnot(all(required %in% names(gbd)))
stopifnot(all(gbd$measure_name == "Deaths"),
          all(gbd$location_name == "Hungary"),
          setequal(unique(gbd$year), years),
          setequal(unique(gbd$sex_name), sexes),
          setequal(unique(gbd$age_name), ages),
          setequal(unique(gbd$metric_name), c("Number", "Rate")),
          length(unique(gbd$cause_name)) == 19L,
          nrow(gbd) == 44L * 2L * 6L * 19L * 2L,
          all(is.finite(gbd$val)), all(gbd$val >= 0))

keys <- c("year", "sex_name", "age_name", "cause_name")
stopifnot(!anyDuplicated(gbd[, c(keys, "metric_name")]))
rates <- gbd[gbd$metric_name == "Rate", c(keys, "val")]
names(rates)[names(rates) == "val"] <- "rate_per_100k"
stopifnot(nrow(rates) == 44L * 2L * 6L * 19L,
          !anyNA(rates[, c(keys, "rate_per_100k")]))
# Population standardization requires AGE-SPECIFIC Rate.
# For additional certainty, check the age labels above: no "Age-standardized".
if (any(grepl("standard", rates$age_name, ignore.case = TRUE))) {
  stop("Rate contains an age-standardized age label; do not standardize twice.")
}

# --- Part B: read independent 1980 population ---
# GBD 2023 Population data may include more years or columns.
# The file must contain year, sex_name, age_name, and val.
pop <- read.csv(population_path, stringsAsFactors = FALSE,
                check.names = FALSE)
pop_required <- c("year", "sex_name", "age_name", "val")
if (!all(pop_required %in% names(pop))) {
  stop("Population CSV must include columns: ",
       paste(pop_required, collapse = ", "))
}
if ("location_name" %in% names(pop)) {
  if (!all(pop$location_name == "Hungary")) {
    stop("Population file contains locations other than Hungary.")
  }
}
if ("measure_name" %in% names(pop)) {
  if (!all(grepl("population", pop$measure_name, ignore.case = TRUE))) {
    stop("Population file measure_name is not Population.")
  }
}
# GBD population must be a COUNT, not a Rate or Percent.
if ("metric_name" %in% names(pop)) {
  if (!all(pop$metric_name == "Number")) {
    stop("Population file must use metric_name = Number.")
  }
}
pop <- pop[pop$year == 1980 & pop$sex_name %in% sexes &
             pop$age_name %in% ages, c(pop_required)]
stopifnot(nrow(pop) == 12L,
          !anyDuplicated(pop[, c("year", "sex_name", "age_name")]),
          setequal(unique(pop$sex_name), sexes),
          setequal(unique(pop$age_name), ages),
          all(is.finite(pop$val)), all(pop$val > 0))

# --- Part C: one pooled standard age structure for both sexes ---
# This is OUR documented convention; the paper only says "1970 mid-year
# population" and does not specify whether sex-specific standards were used.
# Original paper also includes the older 90+ population group.
standard <- aggregate(val ~ age_name, data = pop, FUN = sum)
standard <- standard[match(ages, standard$age_name), ]
names(standard)[names(standard) == "val"] <- "population_1980_both_sexes"
standard$weight <- standard$population_1980_both_sexes /
                   sum(standard$population_1980_both_sexes)
stopifnot(nrow(standard) == 6L, !anyNA(standard$weight),
          abs(sum(standard$weight) - 1) < 1e-12)

# --- Part D: calculate age-standardized rates ---
rates <- merge(rates, standard[, c("age_name", "weight")],
               by = "age_name", all.x = TRUE, sort = FALSE)
stopifnot(nrow(rates) == 44L * 2L * 6L * 19L,
          all(is.finite(rates$weight)))
rates$rate_per_person_year <- rates$rate_per_100k / 100000
rates$weighted_rate <- rates$rate_per_person_year * rates$weight
asmr <- aggregate(weighted_rate ~ year + sex_name + cause_name,
                  data = rates, FUN = sum)
names(asmr)[names(asmr) == "weighted_rate"] <- "asmr_per_person_year"
asmr$asmr_per_100k <- asmr$asmr_per_person_year * 100000
asmr$standard_year <- 1980L
asmr <- asmr[order(asmr$cause_name, asmr$year, asmr$sex_name), ]
stopifnot(nrow(asmr) == 44L * 2L * 19L,
          all(is.finite(asmr$asmr_per_person_year)),
          all(asmr$asmr_per_person_year >= 0))

# A weighted mean must fall between the minimum and maximum age-specific rates.
extremes <- aggregate(rate_per_person_year ~ year + sex_name + cause_name,
                      data = rates, FUN = function(z) c(min(z), max(z)))
# aggregate() with a length-two vector FUN returns a matrix column.
rate_ranges <- extremes$rate_per_person_year
if (!is.matrix(rate_ranges) || ncol(rate_ranges) != 2L) {
  stop("Unexpected aggregate() bounds structure.")
}
mins <- rate_ranges[, 1]
maxs <- rate_ranges[, 2]
bounds <- data.frame(extremes[, c("year", "sex_name", "cause_name")],
                     min_rate = mins, max_rate = maxs)
check <- merge(asmr, bounds, by = c("year", "sex_name", "cause_name"))
stopifnot(nrow(check) == nrow(asmr),
          all(check$asmr_per_person_year >= check$min_rate - 1e-12),
          all(check$asmr_per_person_year <= check$max_rate + 1e-12))

# --- Part E: nine GBD categories selected on 1980-2015 only ---
deaths <- gbd[gbd$metric_name == "Number" & gbd$year <= 2015,
              c("year", "cause_name", "val")]
ranking <- aggregate(val ~ cause_name, data = deaths, FUN = sum)
ranking <- ranking[order(-ranking$val, ranking$cause_name), ]
top9 <- head(ranking$cause_name, 9L)

# --- Part F: save verified tables and original-STYLE 3 x 3 figure ---
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(standard, file.path(outdir, "standard_age_weights_1980.csv"),
          row.names = FALSE)
write.csv(asmr, file.path(outdir, "asmr_19_causes_1980_2023.csv"),
          row.names = FALSE)
write.csv(data.frame(rank = seq_along(top9), cause_name = top9),
          file.path(outdir, "figure2_selected_9_causes.csv"),
          row.names = FALSE)

png(file.path(outdir, "figure2_gbd_asmr_9_panels.png"),
    width = 2250, height = 1950, res = 180)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3.7, 4.7, 4.7, 1.4),
    oma = c(2.7, 1, 3.2, 1))
cols <- grDevices::hcl.colors(9, palette = "Dark 3")
for (i in seq_along(top9)) {
  z <- asmr[asmr$cause_name == top9[i], ]
  male <- z[z$sex_name == "Male", ]
  female <- z[z$sex_name == "Female", ]
  male <- male[order(male$year), ]
  female <- female[order(female$year), ]
  stopifnot(nrow(male) == 44L, nrow(female) == 44L,
            identical(male$year, female$year))
  yspan <- range(c(male$asmr_per_person_year,
                   female$asmr_per_person_year))
  if (diff(yspan) == 0) yspan <- yspan + c(-1, 1) * 1e-10
  ttl <- paste(strwrap(top9[i], width = 26), collapse = "\n")
  plot(male$year, male$asmr_per_person_year, type = "l",
       col = cols[i], lwd = 1.7, lty = 1,
       xlab = "Year", ylab = "ASMR", main = ttl,
       ylim = yspan, xaxt = "n", cex.main = 0.82)
  axis(1, at = c(1980, 2000, 2020))
  lines(female$year, female$asmr_per_person_year,
        col = cols[i], lwd = 1.7, lty = 2)
  legend("topright", legend = c("Males", "Females"),
         col = cols[i], lty = c(1, 2), lwd = 1.7,
         cex = 0.65, bty = "n")
}
mtext("GBD Hungary, 60-89 years (1980-2023)", outer = TRUE,
      side = 3, line = 1.2, font = 2, cex = 1.05)
mtext("1980 pooled population weights; NOT original HCSO Figure 2",
      outer = TRUE, side = 1, line = 1.0, cex = 0.75)
par(oldpar)
dev.off()

cat("Figure 2 method-replication calculation finished.\n")
cat("Source for age weights: independent 1980 GBD Population CSV.\n")
cat("Sum of six weights:", sum(standard$weight), "\n")
cat("Number of ASMR rows:", nrow(asmr), "\n")
cat("Nine GBD categories:", paste(top9, collapse = "; "), "\n")
cat("Files:", outdir, "\n")
