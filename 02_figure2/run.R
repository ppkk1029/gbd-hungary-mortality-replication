# 02 | Figure 2: age-standardized cause-specific mortality rates (ASMR)
# Varga (2025), Equation (5), method replication using GBD 2023
#
# Run from the repository ROOT in RStudio:
# source("02_figure2/run.R")
#
# Input: data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv
# No external R packages required.
# Original paper: HCSO, 1970-2021, age 60+, 1970 standard population.
# Here: GBD 2023, 1980-2023, ages 60-89, 1980 standard population.
#
# IMPORTANT: outputs are not original HCSO data or original-paper values.

input <- "data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv"
if (!file.exists(input)) {
  stop("Missing raw GBD CSV: ", input,
       "\nDownload the file and place it in data/raw/; do not rename columns.")
}

gbd <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("year", "location_name", "measure_name", "sex_name",
              "age_name", "cause_name", "metric_name", "val")
stopifnot(all(required %in% names(gbd)))

ages <- c("60-64 years", "65-69 years", "70-74 years",
          "75-79 years", "80-84 years", "85-89 years")
sexes <- c("Male", "Female")
years <- 1980:2023
keys <- c("year", "sex_name", "age_name", "cause_name")

# Part 1. Validate the EXACT dataset (no implicit age/year/cause mixing).
stopifnot(
  all(gbd$location_name == "Hungary"),
  all(gbd$measure_name == "Deaths"),
  setequal(unique(gbd$year), years),
  setequal(unique(gbd$sex_name), sexes),
  setequal(unique(gbd$age_name), ages),
  setequal(unique(gbd$metric_name), c("Number", "Rate")),
  length(unique(gbd$cause_name)) == 19L,
  nrow(gbd) == 44L * 2L * 6L * 19L * 2L,
  !anyNA(gbd[, c(keys, "metric_name", "val")]),
  !anyDuplicated(gbd[, c(keys, "metric_name")]),
  all(is.finite(gbd$val)), all(gbd$val >= 0)
)
dir.create("outputs/02_figure2", recursive = TRUE, showWarnings = FALSE)
out <- "outputs/02_figure2"

# Part 2. Pair Number and age-specific Rate, matching ALL four keys.
deaths <- subset(gbd, metric_name == "Number", select = c(keys, "val"))
rates  <- subset(gbd, metric_name == "Rate",   select = c(keys, "val"))
names(deaths)[ncol(deaths)] <- "death_number"
names(rates)[ncol(rates)] <- "rate_per_100k"

d <- merge(deaths, rates, by = keys, all = FALSE, sort = FALSE)
stopifnot(nrow(d) == 44L * 2L * 6L * 19L)
stopifnot(all(d$rate_per_100k >= 0))
if (any(xor(d$death_number == 0, d$rate_per_100k == 0))) {
  stop("Number and Rate disagree on zero cells; cannot infer exposures.")
}

# In GBD age-specific death rates, Rate is reported per 100,000 persons.
# If N = deaths, R = deaths / population * 100000, then
# population = N * 100000 / R. Undefined if N = R = 0.
# This is a DERIVED exposure proxy, NOT observed HCSO central exposure.
d$population_implied <- ifelse(
  d$rate_per_100k > 0,
  d$death_number * 100000 / d$rate_per_100k,
  NA_real_
)

# The same sex-age-year population should be implied by every cause.
# Median over available causes protects against rare 0/0 cells.
pop <- aggregate(
  population_implied ~ year + sex_name + age_name,
  data = d, FUN = function(z) median(z[is.finite(z) & z > 0]),
  na.action = na.pass
)
names(pop)[names(pop) == "population_implied"] <- "population_estimated"
stopifnot(nrow(pop) == 44L * 2L * 6L,
          all(is.finite(pop$population_estimated)),
          all(pop$population_estimated > 0))

check <- merge(d, pop, by = c("year", "sex_name", "age_name"), sort = FALSE)
valid <- is.finite(check$population_implied) & check$population_implied > 0
rel_err <- abs(check$population_implied[valid] -
                 check$population_estimated[valid]) /
           check$population_estimated[valid]
max_rel_error <- max(rel_err)
if (!is.finite(max_rel_error) || max_rel_error > 0.01) {
  stop(sprintf(
    "Number/Rate imply inconsistent denominators across causes: max rel error %.6f. Stop and check GBD metric units/hierarchy.",
    max_rel_error
  ))
}
write.csv(pop, file.path(out, "implied_population_by_age_sex_year.csv"),
          row.names = FALSE)

# Part 3. 1980 common standard age distribution.
# Pool Male + Female ONLY to construct ONE standard age weight for both.
# Using sex-specific weights would undermine direct sex comparison.
base <- subset(pop, year == 1980)
standard <- aggregate(population_estimated ~ age_name, data = base, FUN = sum)
standard <- standard[match(ages, standard$age_name), ]
stopifnot(nrow(standard) == 6L, !anyNA(standard$population_estimated))
standard$weight <- standard$population_estimated /
                   sum(standard$population_estimated)
stopifnot(all(standard$weight > 0),
          abs(sum(standard$weight) - 1) < 1e-12)
names(standard)[names(standard) == "population_estimated"] <-
  "population_1980_both_sexes"
write.csv(standard, file.path(out, "standard_age_weights_1980.csv"),
          row.names = FALSE)

# Part 4. Age-standardize each sex-cause-year series.
# ASMR[t, sex, cause] = sum_age (Rate[t,sex,age,cause] / 100000) * w[age].
d <- merge(d, standard[, c("age_name", "weight")],
           by = "age_name", all.x = TRUE, sort = FALSE)
stopifnot(nrow(d) == 44L * 2L * 6L * 19L, !anyNA(d$weight))
d$rate_per_person_year <- d$rate_per_100k / 100000
d$weighted_contribution <- d$rate_per_person_year * d$weight
asmr <- aggregate(
  weighted_contribution ~ year + sex_name + cause_name,
  data = d, FUN = sum
)
names(asmr)[names(asmr) == "weighted_contribution"] <-
  "asmr_per_person_year"
asmr$asmr_per_100k <- asmr$asmr_per_person_year * 100000
asmr$standard_year <- 1980L
asmr <- asmr[order(asmr$cause_name, asmr$year, asmr$sex_name), ]
stopifnot(nrow(asmr) == 44L * 2L * 19L,
          all(is.finite(asmr$asmr_per_person_year)),
          all(asmr$asmr_per_person_year >= 0))

# Sanity check: weighted average must lie between min/max age-specific rates
extremes <- aggregate(rate_per_person_year ~ year + sex_name + cause_name,
  d, FUN = function(z) c(min = min(z), max = max(z)))
extremes <- data.frame(
  extremes[c("year", "sex_name", "cause_name")],
  min_rate = extremes$rate_per_person_year[, "min"],
  max_rate = extremes$rate_per_person_year[, "max"]
)
chk <- merge(asmr, extremes, by = c("year", "sex_name", "cause_name"))
stopifnot(all(chk$asmr_per_person_year >= chk$min_rate - 1e-12),
          all(chk$asmr_per_person_year <= chk$max_rate + 1e-12))

write.csv(asmr, file.path(out, "asmr_19_causes_1980_2023.csv"),
          row.names = FALSE)

# Part 5. Match Figure 1's PRE-TEST selection: the nine highest
# cumulative Number causes during 1980-2015 (pooled sexes/ages).
ranking <- aggregate(death_number ~ cause_name,
  data = subset(d, year <= 2015), FUN = sum)
ranking <- ranking[order(-ranking$death_number, ranking$cause_name), ]
top9 <- head(ranking$cause_name, 9L)
write.csv(data.frame(rank = 1:9, cause_name = top9),
  file.path(out, "figure2_selected_9_causes.csv"), row.names = FALSE)

# Part 6. Reproduce original Figure 2's 3-by-3 layout and line styles.
# One panel/cause: Male = solid, Female = dashed.
png(file.path(out, "figure2_gbd_asmr_9_panels.png"),
    width = 2250, height = 1950, res = 180)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(3.7, 4.7, 4.7, 1.4),
    oma = c(2.7, 1, 3.2, 1))
cols <- grDevices::hcl.colors(9, palette = "Dark 3")
for (i in seq_along(top9)) {
  cause <- top9[i]
  z <- subset(asmr, cause_name == cause)
  male <- subset(z, sex_name == "Male")
  female <- subset(z, sex_name == "Female")
  male <- male[order(male$year), ]
  female <- female[order(female$year), ]
  stopifnot(nrow(male) == 44L, nrow(female) == 44L,
            identical(male$year, female$year))
  yrange <- range(c(male$asmr_per_person_year,
                    female$asmr_per_person_year))
  if (diff(yrange) == 0) yrange <- yrange + c(-1, 1) * 1e-10
  title_text <- paste(strwrap(cause, width = 27), collapse = "\n")
  plot(male$year, male$asmr_per_person_year, type = "l", lwd = 1.7,
       lty = 1, col = cols[i], xlab = "Year", ylab = "ASMR",
       ylim = yrange, xaxt = "n", main = title_text, cex.main = 0.85)
  axis(1, at = c(1980, 2000, 2020))
  lines(female$year, female$asmr_per_person_year,
        lwd = 1.7, lty = 2, col = cols[i])
  legend("topright", legend = c("Males", "Females"),
         col = cols[i], lty = c(1, 2), lwd = 1.7, bty = "n", cex = 0.68)
}
mtext("GBD Hungary, ages 60-89, 1980-2023", outer = TRUE,
      side = 3, line = 1.25, cex = 1.15, font = 2)
mtext("Age weights fixed to the pooled 1980 age distribution; not original HCSO Figure 2",
      outer = TRUE, side = 1, line = 1.0, cex = 0.72)
par(oldpar)
dev.off()

cat("GBD Figure 2 method replication completed.\n")
cat("Age weights sum:", sum(standard$weight), "\n")
cat("Max relative implied-population disagreement:", max_rel_error, "\n")
cat("All-cause-category ASMR rows:", nrow(asmr), "\n")
cat("Selected panels:", paste(top9, collapse = "; "), "\n")
cat("Saved outputs under:", out, "\n")
