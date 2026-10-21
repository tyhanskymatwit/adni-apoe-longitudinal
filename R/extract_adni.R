pkg_data <- file.path("ADNIMERGE2", "data")
out_dir  <- "data"
dir.create(out_dir, showWarnings = FALSE)

load_rda <- function(name) {
  e <- new.env()
  # ADNIMERGE2 objects carry xportr/metacore attributes; the namespace warning
  # on load is harmless.
  suppressWarnings(load(file.path(pkg_data, paste0(name, ".rda")), envir = e))
  as.data.frame(get(name, envir = e))
}

ADSL   <- load_rda("ADSL")
ADQS   <- load_rda("ADQS")
cutoff <- as.Date(load_rda("DATA_DOWNLOADED_DATE")[[1]])

MAX_YEAR <- 6                    # nominal annual visits 0..6
OUTCOMES <- c(adas13 = "ADASTT13",  # ADAS-Cog 13 (0-85, higher = worse)
              mmse   = "MMSCORE",   # MMSE (0-30, higher = better)
              cdrsb  = "CDRSB")     # CDR Sum of Boxes (0-18, higher = worse)


sl <- ADSL[ADSL$ENRLFL %in% "Y", ]

# APOE e4 dose: count the "4" alleles in genotype strings such as "e3/e4"
apoe_chr <- as.character(sl$APOE)
apoe4 <- ifelse(is.na(apoe_chr), NA_integer_,
                lengths(regmatches(apoe_chr, gregexpr("4", apoe_chr))))

enrl_years <- function(d) as.numeric(d - sl$ENRLDT) / 365.25

base <- data.frame(
  rid       = as.integer(sl$SUBJID),
  usubjid   = sl$USUBJID,
  protocol  = as.character(sl$ORIGPROT),        # ADNI1 / ADNIGO / ADNI2 / ADNI3 / ADNI4
  age_bl    = round(sl$AGE, 2),
  female    = as.integer(sl$SEX == "Female"),
  educ      = sl$EDUC,
  dx_bl     = as.character(sl$DX),               # CN / MCI / DEM
  apoe_geno = gsub("ε", "e", apoe_chr),     # e3/e4 etc. (ASCII for Stata/SAS)
  apoe4     = apoe4,                             # 0 / 1 / 2 copies of e4
  amyloid   = ifelse(is.na(sl$AMYSTAT), NA_integer_,
                     as.integer(sl$AMYSTAT == "Elevated")),
  enrl_date = format(sl$ENRLDT, "%Y-%m-%d"),
  # years the participant COULD have been followed by the data cut-off:
  # used to separate administrative censoring from true dropout
  fu_potential = round(as.numeric(cutoff - sl$ENRLDT) / 365.25, 3),
  died      = as.integer(sl$DTHFL %in% "Yes"),
  eos_years = round(enrl_years(sl$EOSDT), 3),    # years to end-of-study (if ended)
  stringsAsFactors = FALSE
)


qs <- ADQS[ADQS$ENRLFL %in% "Y" & ADQS$PARAMCD %in% OUTCOMES &
           !is.na(ADQS$AVISITN) & !is.na(ADQS$AVAL),
           c("USUBJID", "SUBJID", "PARAMCD", "AVISIT", "AVISITN", "ADT", "AVAL")]
qs$years <- qs$AVISITN / 365.25
qs$year  <- round(qs$years)
qs$dist  <- abs(qs$years - qs$year)

write.csv(data.frame(rid = as.integer(qs$SUBJID), paramcd = qs$PARAMCD,
                     avisit = qs$AVISIT, days = qs$AVISITN,
                     years = round(qs$years, 4), aval = qs$AVAL),
          file.path(out_dir, "adni_visits_long.csv"), row.names = FALSE, na = "")

# Visit window: +/- 0.35 years (~4 months) around each anniversary. This
# drops semi-annual visits (m06, m18, ...) and off-cycle protocol-transition
# visits, whose deviation from the annual grid is ~0.5 years.
WINDOW <- 0.35
qs <- qs[qs$year <= MAX_YEAR & qs$dist <= WINDOW, ]
qs <- qs[order(qs$USUBJID, qs$PARAMCD, qs$year, qs$dist), ]
qs <- qs[!duplicated(qs[, c("USUBJID", "PARAMCD", "year")]), ]

# One time variable per nominal year. The ADAS-Cog visit time is used
# (MMSE/CDR are given at the same visits); fall back to MMSE when missing.
pick_time <- function(code) {
  x <- qs[qs$PARAMCD == code, c("USUBJID", "year", "years")]
  names(x)[3] <- code; x
}
tt <- merge(pick_time("ADASTT13"), pick_time("MMSCORE"),
            by = c("USUBJID", "year"), all = TRUE)
tt$t <- ifelse(is.na(tt$ADASTT13), tt$MMSCORE, tt$ADASTT13)

to_wide <- function(df, value, prefix) {
  w <- reshape(df[, c("USUBJID", "year", value)], idvar = "USUBJID",
               timevar = "year", direction = "wide", sep = "_")
  names(w) <- sub(paste0("^", value, "_"), prefix, names(w))
  # guarantee all columns 0..MAX_YEAR exist and are in order
  for (k in 0:MAX_YEAR) {
    nm <- paste0(prefix, k); if (!nm %in% names(w)) w[[nm]] <- NA
  }
  w[, c("USUBJID", paste0(prefix, 0:MAX_YEAR))]
}

wide <- to_wide(data.frame(USUBJID = tt$USUBJID, year = tt$year,
                           t = round(tt$t, 4)), "t", "t_")
for (nm in names(OUTCOMES)) {
  d <- qs[qs$PARAMCD == OUTCOMES[[nm]], c("USUBJID", "year", "AVAL")]
  wide <- merge(wide, to_wide(d, "AVAL", paste0(nm, "_")),
                by = "USUBJID", all = TRUE)
}


adni <- merge(base, wide, by.x = "usubjid", by.y = "USUBJID", all.x = TRUE)
adni <- adni[order(adni$rid), ]
stopifnot(!anyDuplicated(adni$rid))
adni$usubjid <- NULL

write.csv(adni, file.path(out_dir, "adni_wide.csv"), row.names = FALSE, na = "")

codebook <- data.frame(
  variable = c("rid", "protocol", "age_bl", "female", "educ", "dx_bl",
               "apoe_geno", "apoe4", "amyloid", "enrl_date", "fu_potential",
               "died", "eos_years", "t_0..t_6", "adas13_0..adas13_6",
               "mmse_0..mmse_6", "cdrsb_0..cdrsb_6"),
  description = c(
    "ADNI subject ID (RID)",
    "Original ADNI protocol at enrolment (ADNI1, ADNIGO, ADNI2, ADNI3, ADNI4)",
    "Age at baseline (years)",
    "1 = female, 0 = male",
    "Years of education",
    "Baseline diagnosis: CN, MCI, DEM",
    "APOE genotype (e.g. e3/e4)",
    "Number of APOE e4 alleles (0, 1, 2)",
    "Baseline amyloid PET status: 1 = elevated, 0 = not elevated",
    "Enrolment date",
    paste0("Years between enrolment and data cut-off (", cutoff,
           "): maximum possible follow-up"),
    "1 = death recorded",
    "Years from enrolment to end-of-study date (missing if still active)",
    "Actual time (years since baseline) of the visit mapped to nominal year k",
    "ADAS-Cog 13-item total at nominal year k (0-85, higher = worse)",
    "MMSE at nominal year k (0-30, higher = better)",
    "CDR Sum of Boxes at nominal year k (0-18, higher = worse)"))
write.csv(codebook, file.path(out_dir, "codebook.csv"), row.names = FALSE)

cat("Data cut-off:", format(cutoff), "\n")
cat("Participants written:", nrow(adni), "\n")
cat("With baseline ADAS-Cog13:", sum(!is.na(adni$adas13_0)), "\n")
cat("Complete baseline covariates + ADAS13_0:",
    sum(complete.cases(adni[, c("age_bl", "female", "educ", "dx_bl",
                                "apoe4", "adas13_0")])), "\n")
cat("Observed ADAS-Cog13 by nominal year:\n")
print(colSums(!is.na(adni[, paste0("adas13_", 0:MAX_YEAR)])))
