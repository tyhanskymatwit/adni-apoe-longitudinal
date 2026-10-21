# ==============================================================================
# Program : adni_tlf_qc.R
# Purpose : Independent QC (double programming) of the TLF tables produced by
#           adni_tlf.sas. Re-derives Tables 14.1.1, 14.1.2, 14.2.1 and Listing
#           16.2.1 from adni_wide.csv in base R, without reading any SAS
#           intermediate data, then compares every numeric cell with the SAS
#           exports in output/TLF/qc_*.csv.
# Run     : Rscript tlf/adni_tlf_qc.R   (from the ADNI project folder)
# Output  : output/TLF/r_t*.csv          R-derived tables
#           output/TLF/qc_report.txt     cell-by-cell comparison verdict
# ==============================================================================

args   <- commandArgs(trailingOnly = TRUE)
root   <- if (length(args)) args[1] else "."
tlfdir <- file.path(root, "output", "TLF")
dir.create(tlfdir, showWarnings = FALSE, recursive = TRUE)
# SAS exports: output/TLF, or tlf/output if they were downloaded there instead
sasdir <- if (length(args) > 1) args[2] else
  Find(function(d) file.exists(file.path(d, "qc_t1411.csv")),
       c(tlfdir, file.path(root, "tlf", "output")), nomatch = tlfdir)
window <- 0.35

# SAS put(x, w.d) rounds half away from zero; match it so strings agree
fmt <- function(x, d) {
  ifelse(is.na(x), ".",
         formatC(sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-9) / 10^d, format = "f", digits = d))
}

# ---- analysis set (same rule as adni_tlf.sas) --------------------------------
w <- read.csv(file.path(root, "data", "adni_wide.csv"), stringsAsFactors = FALSE)
w$dxn <- match(w$dx_bl, c("CN", "MCI", "DEM"))
adsl  <- w[complete.cases(w[, c("age_bl", "female", "educ", "dxn", "apoe4", "adas13_0")]), ]
adsl$trtn  <- adsl$apoe4
adsl$protn <- match(adsl$protocol, c("ADNI1", "ADNIGO", "ADNI2", "ADNI3", "ADNI4"))

y  <- as.matrix(adsl[, paste0("adas13_", 0:6)])
el <- sapply(0:6, function(k) adsl$fu_potential >= k + window | !is.na(y[, k + 1]))
el[is.na(el)] <- FALSE
adsl$lastobs  <- apply(y, 1, function(r) max(which(!is.na(r))) - 1)
adsl$lastelig <- apply(el, 1, function(r) max(which(r)) - 1)
adsl$lastval  <- y[cbind(seq_len(nrow(y)), adsl$lastobs + 1)]
adsl$dropout  <- as.integer(adsl$lastobs < adsl$lastelig)
adsl$dropgrp  <- with(adsl, ifelse(dropout == 0, 0, ifelse(lastobs <= 2, 1, 2)))

arms  <- c(0, 1, 2, 3)                                   # 3 = Total
inarm <- function(d, a) if (a == 3) d else d[d$trtn == a, ]
bign  <- sapply(arms, function(a) nrow(inarm(adsl, a)))
cat("Big N:", bign, "\n")

# ---- row builders -------------------------------------------------------------
cont_rows <- function(var, ord, dec) {
  cells <- sapply(arms, function(a) {
    x <- inarm(adsl, a)[[var]]; x <- x[!is.na(x)]
    c(as.character(length(x)),
      paste0(fmt(mean(x), dec + 1), " (", fmt(sd(x), dec + 2), ")"),
      fmt(median(x), dec + 1),
      paste0(fmt(min(x), dec), ", ", fmt(max(x), dec)))
  })
  data.frame(ord = ord, sub = 1:4, cells, check.names = FALSE)
}
cat_rows <- function(var, ord, levels) {
  x_all <- adsl[[var]]
  lv <- levels
  if (any(is.na(x_all))) lv <- c(lv, NA)
  cells <- sapply(arms, function(a) {
    d <- inarm(adsl, a); x <- d[[var]]
    sapply(lv, function(l) {
      n <- if (is.na(l)) sum(is.na(x)) else sum(x == l, na.rm = TRUE)
      paste0(n, " (", fmt(100 * n / nrow(d), 1), ")")
    })
  })
  cells <- matrix(cells, nrow = length(lv))
  data.frame(ord = ord, sub = ifelse(is.na(lv), 99, lv), cells, check.names = FALSE)
}
name_cols <- function(df, prefix = "trt") { names(df)[-(1:2)] <- paste0(prefix, arms); df }

# ---- Table 14.1.1 -------------------------------------------------------------
t1411 <- do.call(rbind, lapply(list(
  cont_rows("age_bl", 1, 1), cat_rows("female", 2, 0:1), cont_rows("educ", 3, 0),
  cat_rows("dxn", 4, 1:3), cat_rows("amyloid", 5, 0:1), cat_rows("protn", 6, 1:5),
  cont_rows("adas13_0", 7, 2), cont_rows("mmse_0", 8, 0), cont_rows("cdrsb_0", 9, 1)),
  name_cols))

# ---- Table 14.1.2 -------------------------------------------------------------
ret <- do.call(rbind, lapply(1:6, function(k) {
  cells <- sapply(arms, function(a) {
    idx <- if (a == 3) TRUE else adsl$trtn == a
    e <- el[idx, k + 1]; o <- !is.na(y[idx, k + 1]) & e
    paste0(sum(o), "/", sum(e), " (", fmt(100 * sum(o) / sum(e), 1), ")")
  })
  data.frame(ord = 3, sub = k, t(cells))
}))
t1412 <- do.call(rbind, lapply(list(cat_rows("dropgrp", 1, 0:2), cat_rows("died", 2, 0:1), ret),
                               name_cols))

# ---- Table 14.2.1 -------------------------------------------------------------
t1421 <- do.call(rbind, lapply(0:6, function(k) {
  cols <- list()
  for (a in 0:2) for (p in 1:2) {
    idx <- adsl$trtn == a & !is.na(y[, k + 1])
    v <- if (p == 1) y[idx, k + 1] else y[idx, k + 1] - y[idx, 1]
    cols[[paste0("c", a, p)]] <- if (p == 2 && k == 0 || length(v) == 0) rep(NA, 4) else
      c(as.character(length(v)), paste0(fmt(mean(v), 3), " (", fmt(sd(v), 4), ")"),
        fmt(median(v), 3), paste0(fmt(min(v), 2), ", ", fmt(max(v), 2)))
  }
  data.frame(ord = k, sub = 1:4, cols)
}))

# ---- Listing 16.2.1 -----------------------------------------------------------
l1621 <- adsl[adsl$apoe4 == 2 & adsl$dropout == 1,
              c("rid", "protocol", "age_bl", "female", "dxn", "adas13_0",
                "lastobs", "lastelig", "lastval", "died")]
l1621 <- l1621[order(l1621$rid), ]

for (nm in c("t1411", "t1412", "t1421", "l1621"))
  write.csv(get(nm), file.path(tlfdir, paste0("r_", nm, ".csv")), row.names = FALSE, na = "")

# ---- compare with SAS ---------------------------------------------------------
nums <- function(s) {
  s <- ifelse(is.na(s), "", as.character(s))
  m <- regmatches(s, gregexpr("-?[0-9]+(\\.[0-9]+)?", s))
  lapply(m, function(v) list(val = as.numeric(v),
                             dec = nchar(sub("^[^.]*\\.?", "", v))))
}
compare_tab <- function(r, sasfile, keys) {
  s <- read.csv(sasfile, stringsAsFactors = FALSE, colClasses = "character")
  for (k in keys) s[[k]] <- as.numeric(s[[k]])
  vcols <- setdiff(names(r), keys)
  m <- merge(r, s, by = keys, all = TRUE, suffixes = c(".r", ".sas"))
  issues <- character()
  if (nrow(m) != nrow(r) || nrow(m) != nrow(s))
    issues <- c(issues, sprintf("row count R=%d SAS=%d merged=%d", nrow(r), nrow(s), nrow(m)))
  for (v in vcols) {
    a <- nums(m[[paste0(v, ".r")]]); b <- nums(m[[paste0(v, ".sas")]])
    for (i in seq_len(nrow(m))) {
      ok <- length(a[[i]]$val) == length(b[[i]]$val) &&
        all(abs(a[[i]]$val - b[[i]]$val) <= 1.01 * 10^-pmax(a[[i]]$dec, b[[i]]$dec))
      if (!ok) issues <- c(issues, sprintf("%s=%s %s=%s | %s: R='%s' SAS='%s'",
        keys[1], m[i, keys[1]], keys[2], m[i, keys[2]], v,
        m[i, paste0(v, ".r")], m[i, paste0(v, ".sas")]))
    }
  }
  issues
}

report <- c("ADNI TLF independent QC (R vs SAS)", paste("Run:", format(Sys.time())),
            paste("SAS exports read from:", normalizePath(sasdir)),
            paste("Big N (0/1/2/Total e4):", paste(bign, collapse = " / ")), "")
checks <- list(t1411 = c("ord", "sub"), t1412 = c("ord", "sub"), t1421 = c("ord", "sub"))
for (nm in names(checks)) {
  f <- file.path(sasdir, paste0("qc_", nm, ".csv"))
  if (!file.exists(f)) { report <- c(report, sprintf("%-6s NOT RUN (SAS export %s not found)", nm, basename(f))); next }
  iss <- compare_tab(get(nm), f, checks[[nm]])
  report <- c(report, sprintf("%-6s %s (%d issues)", nm, if (length(iss)) "MISMATCH" else "MATCH", length(iss)),
              if (length(iss)) paste("   ", head(iss, 25)))
}
f <- file.path(sasdir, "qc_l1621.csv")
if (file.exists(f)) {
  s <- read.csv(f); s <- s[order(s$rid), ]
  same <- nrow(s) == nrow(l1621) && all(s$rid == l1621$rid) &&
    isTRUE(all.equal(s$lastval, l1621$lastval)) && all(s$lastobs == l1621$lastobs)
  report <- c(report, sprintf("l1621  %s (SAS %d rows, R %d rows)", if (same) "MATCH" else "MISMATCH",
                              nrow(s), nrow(l1621)))
} else report <- c(report, "l1621  NOT RUN (SAS export qc_l1621.csv not found)")

writeLines(report, file.path(tlfdir, "qc_report.txt"))
cat(report, sep = "\n")
