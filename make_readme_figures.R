tr <- read.csv("results/model/sas_trajectories_apoe4.csv")
sl <- read.csv("results/model/sas_slopes_by_apoe4.csv")
col <- c("#3987E5", "#D95926", "#199E70")
png("figures/trajectories_by_apoe4.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 2.5, 1), family = "sans")
plot(NA, xlim = c(0, 6), ylim = range(tr$Lower, tr$Upper), xlab = "Years since baseline",
     ylab = "Predicted ADAS-Cog 13 (higher = worse)", las = 1, bty = "l",
     main = "Model-predicted ADAS-Cog 13 by APOE e4 allele count")
for (a in 0:2) {
  d <- tr[tr$apoe4 == a, ]
  polygon(c(d$time, rev(d$time)), c(d$Lower, rev(d$Upper)), col = adjustcolor(col[a + 1], 0.2), border = NA)
  lines(d$time, d$Estimate, col = col[a + 1], lwd = 3); points(d$time, d$Estimate, col = col[a + 1], pch = 19)
}
legend("topleft", legend = c("0 e4 alleles", "1 e4 allele", "2 e4 alleles"), col = col, lwd = 3, pch = 19, bty = "n")
dev.off()

d <- sl[grepl("^Slope: ", sl$Label), ]
d$a <- as.integer(substr(sub("Slope: ", "", d$Label), 1, 1))
an <- c("MAR (main)", "Complete-case", "Pattern-mixture")
png("figures/slopes_sensitivity.png", width = 1600, height = 900, res = 200)
par(mar = c(4.5, 9, 2.5, 1))
y <- 0
plot(NA, xlim = c(0, max(d$Upper) * 1.05), ylim = c(0.5, 11.5), yaxt = "n", xlab = "Annual change in ADAS-Cog 13 (points/year, 95% CI)",
     ylab = "", bty = "l", main = "APOE e4 slopes under three missing-data assumptions")
labs <- c(); at <- c()
for (i in seq_along(an)) for (a in 0:2) {
  r <- d[d$analysis == an[i] & d$a == a, ]; yy <- (3 - i) * 4 + (3 - a)
  segments(r$Lower, yy, r$Upper, yy, col = col[a + 1], lwd = 3); points(r$Estimate, yy, pch = 19, col = col[a + 1], cex = 1.3)
}
axis(2, at = c(10, 6, 2), labels = an, las = 1, tick = FALSE)
legend("topright", legend = c("0 e4", "1 e4", "2 e4"), col = col, pch = 19, lwd = 3, bty = "n")
dev.off()
