# New portable plotting wrapper around existing, recomputed primary estimates.
# These four pages are not a facsimile regeneration of the original 60-page PDF.
source("R/bootstrap.R")
out <- "results/figures/main"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
pd <- "results/plot_data/main"; dir.create(pd, recursive = TRUE, showWarnings = FALSE)
ef <- read.csv("work/sensitivity_20261112/effects.csv")
ef <- ef[ef$type == "Primary" & ef$primary, ]
rk <- read.csv("work/sensitivity_20261112/rankings.csv"); rk <- rk[rk$type == "Primary", ]
co <- read.csv("work/sensitivity_20261112/components.csv"); co <- co[co$type == "Primary", ]
keys <- c("CAIT", "FAAM_ADL", "Sports", "Anterior", "Posteromedial", "Posterolateral")
display <- c(CAIT = "CAIT", FAAM_ADL = "FAAM-ADL", Sports = "Sports function", Anterior = "Anterior balance", Posteromedial = "Posteromedial balance", Posterolateral = "Posterolateral balance")
blue <- "#2D6484"; grey <- "#52606D"; orange <- "#B96834"
setup <- function() par(mfrow = c(3, 2), oma = c(1, 1, 1, 1), family = "Helvetica", ps = 10, mgp = c(1.9, .5, 0), tcl = -.2)
forest <- function(component = FALSE) {
  setup()
  for (key in keys) {
    a <- if (component) co[co$outcome == key, ] else ef[ef$outcome == key, ]
    labs <- if (component) gsub(":", " x ", a$term, fixed = TRUE) else a$treatment
    if (!component) {
      r <- rk[rk$outcome == key, ]; ord <- order(r$rank[match(a$treatment, r$treatment)])
      a <- a[ord, ]; labs <- labs[ord]
    }
    y <- rev(seq_len(nrow(a))); xr <- range(c(0, a$lower, a$upper), na.rm = TRUE)
    xr <- xr + c(-.06, .10) * max(.1, diff(xr))
    par(mar = c(3.3, 8.2, 2.5, 1.2))
    plot(NA, xlim = xr, ylim = c(.3, nrow(a) + .7), xlab = paste(if (key %in% keys[1:2]) "MD" else "Hedges g", "(95% CI)"), ylab = "", yaxt = "n", bty = "n", cex.axis = .8, cex.lab = .85)
    segments(xr[1], y, xr[2], y, col = "#E9EDF1", lwd = .6); abline(v = 0, lty = 2, col = "#8A929C", lwd = .8)
    axis(2, at = y, labels = labs, las = 1, tick = FALSE, cex.axis = .78)
    col <- if (component) ifelse(grepl(":", a$term), orange, blue) else blue
    segments(a$lower, y, a$upper, y, col = col, lwd = 1.2); points(a$estimate, y, pch = 16, col = col, cex = .78)
    mtext(display[key], side = 3, line = 1, adj = 0, font = 2, cex = 1)
  }
}
ranking <- function() {
  setup()
  for (key in keys) {
    a <- rk[rk$outcome == key, ]; a <- a[order(a$rank), ]; y <- rev(seq_len(nrow(a)))
    par(mar = c(3.3, 8.2, 2.5, 1.2))
    plot(NA, xlim = c(0, 1.12), ylim = c(.3, nrow(a) + .7), xlab = "P-score (higher ranks better)", ylab = "", yaxt = "n", xaxt = "n", bty = "n", cex.lab = .85)
    axis(1, at = seq(0, 1, .25), cex.axis = .8); axis(2, at = y, labels = a$treatment, las = 1, tick = FALSE, cex.axis = .78)
    segments(0, y, a$Pscore, y, col = "#DDE8EE", lwd = 5); points(a$Pscore, y, pch = 16, col = blue, cex = .8)
    text(1.03, y, sprintf("%.3f", a$Pscore), adj = 0, cex = .68, col = grey)
    mtext(display[key], side = 3, line = 1, adj = 0, font = 2, cex = 1)
  }
}
edges <- list(); nodes <- list()
network <- function() {
  setup()
  for (key in keys) {
    a <- read.csv(file.path("work/cinema_20261112/diagnostics", key, "analysed_arms.csv"))
    tt <- sort(unique(a$t)); theta <- seq(0, 2*pi, length.out = length(tt)+1)[-1] + pi/2
    x <- cos(theta); y <- sin(theta); n <- vapply(tt, function(t) sum(a$n[a$t == t]), numeric(1))
    par(mar = c(1.2, 1.2, 2.5, 1.2)); plot(NA, xlim = c(-1.5, 1.5), ylim = c(-1.4, 1.4), asp = 1, xlab = "", ylab = "", axes = FALSE)
    e <- list()
    for (ab in combn(tt, 2, simplify = FALSE)) {
      ids <- intersect(a$id[a$t == ab[1]], a$id[a$t == ab[2]])
      if (!length(ids)) next
      ij <- match(ab, tt); segments(x[ij[1]], y[ij[1]], x[ij[2]], y[ij[2]], col = "#B9C9D3", lwd = .8 + sqrt(length(ids)))
      text(mean(x[ij]), mean(y[ij]), length(ids), cex = .65, col = grey)
      e[[paste(ab, collapse = "|")]] <- data.frame(outcome = key, treatment1 = ab[1], treatment2 = ab[2], report_id_count = length(ids))
    }
    points(x, y, pch = 21, bg = ifelse(tt == "NAT", orange, blue), col = "white", lwd = 1, cex = 2 + 2*sqrt(n/max(n)))
    text(1.21*x, 1.21*y, tt, cex = .7, col = "#243B53")
    mtext(display[key], side = 3, line = 1, adj = 0, font = 2, cex = 1)
    edges[[key]] <<- do.call(rbind, e); nodes[[key]] <<- data.frame(outcome = key, treatment = tt, analysis_n_sum = n, x = x, y = y)
  }
}
pdf(file.path(out, "primary_NMA_CNMA_networks_rankings.pdf"), width = 12, height = 11.8, family = "Helvetica", pointsize = 10)
network(); forest(FALSE); ranking(); forest(TRUE); dev.off()
write.csv(ef, file.path(pd, "NMA_forest_values.csv"), row.names = FALSE)
write.csv(co, file.path(pd, "CNMA_forest_values.csv"), row.names = FALSE)
write.csv(rk, file.path(pd, "Pscore_values.csv"), row.names = FALSE)
write.csv(do.call(rbind, edges), file.path(pd, "network_edges.csv"), row.names = FALSE)
write.csv(do.call(rbind, nodes), file.path(pd, "network_nodes.csv"), row.names = FALSE)
cat("Primary plots: networks, NMA forests, P-scores, fixed-specification CNMA forests. Plot data exported.\n")
