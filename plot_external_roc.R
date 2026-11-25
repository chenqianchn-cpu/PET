# Run from this directory after external_validation.R.
performance <- read.delim("results/performance.tsv")
roc_all <- read.delim("results/GSE14520_paired_ROC_coordinates.tsv")
roc_specific <- read.delim("results/GSE14520_paired_ROC_without_marker_x_at.tsv")
main <- performance[performance$analysis == "GSE14520 paired primary", ]
sensitivity <- performance[
  performance$analysis == "GSE14520 paired sensitivity without marker x_at probes", ]
stopifnot(nrow(main) == 1L, nrow(sensitivity) == 1L,
          all(diff(roc_all$fpr) >= 0), all(diff(roc_specific$fpr) >= 0),
          main$patients == sensitivity$patients)
summary_line <- function(row) {
  sprintf("%.3f (%.3f to %.3f)", row$estimate, row$ci_low, row$ci_high)
}
primary_color <- "#006B79"
sensitivity_color <- "#855373"

draw_roc <- function() {
  par(mar = c(3.65, 3.7, 3.0, 0.9), mgp = c(2.3, 0.65, 0),
      tcl = -0.24, xaxs = "i", yaxs = "i", las = 1, fg = "#252C34")
  plot(NA_real_, NA_real_, xlim = c(0, 1), ylim = c(0, 1),
       xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  abline(0, 1, lty = 3, col = "#AAB0B7", lwd = 0.85)
  lines(roc_specific$fpr, roc_specific$tpr, col = sensitivity_color, lwd = 1.6, lty = 2)
  lines(roc_all$fpr, roc_all$tpr, col = primary_color, lwd = 1.8)
  ticks <- seq(0, 1, 0.25)
  labels <- c("0", "0.25", "0.50", "0.75", "1.00")
  axis(1, at = ticks, labels = FALSE, lwd = 0.65)
  axis(2, at = ticks, labels = FALSE, lwd = 0.65)
  # Separate the label drawing operations to preserve individual PDF text boxes.
  for (i in seq_along(ticks)) {
    axis(1, at = ticks[i], labels = labels[i], cex.axis = 0.87, tick = FALSE,
         col.axis = "#252C34")
    axis(2, at = ticks[i], labels = labels[i], cex.axis = 0.87, tick = FALSE,
         col.axis = "#262D35")
  }
  title(xlab = "False-positive rate", ylab = "True-positive rate", cex.lab = 0.95)
  mtext("GSE14520 external validation", side = 3, line = 1.65, adj = 0, font = 2, cex = 1.08)
  mtext(sprintf("%d patient-matched pairs", main$patients), side = 3, line = 0.55,
        adj = 0, cex = 0.84)
  x <- 0.47
  text(x, 0.38, "AUROC (95% CI)", adj = 0, cex = 0.80, font = 2)
  segments(x, 0.295, x + 0.075, 0.295, col = primary_color, lwd = 1.8)
  text(x + 0.10, 0.295, "All mapped probes", adj = 0, cex = 0.81, col = primary_color)
  text(x + 0.10, 0.224, summary_line(main), adj = 0, cex = 0.78)
  segments(x, 0.148, x + 0.075, 0.148, col = sensitivity_color, lwd = 1.6, lty = 2)
  text(x + 0.10, 0.148, "Exclude marker _x_at", adj = 0, cex = 0.81, col = sensitivity_color)
  text(x + 0.10, 0.078, summary_line(sensitivity), adj = 0, cex = 0.78)
}
width_mm <- 90
w <- width_mm / 25.4
h <- w
cairo_pdf("results/GSE14520_external_ROC.pdf", width = w, height = h,
          pointsize = 9, family = "Arial")
draw_roc()
invisible(dev.off())
png("results/GSE14520_external_ROC.png", width = w, height = h, units = "in",
    res = 600, pointsize = 9, type = "cairo", family = "Arial", bg = "white")
draw_roc()
invisible(dev.off())
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite("results/GSE14520_external_ROC.svg", width = w, height = h,
                   pointsize = 9, system_fonts = list(sans = "Arial"))
  draw_roc()
  invisible(dev.off())
}
