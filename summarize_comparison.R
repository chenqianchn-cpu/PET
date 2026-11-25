# Run after both feature-panel runs, from the parent external-validation directory.
proxy_dir <- "ML15/results"
strict_dir <- "ML15/results_strict14"
rank_summary <- read.delim(file.path(proxy_dir, "gene_rank_summary.tsv"))
strict <- read.delim(file.path(strict_dir, "gene_rank_summary.tsv"))
modes <- c("TCGA grouped CV", "GSE14520 grouped CV", "TCGA to GSE14520",
           "TCGA to GSE14520 without PPIA x_at")
stopifnot(setequal(rank_summary$analysis, modes), setequal(strict$analysis, modes),
          length(unique(rank_summary$gene)) == 15L, length(unique(strict$gene)) == 14L)
collect_columns <- function(data, tag) {
  do.call(cbind, lapply(seq_along(modes), function(i) {
    mode <- data[data$analysis == modes[i], ]
    index <- match(gene_list, mode$gene)
    result <- mode[index, c("mean_rank", "top_5", "models")]
    names(result) <- paste(tag, c("TCGA_CV", "GSE_CV", "Transfer", "Transfer_without_PPIA_x_at")[i],
                           names(result), sep = "_")
    result
  }))
}
gene_list <- sort(unique(rank_summary$gene))
output <- cbind(data.frame(gene = gene_list), collect_columns(rank_summary, "proxy15"),
                collect_columns(strict, "strict14"))
output <- output[order(output$proxy15_GSE_CV_mean_rank, output$gene), ]
write.table(output, "ML15/combined_gene_ranks.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

ranks <- read.delim(file.path(proxy_dir, "feature_ranks_by_model.tsv"))
gse <- ranks[ranks$analysis == "GSE14520 grouped CV", ]
models <- unique(gse$model)
display_genes <- output$gene
stopifnot(length(models) == 10L, nrow(gse) == length(display_genes) * length(models),
          !anyDuplicated(gse[c("gene", "model")]))
values <- xtabs(rank ~ gene + model, data = gse)[display_genes, models]
palette <- grDevices::colorRampPalette(c("#075C63", "#EBF3F1"))(141)
rank_color <- function(rank) palette[1L + round((rank - 1) * 10)]
draw_heatmap <- function() {
  par(mar = rep(0.45, 4), family = "Arial", xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  x <- seq(0.131, 0.875, length.out = length(models) + 1L)
  y <- seq(0.850, 0.160, length.out = length(display_genes) + 1L)
  for (i in seq_along(display_genes)) {
    rect(head(x, -1L), y[i + 1L], tail(x, -1L), y[i],
         col = rank_color(values[i, ]), border = "white", lwd = 0.45)
  }
  axis_ink <- c("#252C34", "#262D35")
  for (i in seq_along(display_genes)) {
    label <- if (display_genes[i] == "ANG") "ANG proxy" else display_genes[i]
    text(x[1L] - 0.013, mean(y[i + 0:1]), label, adj = c(1, 0.5),
         cex = 0.94, col = axis_ink[1L + (i %% 2L)])
  }
  model_labels <- models
  model_labels[models == "NaiveBayes"] <- "Naive Bayes"
  model_labels[models == "ElasticNet"] <- "Elastic net"
  for (i in seq_along(models)) {
    text(mean(x[i + 0:1]), 0.121, model_labels[i], adj = c(0.5, 0.5),
         cex = 0.88, col = axis_ink[1L + (i %% 2L)])
  }
  text(x[1L], 0.96, "Feature ranks in GSE14520", adj = 0,
       font = 2, cex = 1.40, col = "#252C34")
  text(x[1L], 0.909, "Patient-grouped cross-validation, 214 matched pairs", adj = 0,
       cex = 0.92, col = "#252C34")
  bar_x <- c(0.898, 0.923)
  bar_y <- seq(0.33, 0.70, length.out = length(palette) + 1L)
  rect(bar_x[1L], head(bar_y, -1L), bar_x[2L], tail(bar_y, -1L),
       col = rev(palette), border = NA)
  text(bar_x[1L], max(bar_y) + 0.043, "Rank", adj = 0, font = 2, cex = 0.88)
  ticks <- c(1L, 5L, 10L, 15L)
  for (i in seq_along(ticks)) {
    pos <- max(bar_y) - (ticks[i] - 1L) * diff(range(bar_y)) / 14
    text(bar_x[2L] + 0.017, pos, ticks[i], adj = 0, cex = 0.84,
         col = axis_ink[1L + (i %% 2L)])
  }
}
width_mm <- 183
height_mm <- 113
figure_name <- file.path(proxy_dir, "GSE14520_rank_heatmap")
w <- width_mm / 25.4
h <- height_mm / 25.4
grDevices::cairo_pdf(paste0(figure_name, ".pdf"), width = w, height = h,
                    pointsize = 9, family = "Arial")
draw_heatmap()
invisible(dev.off())
png(paste0(figure_name, ".png"), width = w, height = h, units = "in",
    res = 600, type = "cairo", pointsize = 9, bg = "white")
draw_heatmap()
invisible(dev.off())
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(paste0(figure_name, ".svg"), width = w, height = h,
                   pointsize = 9, bg = "white", system_fonts = list(sans = "Arial"))
  draw_heatmap()
  invisible(dev.off())
}

print(output[, c("gene", "proxy15_TCGA_CV_mean_rank", "proxy15_GSE_CV_mean_rank",
                 "proxy15_Transfer_mean_rank", "strict14_GSE_CV_mean_rank")], row.names = FALSE)
