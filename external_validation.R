# Run from this directory with Rscript --vanilla external_validation.R.
# Only base R is required. The source URLs and checksums are in source_manifest.tsv.
source("prepare_expression.R")
source("validation_metrics.R")

# Fixed, label-blind reference panel. Missing measurements are checked using TCGA only.
# GEO expression is not used to choose reference genes, fit the model or set a threshold.
complete_tcga <- rowSums(!is.finite(tcga)) == 0L
reference <- sort(setdiff(
  intersect(rownames(tcga)[complete_tcga], names(probe_groups)), genes
), method = "radix")
stopifnot(
  length(reference) > 1000L, all(genes %in% names(probe_groups)),
  all(genes %in% rownames(tcga)[complete_tcga]),
  all(is.finite(geo_by_gene[c(reference, genes), ]))
)
sample_percentiles <- function(values) {
  ref <- values[reference, , drop = FALSE]
  result <- vapply(genes, function(gene) {
    target <- values[gene, ]
    (colSums(sweep(ref, 2L, target, "<")) +
       0.5 * colSums(sweep(ref, 2L, target, "=="))) / nrow(ref)
  }, FUN.VALUE = numeric(ncol(values)))
  rownames(result) <- colnames(values)
  stopifnot(all(is.finite(result)), all(result >= 0), all(result <= 1))
  as.data.frame(result)
}
make_features <- function(meta, values) {
  stopifnot(identical(meta$sample, colnames(values)))
  cbind(meta, sample_percentiles(values))
}
tcga_data <- make_features(tcga_meta, tcga)
geo_data <- make_features(geo_meta, geo_by_gene)
fit <- glm(tumor ~ PPIA + IMPDH2, family = binomial(), data = tcga_data)
stopifnot(isTRUE(fit$converged), all(is.finite(coef(fit))))
weights <- coef(fit)
linear_score <- function(ppia, impdh2) {
  unname(weights["(Intercept)"] + weights["PPIA"] * ppia + weights["IMPDH2"] * impdh2)
}
tcga_data$score <- linear_score(tcga_data$PPIA, tcga_data$IMPDH2)
geo_data$score <- linear_score(geo_data$PPIA, geo_data$IMPDH2)
stopifnot(isTRUE(all.equal(geo_data$score, unname(predict(fit, newdata = geo_data, type = "link")))))

# Restrict the marker measurement to probes without Affymetrix's cross-hybridization
# suffix as a sensitivity check. The original reference and fitted model stay fixed.
specific_mapping <- mapping[mapping$gene %in% genes & !endsWith(mapping$probe, "_x_at"), ]
specific_probes <- split(specific_mapping$probe, specific_mapping$gene)
stopifnot(setequal(names(specific_probes), genes))
geo_marker_sensitivity <- geo_by_gene
for (gene in genes) {
  probes <- specific_probes[[gene]]
  geo_marker_sensitivity[gene, ] <- apply(geo[probes, , drop = FALSE], 2L, median)
}
geo_marker_sensitivity <- make_features(geo_meta, geo_marker_sensitivity)
geo_marker_sensitivity$score <- linear_score(geo_marker_sensitivity$PPIA, geo_marker_sensitivity$IMPDH2)
stopifnot(identical(geo_data$sample, geo_marker_sensitivity$sample))

geo_pairs <- align_pairs(geo_data)
tcga_pairs <- align_pairs(tcga_data)

paired_concordance <- function(neg, pos) mean((pos > neg) + 0.5 * (pos == neg))

# Perturb one feature at a time by randomly swapping its values within each pair.
# Model coefficients and the other feature stay fixed. Use the same swap draws for both genes.
paired_permutation_losses <- function(pairs, bootstrap_n = 2000L, draws = 512L,
                                      bootstrap_draws = 64L, seed = 14520L) {
  set.seed(seed)
  score_pos <- pairs$pos$score
  score_neg <- pairs$neg$score
  pos_features <- as.matrix(pairs$pos[genes])
  neg_features <- as.matrix(pairs$neg[genes])
  n <- length(score_pos)
  get_losses <- function(idx, n_draws) {
    pos <- score_pos[idx]
    neg <- score_neg[idx]
    observed <- auc(neg, pos)
    sw <- matrix(sample(c(FALSE, TRUE), length(idx) * n_draws, replace = TRUE),
                 nrow = length(idx))
    pos_markers <- pos_features[idx, , drop = FALSE]
    neg_markers <- neg_features[idx, , drop = FALSE]
    observed - vapply(seq_along(genes), function(gene) {
      mean(vapply(seq_len(n_draws), function(j) {
        pos_permuted <- pos_markers
        neg_permuted <- neg_markers
        swap <- sw[, j]
        pos_permuted[swap, gene] <- neg_markers[swap, gene]
        neg_permuted[swap, gene] <- pos_markers[swap, gene]
        auc(linear_score(neg_permuted[, 1L], neg_permuted[, 2L]),
            linear_score(pos_permuted[, 1L], pos_permuted[, 2L]))
      }, numeric(1)))
    }, numeric(1))
  }
  point <- get_losses(seq_len(n), draws)
  boot <- replicate(bootstrap_n, get_losses(sample.int(n, n, replace = TRUE), bootstrap_draws))
  intervals <- t(apply(boot, 1L, ci95))
  data.frame(
    gene = genes, delta_auc = point, ci_low = intervals[, 1L], ci_high = intervals[, 2L],
    pairs = n, bootstrap_replicates = bootstrap_n, permutation_draws = draws,
    permutations_per_bootstrap = bootstrap_draws
  )
}

boot_n <- 4000L
set.seed(20261124L)
n <- nrow(geo_pairs$pos)
pos <- geo_pairs$pos$score
neg <- geo_pairs$neg$score
paired_boot <- replicate(boot_n, {
  idx <- sample.int(n, n, replace = TRUE)
  c(auc(neg[idx], pos[idx]), paired_concordance(neg[idx], pos[idx]))
})
# All platform samples sensitivity analysis: resample patients, keep all their tissues.
clusters <- split(seq_len(nrow(geo_data)), geo_data$patient)
set.seed(20261125L)
full_boot <- replicate(boot_n, {
  idx <- unlist(sample(clusters, length(clusters), replace = TRUE), use.names = FALSE)
  y <- geo_data$tumor[idx]
  auc(geo_data$score[idx[y == 0L]], geo_data$score[idx[y == 1L]])
})
specific_pairs <- align_pairs(geo_marker_sensitivity)
stopifnot(identical(specific_pairs$pos$patient, geo_pairs$pos$patient))
set.seed(20261126L)
specific_boot <- replicate(boot_n, {
  idx <- sample.int(n, n, replace = TRUE)
  auc(specific_pairs$neg$score[idx], specific_pairs$pos$score[idx])
})
specific_auc <- auc(specific_pairs$neg$score, specific_pairs$pos$score)
cohort_auc <- function(data) auc(data$score[data$tumor == 0L], data$score[data$tumor == 1L])
primary_auc <- auc(neg, pos)
roc_curve <- roc(neg, pos)
trapezoid_auc <- sum(diff(roc_curve$fpr) *
                      (head(roc_curve$tpr, -1L) + tail(roc_curve$tpr, -1L)) / 2)
stopifnot(isTRUE(all.equal(primary_auc, trapezoid_auc, tolerance = 1e-12)))
metrics <- data.frame(
  analysis = c("TCGA-LIHC apparent training", "GSE14520 paired primary",
               "GSE14520 all samples sensitivity", "GSE14520 within-patient concordance"),
  metric = c("AUROC", "AUROC", "AUROC", "Concordance"),
  estimate = c(cohort_auc(tcga_data), primary_auc, cohort_auc(geo_data),
               paired_concordance(neg, pos)),
  rbind(c(NA_real_, NA_real_), ci95(paired_boot[1L, ]), ci95(full_boot),
        ci95(paired_boot[2L, ])),
  patients = c(length(unique(tcga_data$patient)), n, length(clusters), n),
  tumor = c(sum(tcga_data$tumor), n, sum(geo_data$tumor), n),
  non_tumor = c(sum(tcga_data$tumor == 0L), n, sum(geo_data$tumor == 0L), n),
  bootstrap_replicates = c(0L, boot_n, boot_n, boot_n)
)
names(metrics)[4:5] <- c("ci_low", "ci_high")
metrics <- rbind(metrics, data.frame(
  analysis = "GSE14520 paired sensitivity without marker x_at probes",
  metric = "AUROC", estimate = specific_auc, ci_low = ci95(specific_boot)[1L],
  ci_high = ci95(specific_boot)[2L], patients = n, tumor = n, non_tumor = n,
  bootstrap_replicates = boot_n
))
geo_importance <- paired_permutation_losses(geo_pairs)
tcga_importance <- paired_permutation_losses(tcga_pairs, seed = 14521L)
importance <- rbind(
  cbind(cohort = "GSE14520 external", geo_importance),
  cbind(cohort = "TCGA-LIHC apparent training", tcga_importance)
)
qc <- data.frame(
  measure = c("TCGA sourced samples", "TCGA recurrent tumors excluded", "TCGA matched pairs",
              "TCGA measured genes", "GEO GPL3921 samples", "GEO measured probes",
              "GEO uniquely mapped probe rows used", "Shared reference genes",
              "GEO tumor-only patients", "GEO non-tumor-only patients",
              "GEO complete pairs"),
  count = c(length(tcga_ids), sum(tcga_type == "02"), nrow(tcga_pairs$pos),
            nrow(tcga), nrow(geo_meta), nrow(geo), nrow(mapping), length(reference),
            sum(geo_data$tumor) - n, sum(geo_data$tumor == 0L) - n, n)
)
geo_data$paired_primary <- geo_data$patient %in% geo_pairs$pos$patient
geo_data$score_without_marker_x_at <- geo_marker_sensitivity$score
mapping$role <- ifelse(mapping$gene %in% genes, "marker",
                       ifelse(mapping$gene %in% reference, "reference", "unused"))
write_tsv <- function(table, filename) {
  write.table(table, file.path(result_dir, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
write_tsv(metrics, "performance.tsv")
write_tsv(importance, "permutation_importance.tsv")
write_tsv(qc, "sample_and_mapping_qc.tsv")
write_tsv(data.frame(term = names(coef(fit)), coefficient = unname(coef(fit))),
          "coefficients.tsv")
write_tsv(geo_data, "GSE14520_predictions.tsv")
write_tsv(roc_curve, "GSE14520_paired_ROC_coordinates.tsv")
write_tsv(roc(specific_pairs$neg$score, specific_pairs$pos$score),
          "GSE14520_paired_ROC_without_marker_x_at.tsv")
write_tsv(mapping[mapping$role != "unused", ], "probe_mapping.tsv")
write_tsv(data.frame(gene = reference), "reference_genes.tsv")
saveRDS(
  list(coefficient = coef(fit), features = genes, reference_genes = reference,
       gpl3921_mapping = mapping[mapping$role != "unused", ],
       model = "binomial logistic with within-sample midrank percentiles"),
  file.path(result_dir, "locked_TCGA_model.rds")
)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "R_session.txt"))
print(qc, row.names = FALSE)
print(metrics, row.names = FALSE, digits = 6)
print(importance, row.names = FALSE, digits = 6)
