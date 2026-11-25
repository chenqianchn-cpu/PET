# Run from this directory with Rscript --vanilla candidate_ridge_validation.R.
# Requires glmnet. The 334 input symbols are the case-sensitive intersection
# of the PET union (435 symbols, workbook tab Sheet1 column G) and
# HCC union (9495 symbols, workbook tab Sheet1 column E).
source("prepare_expression.R")
source("validation_metrics.R")
stopifnot(requireNamespace("glmnet", quietly = TRUE))
candidate_symbols <- read.delim(file.path(source_dir, "candidate_genes.tsv"))$gene
stopifnot(length(candidate_symbols) == 334L, !anyDuplicated(candidate_symbols),
          all(genes %in% candidate_symbols))

# Candidate availability is determined by the official array annotation and
# TCGA expression completeness. The reference gene set excludes all candidates.
complete_tcga <- rowSums(!is.finite(tcga)) == 0L
reference <- sort(setdiff(
  intersect(rownames(tcga)[complete_tcga], rownames(geo_by_gene)), candidate_symbols
), method = "radix")
available <- intersect(candidate_symbols, intersect(rownames(tcga)[complete_tcga], rownames(geo_by_gene)))
stopifnot(all(genes %in% available), length(reference) > 1000L)

relative_expression <- function(values, features) {
  reference_values <- values[reference, , drop = FALSE]
  marker_values <- values[features, , drop = FALSE]
  stopifnot(all(is.finite(reference_values)), all(is.finite(marker_values)))
  result <- t(vapply(seq_len(ncol(values)), function(j) {
    ref <- sort(reference_values[, j], method = "radix")
    marker <- marker_values[, j]
    lower <- findInterval(marker, ref, left.open = TRUE)
    upper <- findInterval(marker, ref)
    (lower + upper) / (2 * length(ref))
  }, numeric(length(features))))
  dimnames(result) <- list(colnames(values), features)
  result
}
x_tcga <- relative_expression(tcga, available)
variable_tcga <- vapply(as.data.frame(x_tcga), function(x) diff(range(x)) > 0, logical(1))
features <- available[variable_tcga]
stopifnot(all(genes %in% features))
x_tcga <- x_tcga[, features, drop = FALSE]
x_geo <- relative_expression(geo_by_gene, features)
specific_map <- mapping[mapping$gene %in% genes & !endsWith(mapping$probe, "_x_at"), ]
specific_gene_values <- geo_by_gene
for (gene in genes) {
  probes <- specific_map$probe[specific_map$gene == gene]
  stopifnot(length(probes) > 0L)
  specific_gene_values[gene, ] <- apply(geo[probes, , drop = FALSE], 2L, median)
}
x_geo_specific <- x_geo
x_geo_specific[, genes] <- relative_expression(specific_gene_values, genes)

# Fix alpha=0 (ridge) in advance. Use patient-grouped stratified CV only in TCGA
# and select lambda using the one-standard-error rule for binomial deviance.
set.seed(14520334L)
groups <- split(seq_len(nrow(tcga_meta)), tcga_meta$patient)
group_profiles <- vapply(groups, function(i) {
  paste(sort(tcga_meta$tumor[i]), collapse = "")
}, character(1))
folds <- 10L
group_fold <- integer(length(groups))
for (profile in unique(group_profiles)) {
  group_index <- which(group_profiles == profile)
  group_fold[sample(group_index)] <- rep(seq_len(folds), length.out = length(group_index))
}
fold_id <- integer(nrow(tcga_meta))
for (j in seq_along(groups)) {
  fold_id[groups[[j]]] <- group_fold[j]
}
stopifnot(length(unique(fold_id)) == folds, all(table(fold_id, tcga_meta$tumor) > 0),
          all(vapply(groups, function(idx) length(unique(fold_id[idx])) == 1L, logical(1))))
fit <- glmnet::cv.glmnet(
  x = x_tcga, y = tcga_meta$tumor, family = "binomial", alpha = 0,
  foldid = fold_id, type.measure = "deviance", standardize = TRUE
)
fit_weights <- as.matrix(coef(fit, s = "lambda.1se"))[, 1L]
intercept <- unname(fit_weights["(Intercept)"])
beta <- fit_weights[features]
stopifnot(all(is.finite(beta)), is.finite(intercept))
geo_meta$score <- drop(intercept + x_geo %*% beta)
geo_meta$score_without_marker_x_at <- drop(intercept + x_geo_specific %*% beta)
stopifnot(isTRUE(all.equal(geo_meta$score, as.vector(predict(fit, x_geo, s = "lambda.1se")))))
patients <- align_pairs(geo_meta)
pos <- match(patients$pos$sample, geo_meta$sample)
neg <- match(patients$neg$sample, geo_meta$sample)
n <- length(pos)

# Expected AUROC after a randomly selected single-gene swap within each pair.
# For different patients the two swaps are independent. For the same patient
# they share a draw, so correct the diagonal of the independent-swap AUROC.
set.seed(14520335L)
boot_n <- 2000L
boot_idx <- cbind(seq_len(n), replicate(boot_n, sample.int(n, n, replace = TRUE)))
boot_weights <- vapply(seq_len(ncol(boot_idx)), function(j) {
  tabulate(boot_idx[, j], nbins = n)
}, integer(n))
interval <- function(x) as.numeric(quantile(x, c(0.025, 0.975), names = FALSE))
win <- function(positive, negative) (positive > negative) + 0.5 * (positive == negative)

weighted_aucs <- function(negative, positive, negative_ids, positive_ids) {
  index <- order(negative, method = "radix")
  sorted_negative <- negative[index]
  negative_patients <- negative_ids[index]
  less <- findInterval(positive, sorted_negative, left.open = TRUE) + 1L
  leq <- findInterval(positive, sorted_negative) + 1L
  denominator <- 2 * length(negative) * length(positive)
  vapply(seq_len(ncol(boot_weights)), function(j) {
    weights <- boot_weights[, j]
    cumulative <- c(0L, cumsum(weights[negative_patients]))
    sum(weights[positive_ids] * (cumulative[less] + cumulative[leq])) / denominator
  }, numeric(1))
}
paired_importance <- function(features_matrix, scores) {
  score_pos <- scores[pos]
  score_neg <- scores[neg]
  patients <- seq_len(n)
  both <- rep(patients, 2L)
  observed <- weighted_aucs(score_neg, score_pos, patients, patients)
  stopifnot(isTRUE(all.equal(observed[1L], auc(score_neg, score_pos), tolerance = 1e-12)))
  contributions <- sweep(features_matrix[pos, , drop = FALSE] -
                           features_matrix[neg, , drop = FALSE], 2L, beta, "*")
  keep_both <- win(score_pos, score_neg)
  vapply(seq_along(features), function(j) {
    shift <- contributions[, j]
    swap_pos <- score_pos - shift
    swap_neg <- score_neg + shift
    pool_expected <- weighted_aucs(c(score_neg, swap_neg), c(score_pos, swap_pos), both, both)
    diagonal <- (keep_both + win(swap_pos, swap_neg) -
                   win(score_pos, swap_neg) - win(swap_pos, score_neg)) / 4
    correction <- drop(crossprod(diagonal, boot_weights)) / n^2
    round(observed - (pool_expected + correction), digits = 13L)
  }, numeric(ncol(boot_weights)))
}
tabulate_importance <- function(features_matrix, scores, probe_policy) {
  losses <- paired_importance(features_matrix, scores)
  point <- losses[1L, ]
  boot <- t(losses[-1L, , drop = FALSE])
  # Tied genes share the smallest rank that they occupy. Report ties explicitly.
  bootstrap_ranks <- apply(boot, 2L, function(x) rank(-x, ties.method = "min"))
  loss_ci <- t(apply(boot, 1L, interval))
  rank_ci <- t(apply(bootstrap_ranks, 1L, function(x) {
    as.integer(quantile(x, c(0.025, 0.975), names = FALSE, type = 1L))
  }))
  output <- data.frame(
    probe_policy = probe_policy, gene = features, delta_auc = point,
    ci_low = loss_ci[, 1L], ci_high = loss_ci[, 2L],
    rank = rank(-point, ties.method = "min"),
    genes_tied_at_point_rank = vapply(point, function(x) sum(point == x), integer(1)),
    rank_ci_low = rank_ci[, 1L], rank_ci_high = rank_ci[, 2L],
    bootstrap_fraction_rank_le_10 = rowMeans(bootstrap_ranks <= 10L),
    pairs = n, candidate_features = length(features)
  )
  output[order(output$rank, output$gene), ]
}
importance <- rbind(
  tabulate_importance(x_geo, geo_meta$score, "All uniquely annotated marker probes"),
  tabulate_importance(x_geo_specific, geo_meta$score_without_marker_x_at,
                      "Omit marker probes ending in _x_at")
)

performance <- do.call(rbind, lapply(c("score", "score_without_marker_x_at"), function(column) {
  score <- geo_meta[[column]]
  values <- apply(boot_idx[, -1L, drop = FALSE], 2L, function(idx) auc(score[neg[idx]], score[pos[idx]]))
  range <- interval(values)
  data.frame(score_definition = column, auc = auc(score[neg], score[pos]),
             ci_low = range[1L], ci_high = range[2L], pairs = n,
             reference_genes = length(reference), candidate_input = length(candidate_symbols),
             candidate_features = length(features), lambda = fit$lambda.1se,
             cv_folds = folds, bootstrap_replicates = boot_n)
}))
write_tsv <- function(tab, filename) {
  write.table(tab, file.path(result_dir, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
write_tsv(performance, "candidate_ridge_performance.tsv")
write_tsv(importance, "candidate_ridge_GSE14520_feature_ranks.tsv")
write_tsv(data.frame(gene = candidate_symbols, included = candidate_symbols %in% features,
                     in_TCGA = candidate_symbols %in% rownames(tcga),
                     TCGA_complete = candidate_symbols %in% rownames(tcga)[complete_tcga],
                     on_GPL3921 = candidate_symbols %in% rownames(geo_by_gene)),
          "candidate_ridge_gene_qc.tsv")
write_tsv(data.frame(sample = geo_meta$sample, patient = geo_meta$patient,
                     tumor = geo_meta$tumor, paired = seq_len(nrow(geo_meta)) %in% c(pos, neg),
                     score = geo_meta$score,
                     score_without_marker_x_at = geo_meta$score_without_marker_x_at),
          "candidate_ridge_GSE14520_predictions.tsv")
write_tsv(data.frame(patient = tcga_meta$patient, sample = tcga_meta$sample, fold = fold_id),
          "candidate_ridge_TCGA_cv_folds.tsv")
saveRDS(list(coefficient = fit_weights, features = features, reference_genes = reference,
             gpl3921_mapping = mapping[mapping$gene %in% c(features, reference), ],
             lambda = fit$lambda.1se, cv_folds = folds,
             glmnet_version = as.character(utils::packageVersion("glmnet"))),
        file.path(result_dir, "locked_candidate_ridge_TCGA_model.rds"))
print(performance, row.names = FALSE, digits = 6)
print(importance[importance$gene %in% genes, ], row.names = FALSE, digits = 6)
