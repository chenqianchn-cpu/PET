# Run from the GSE14520 external-validation directory:
# Rscript --vanilla model_comparison.R
suppressPackageStartupMessages(library(caret))
source("prepare_expression.R")
source("validation_metrics.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L, all(args %in% c("proxy15", "strict14")))
proxy <- !length(args) || args == "proxy15"
output_dir <- if (proxy) "ML15/results" else "ML15/results_strict14"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
panel15 <- c("ANG", "ATIC", "CDC25B", "EPHX2", "FKBP1A", "GNPDA1", "HRAS", "IMPDH2",
           "LGALS3", "MIF", "MMP9", "NQO1", "PPIA", "RAP2A", "SDS")
panel <- if (proxy) panel15 else setdiff(panel15, "ANG")
# The official platform annotation assigns ANG's only probe to RNASE4 and ANG.
# Retain the probe solely as a declared proxy in the requested 15-feature run.
proxy_probe <- "205141_at"
stopifnot(proxy_probe %in% rownames(geo),
          any(annotation$ID == proxy_probe &
                annotation[["Gene symbol"]] == "RNASE4///ANG", na.rm = TRUE),
          !("ANG" %in% rownames(geo_by_gene)))
if (proxy) {
  geo_by_gene <- rbind(geo_by_gene, ANG = geo[proxy_probe, ])
}
complete_tcga <- rowSums(!is.finite(tcga)) == 0L
reference <- sort(setdiff(
  intersect(rownames(tcga)[complete_tcga], rownames(geo_by_gene)), c(panel15, "RNASE4")
), method = "radix")
stopifnot(all(panel %in% rownames(tcga)[complete_tcga]), all(panel %in% rownames(geo_by_gene)),
          length(reference) > 1000L)

panel_percentiles <- function(expression) {
  ref <- expression[reference, , drop = FALSE]
  marks <- expression[panel, , drop = FALSE]
  stopifnot(all(is.finite(ref)), all(is.finite(marks)))
  answer <- t(vapply(seq_len(ncol(expression)), function(j) {
    sorted_ref <- sort(ref[, j], method = "radix")
    (findInterval(marks[, j], sorted_ref, left.open = TRUE) +
       findInterval(marks[, j], sorted_ref)) / (2 * length(sorted_ref))
  }, numeric(length(panel))))
  dimnames(answer) <- list(colnames(expression), panel)
  as.data.frame(answer)
}
x_tcga <- panel_percentiles(tcga)
x_geo <- panel_percentiles(geo_by_gene)
ppia <- mapping$probe[mapping$gene == "PPIA"]
ppia_specific <- ppia[!endsWith(ppia, "_x_at")]
stopifnot(length(ppia_specific) > 0L, length(ppia_specific) < length(ppia))
geo_specific <- geo_by_gene
geo_specific["PPIA", ] <- apply(geo[ppia_specific, , drop = FALSE], 2L, median)
x_geo_specific <- panel_percentiles(geo_specific)
stopifnot(identical(x_geo[setdiff(panel, "PPIA")], x_geo_specific[setdiff(panel, "PPIA")]),
          identical(rownames(x_tcga), tcga_meta$sample), identical(rownames(x_geo), geo_meta$sample))
tcga_groups <- align_pairs(tcga_meta)
geo_groups <- align_pairs(geo_meta)

# The historical exports give algorithm labels but no method IDs or grids.
# LogitBoost is the declared boosting alternative for the ambiguous historical Logit label.
specs <- list(
  PLS = list(method = "pls", grid = data.frame(ncomp = c(2L, 4L, 6L))),
  NaiveBayes = list(method = "nb", grid = expand.grid(fL = 0, usekernel = c(FALSE, TRUE), adjust = 1)),
  GBM = list(method = "gbm", grid = expand.grid(n.trees = c(100L, 200L),
              interaction.depth = c(1L, 2L), shrinkage = 0.05, n.minobsinnode = 5L),
              extra = list(verbose = FALSE)),
  GLM = list(method = "glm", grid = data.frame(parameter = "none"), extra = list(family = binomial())),
  RF = list(method = "rf", grid = data.frame(mtry = c(3L, 5L, 8L)), extra = list(ntree = 300L)),
  SVM = list(method = "svmRadial", grid = expand.grid(sigma = c(0.04, 0.12), C = c(0.5, 2, 8))),
  stepLDA = list(method = "stepLDA", grid = data.frame(maxvar = c(4L, 8L, 12L), direction = "forward")),
  ElasticNet = list(method = "glmnet", grid = expand.grid(alpha = c(0.25, 0.5, 0.75),
                                                        lambda = c(0.005, 0.02, 0.08))),
  KNN = list(method = "knn", grid = data.frame(k = c(5L, 11L, 21L))),
  LogitBoost = list(method = "LogitBoost", grid = data.frame(nIter = c(25L, 50L, 100L)))
)
for (name in names(specs)) {
  spec <- specs[[name]]
  model <- getModelInfo(spec$method, regex = FALSE)[[1L]]
  stopifnot(!is.null(model), setequal(names(spec$grid), model$parameters$parameter))
  missing <- model$library[!vapply(model$library, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) stop(name, " requires packages: ", paste(missing, collapse = ", "))
}

# A patient is always wholly in a single validation fold. Stratify by the
# type of sample(s) available for that patient (paired, tumor-only, or normal-only).
patient_folds <- function(meta, k, seed) {
  groups <- split(seq_len(nrow(meta)), meta$patient)
  profiles <- vapply(groups, function(i) paste(sort(meta$tumor[i]), collapse = ""), character(1))
  set.seed(seed)
  fold <- integer(nrow(meta))
  for (indices in split(seq_along(groups), profiles)) {
    shuffled <- indices[sample.int(length(indices))]
    for (i in seq_along(shuffled)) {
      fold[groups[[shuffled[i]]]] <- (i - 1L) %% k + 1L
    }
  }
  stopifnot(all(table(fold, meta$tumor) > 0),
            all(vapply(groups, function(i) length(unique(fold[i])) == 1L, logical(1))))
  fold
}
fit_model <- function(name, x, meta, folds, seed) {
  spec <- specs[[name]]
  training_indices <- lapply(sort(unique(folds)), function(i) which(folds != i))
  testing_indices <- lapply(sort(unique(folds)), function(i) which(folds == i))
  control <- trainControl(
    method = "cv", number = length(testing_indices), classProbs = TRUE,
    summaryFunction = twoClassSummary, index = training_indices, indexOut = testing_indices,
    allowParallel = FALSE, returnResamp = "final"
  )
  target <- factor(ifelse(meta$tumor == 1L, "Tumor", "NonTumor"),
                   levels = c("Tumor", "NonTumor"))
  set.seed(seed)
  fit <- do.call(caret::train, c(list(
    x = x, y = target, method = spec$method, tuneGrid = spec$grid,
    preProcess = c("center", "scale"), metric = "ROC", trControl = control
  ), spec$extra))
  stopifnot(nrow(fit$bestTune) == 1L)
  fit
}

# klaR normalizes exp(log-likelihoods), which can round both classes to zero
# for a shifted external sample. Use the same fitted densities and threshold
# and normalize on the log scale for affected rows.
stable_nb_score <- function(fit, data) {
  nb <- fit$finalModel
  x <- predict(fit$preProcess, data)[nb$varnames]
  classes <- nb$levels
  tumor <- match("Tumor", classes)
  stopifnot(!is.na(tumor), all(vapply(x, is.numeric, logical(1))))
  kernel_density <- getFromNamespace("dkernel", "klaR")
  vapply(seq_len(nrow(x)), function(i) {
    log_density <- vapply(seq_along(x), function(j) {
      distribution <- nb$tables[[j]]
      value <- x[[j]][i]
      density <- if (nb$usekernel) {
        vapply(distribution, function(d) kernel_density(x = value, kernel = d), numeric(1))
      } else {
        dnorm(value, distribution[, 1L], distribution[, 2L])
      }
      density[density == 0] <- 0.001
      log(density)
    }, numeric(length(classes)))
    joint <- log(as.numeric(nb$apriori)) + rowSums(log_density)
    scaled <- exp(joint - max(joint))
    scaled[tumor] / sum(scaled)
  }, numeric(1))
}

score_model <- function(fit, data) {
  score <- unname(predict(fit, newdata = data, type = "prob")$Tumor)
  invalid <- !is.finite(score)
  if (fit$method == "nb" && any(invalid)) {
    score[invalid] <- stable_nb_score(fit, data[invalid, , drop = FALSE])
  }

  stopifnot(length(score) == nrow(data), all(is.finite(score)), all(score >= 0), all(score <= 1))
  score
}
# Every complete patient pair stays together. For each gene, predict both tissues
# after exchanging just that gene's values. This allows exact reuse of predictions
# across subsequent random swap patterns, without refitting or predicting again.
score_and_swap <- function(fit, data, meta) {
  pairs <- align_pairs(meta)
  pos <- match(pairs$pos$sample, meta$sample)
  neg <- match(pairs$neg$sample, meta$sample)
  original <- score_model(fit, data)
  paired <- data[c(pos, neg), , drop = FALSE]
  n <- length(pos)
  both_swapped <- vapply(panel, function(gene) {
    changed <- paired
    changed[[gene]] <- paired[[gene]][c(n + seq_len(n), seq_len(n))]
    score_model(fit, changed)
  }, numeric(2L * n))
  list(score = original, pos = pos, neg = neg,
       swapped_pos = both_swapped[seq_len(n), , drop = FALSE],
       swapped_neg = both_swapped[n + seq_len(n), , drop = FALSE])
}
tune_record <- function(fit, cohort, name, stage) {
  vals <- paste(names(fit$bestTune), unlist(fit$bestTune), sep = "=", collapse = ", ")
  data.frame(cohort = cohort, model = name, fit = stage, method = fit$method, parameters = vals)
}
crossfit <- function(name, x, meta, outer, seed, cohort) {
  pairs <- align_pairs(meta)
  n <- nrow(pairs$pos)
  result <- list(score = rep(NA_real_, nrow(meta)), pos = match(pairs$pos$sample, meta$sample),
                 neg = match(pairs$neg$sample, meta$sample),
                 swapped_pos = matrix(NA_real_, n, length(panel), dimnames = list(NULL, panel)),
                 swapped_neg = matrix(NA_real_, n, length(panel), dimnames = list(NULL, panel)))
  tune <- vector("list", length(unique(outer)))
  for (i in sort(unique(outer))) {
    held <- which(outer == i)
    train_meta <- meta[-held, , drop = FALSE]
    inner <- patient_folds(train_meta, 3L, seed + i)
    fit <- fit_model(name, x[-held, , drop = FALSE], train_meta, inner, seed + 100L + i)
    predicted <- score_and_swap(fit, x[held, , drop = FALSE], meta[held, , drop = FALSE])
    result$score[held] <- predicted$score
    pair_rows <- match(held[predicted$pos], result$pos)
    stopifnot(!anyNA(pair_rows), identical(held[predicted$neg], result$neg[pair_rows]))
    result$swapped_pos[pair_rows, ] <- predicted$swapped_pos
    result$swapped_neg[pair_rows, ] <- predicted$swapped_neg
    tune[[i]] <- tune_record(fit, cohort, name, paste0("outer_fold_", i))
  }
  stopifnot(all(is.finite(result$score)), all(is.finite(result$swapped_pos)),
            all(is.finite(result$swapped_neg)))
  list(predicted = result, tune = do.call(rbind, tune))
}

assess <- function(scored, meta, model, analysis, seed) {
  pos <- scored$score[scored$pos]
  neg <- scored$score[scored$neg]
  n <- length(pos)
  set.seed(seed)
  boot <- replicate(2000L, sample.int(n, n, replace = TRUE))
  boot_auc <- apply(boot, 2L, function(i) auc(neg[i], pos[i]))
  bounds <- ci95(boot_auc)
  auc_all <- auc(scored$score[meta$tumor == 0L], scored$score[meta$tumor == 1L])
  paired_auc <- auc(neg, pos)
  curve <- roc(neg, pos)
  curve_auc <- sum(diff(curve$fpr) * (head(curve$tpr, -1L) + tail(curve$tpr, -1L)) / 2)
  stopifnot(isTRUE(all.equal(curve_auc, paired_auc, tolerance = 1e-12)))
  observed_error <- (1 - pos)^2 + neg^2
  change <- (1 - scored$swapped_pos)^2 + scored$swapped_neg^2 - observed_error
  n_draws <- 1024L
  masks <- matrix(sample(c(FALSE, TRUE), n * n_draws, replace = TRUE), nrow = n)
  baseline <- sqrt(mean(observed_error) / 2)
  permuted_mse <- (sum(observed_error) + crossprod(masks, change)) / (2 * n)
  stopifnot(is.matrix(permuted_mse), all(permuted_mse >= -1e-12))
  delta <- colMeans(sqrt(pmax(permuted_mse, 0))) - baseline
  ranks <- rank(-delta, ties.method = "average")
  list(
    performance = data.frame(analysis = analysis, model = model, pairs = n, samples = nrow(meta),
        tumor = sum(meta$tumor), non_tumor = sum(meta$tumor == 0L),
        AUROC_paired = paired_auc, CI_low = bounds[1L], CI_high = bounds[2L],
        AUROC_all_samples = auc_all, RMSE_paired = baseline),
    ranks = data.frame(analysis = analysis, model = model, gene = panel, delta_RMSE = delta, rank = ranks),
    scores = data.frame(analysis = analysis, model = model, sample = meta$sample,
                        patient = meta$patient, tumor = meta$tumor, tumor_probability = scored$score),
    curves = cbind(analysis = analysis, model = model, curve)
  )
}
write_tsv <- function(data, filename) {
  write.table(data, file.path(output_dir, filename), sep = "\t", row.names = FALSE, quote = FALSE)
}

outer_tcga <- patient_folds(tcga_meta, 5L, 1452015L)
outer_geo <- patient_folds(geo_meta, 5L, 1452016L)
full_folds <- patient_folds(tcga_meta, 5L, 1452017L)
make_folds <- function(meta, fold, cohort) {
  data.frame(cohort = cohort, meta[c("sample", "patient", "tumor")], fold = fold)
}
write_tsv(rbind(make_folds(tcga_meta, outer_tcga, "TCGA_outer"),
                make_folds(geo_meta, outer_geo, "GSE14520_outer"),
                make_folds(tcga_meta, full_folds, "TCGA_final_tuning")), "patient_folds.tsv")
write_tsv(data.frame(gene = panel), "gene_panel.tsv")
write_tsv(data.frame(gene = reference), "reference_genes.tsv")
gene_probes <- mapping[mapping$gene %in% panel, ]
gene_probes$mapping_type <- "unique symbol"
if (proxy) {
  gene_probes <- rbind(gene_probes, data.frame(
    probe = proxy_probe, gene = "ANG", mapping_type = "shared ANG/RNASE4 proxy"
  ))
}
gene_probes$retained_in_PPIA_sensitivity <- (gene_probes$gene != "PPIA" |
                                              gene_probes$probe %in% ppia_specific)
write_tsv(gene_probes, "GSE14520_panel_probe_mapping.tsv")

results <- setNames(vector("list", length(specs)), names(specs))
tuning <- results
failures <- list()
warning_log <- list()
model_dir <- file.path(output_dir, "locked_TCGA_models")
dir.create(model_dir, showWarnings = FALSE)
checkpoint <- Sys.getenv("PET_ML15_CHECKPOINT")
if (nzchar(checkpoint)) {
  checkpoint <- file.path(checkpoint, if (proxy) "proxy15" else "strict14")
  dir.create(checkpoint, recursive = TRUE, showWarnings = FALSE)
}
for (j in seq_along(specs)) {
  name <- names(specs)[j]
  seed <- 781000L + 1000L * j
  for (stage in c("TCGA", "GSE14520", "transfer")) {
    warnings_found <- character()
    cache_file <- file.path(checkpoint, paste0(name, "-", stage, ".rds"))
    cached <- nzchar(checkpoint) && file.exists(cache_file)
    outcome <- if (cached) readRDS(cache_file) else tryCatch(withCallingHandlers({
      if (stage == "TCGA") {
        tc <- crossfit(name, x_tcga, tcga_meta, outer_tcga, seed, stage)
        list(metrics = list(assess(tc$predicted, tcga_meta, name, "TCGA grouped CV", 5091L)),
             tune = tc$tune)
      } else if (stage == "GSE14520") {
        gs <- crossfit(name, x_geo, geo_meta, outer_geo, seed + 300L, stage)
        list(metrics = list(assess(gs$predicted, geo_meta, name, "GSE14520 grouped CV", 5092L)),
             tune = gs$tune)
      } else {
        frozen <- fit_model(name, x_tcga, tcga_meta, full_folds, seed + 600L)
        transferred <- score_and_swap(frozen, x_geo, geo_meta)
        specific <- score_and_swap(frozen, x_geo_specific, geo_meta)
        frozen$trainingData <- NULL
        saveRDS(frozen, file.path(model_dir, paste0(name, ".rds")))
        list(metrics = list(assess(transferred, geo_meta, name, "TCGA to GSE14520", 5092L),
                            assess(specific, geo_meta, name, "TCGA to GSE14520 without PPIA x_at", 5092L)),
             tune = tune_record(frozen, "TCGA", name, "full_training"))
      }
    }, warning = function(w) {
      warnings_found <<- c(warnings_found, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    if (!cached && !inherits(outcome, "error") && nzchar(checkpoint)) saveRDS(outcome, cache_file)
    if (length(warnings_found)) {
      distinct <- table(warnings_found)
      warning_log[[paste(name, stage)]] <- data.frame(model = name, stage = stage,
        count = as.integer(distinct), message = gsub("[\r\n\t]+", " ", names(distinct)))
    }
    if (inherits(outcome, "error")) {
      failures[[paste(name, stage)]] <- data.frame(model = name, stage = stage,
                                                  error = conditionMessage(outcome))
      message(name, " / ", stage, " failed: ", conditionMessage(outcome))
    } else {
      results[[name]] <- c(results[[name]], outcome$metrics)
      tuning[[name]] <- rbind(tuning[[name]], outcome$tune)
      message(name, " / ", stage, " complete")
    }
  }
}
if (length(warning_log)) write_tsv(do.call(rbind, warning_log), "model_warnings.tsv")
if (length(failures)) write_tsv(do.call(rbind, failures), "model_failures.tsv")
all_results <- unlist(results, recursive = FALSE)
if (!length(all_results)) stop("No model assessment completed.")
extract <- function(x) do.call(rbind, lapply(all_results, "[[", x))
performance <- extract("performance")
ranks <- extract("ranks")
groups <- split(ranks, list(ranks$analysis, ranks$gene), drop = TRUE)
rank_summary <- do.call(rbind, lapply(groups, function(g) data.frame(
  analysis = g$analysis[1L], gene = g$gene[1L], models = nrow(g),
  mean_rank = mean(g$rank), median_rank = median(g$rank),
  best_rank = min(g$rank), worst_rank = max(g$rank),
  top_5 = sum(g$rank <= 5), top_10 = sum(g$rank <= 10), positive_loss = sum(g$delta_RMSE > 0)
)))
rank_summary <- rank_summary[order(rank_summary$analysis, rank_summary$mean_rank, rank_summary$gene), ]
write_tsv(performance, "model_performance.tsv")
write_tsv(ranks, "feature_ranks_by_model.tsv")
write_tsv(rank_summary, "gene_rank_summary.tsv")
write_tsv(extract("scores"), "out_of_fold_and_external_scores.tsv")
write_tsv(extract("curves"), "paired_ROC_coordinates.tsv")
write_tsv(do.call(rbind, tuning), "selected_hyperparameters.tsv")
grids <- do.call(rbind, lapply(names(specs), function(name) {
  grid <- specs[[name]]$grid
  values <- apply(grid, 1L, function(row) paste(names(grid), row, sep = "=", collapse = ", "))
  data.frame(model = name, method = specs[[name]]$method, parameters = values)
}))
write_tsv(grids, "tuning_grids.tsv")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_session.txt"))
print(performance[, c("analysis", "model", "AUROC_paired", "CI_low", "CI_high")], row.names = FALSE)
print(rank_summary[rank_summary$gene %in% c("PPIA", "IMPDH2"), ], row.names = FALSE)
if (length(failures)) quit(status = 1L)
