align_pairs <- function(data) {
  pos <- data[data$tumor == 1L, ]
  neg <- data[data$tumor == 0L, ]
  patients <- sort(intersect(pos$patient, neg$patient), method = "radix")
  pairs <- list(pos = pos[match(patients, pos$patient), ], neg = neg[match(patients, neg$patient), ])
  stopifnot(length(patients) > 1L, identical(pairs$pos$patient, pairs$neg$patient))
  pairs
}
auc <- function(negative, positive) {
  n <- length(negative)
  m <- length(positive)
  stopifnot(n > 0L, m > 0L)
  ranks <- rank(c(negative, positive), ties.method = "average")
  (sum(ranks[n + seq_len(m)]) - m * (m + 1L) / 2) / (m * n)
}
roc <- function(negative, positive) {
  y <- c(rep(0L, length(negative)), rep(1L, length(positive)))
  scores <- c(negative, positive)
  order_idx <- order(scores, decreasing = TRUE)
  scores <- scores[order_idx]
  y <- y[order_idx]
  block_end <- which(!duplicated(scores, fromLast = TRUE))
  data.frame(
    threshold = c(Inf, scores[block_end]),
    fpr = c(0, cumsum(y == 0L)[block_end] / length(negative)),
    tpr = c(0, cumsum(y == 1L)[block_end] / length(positive))
  )
}
ci95 <- function(x) as.numeric(quantile(x, c(0.025, 0.975), names = FALSE))
