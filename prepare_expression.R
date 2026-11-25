source_dir <- "source"
result_dir <- "results"
dir.create(result_dir, showWarnings = FALSE)
genes <- c("PPIA", "IMPDH2")

source_names <- c(
  "TCGA-LIHC_Xena_HiSeqV2.tsv.gz",
  "GSE14520-GPL3921_series_matrix.txt.gz",
  "GPL3921.annot.gz"
)
input <- file.path(source_dir, source_names)
stopifnot(all(file.info(input)$size == c(22166005L, 21301282L, 4551278L)))

read_gzip_lines <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  readLines(con, warn = TRUE, encoding = "UTF-8")
}

table_bounds <- function(lines, markers) {
  indices <- match(markers, lines)
  stopifnot(all(is.finite(indices)), indices[1] + 1L < indices[2])
  seq.int(indices[1] + 1L, indices[2] - 1L)
}

expression_matrix <- function(tab) {
  stopifnot(!anyDuplicated(tab[[1L]]), !anyDuplicated(names(tab)[-1L]))
  values <- as.matrix(tab[-1L])
  stopifnot(is.numeric(values))
  dimnames(values) <- list(tab[[1L]], names(tab)[-1L])
  values
}

# Cohort labels come from the source metadata and must agree with independent tags.
tcga <- expression_matrix(read.delim(
  gzfile(input[1L]), check.names = FALSE, quote = "", comment.char = ""
))
tcga_ids <- colnames(tcga)
stopifnot(all(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$", tcga_ids)))
tcga_type <- substr(tcga_ids, 14L, 15L)
stopifnot(all(tcga_type %in% c("01", "02", "11")))
tcga_keep <- tcga_type %in% c("01", "11")
tcga <- tcga[, tcga_keep, drop = FALSE]
tcga_meta <- data.frame(
  sample = tcga_ids[tcga_keep], patient = substr(tcga_ids[tcga_keep], 1L, 12L),
  tumor = as.integer(tcga_type[tcga_keep] == "01")
)
stopifnot(!anyDuplicated(tcga_meta[c("patient", "tumor")]))

geo_lines <- read_gzip_lines(input[2L])
geo_bounds <- table_bounds(geo_lines, c("!series_matrix_table_begin", "!series_matrix_table_end"))
geo_meta_lines <- strsplit(geo_lines[seq_len(min(geo_bounds) - 1L)], "\t", fixed = TRUE)
geo_field <- function(key, prefix = NULL) {
  hits <- Filter(function(row) row[1L] == key, geo_meta_lines)
  if (!is.null(prefix)) {
    hits <- Filter(function(row) startsWith(tolower(row[2L]), paste0('"', prefix)), hits)
  }
  stopifnot(length(hits) == 1L)
  sub('^"(.*)"$', "\\1", hits[[1L]][-1L])
}
geo_ids <- geo_field("!Sample_geo_accession")
sample_title <- geo_field("!Sample_title")
source_id <- geo_field("!Sample_source_name_ch1")
tissue <- tolower(geo_field("!Sample_characteristics_ch1", "tissue:"))
disease <- tolower(geo_field("!Sample_characteristics_ch1", "disease state:"))
individual <- sub("^[^:]+: ", "", geo_field("!Sample_characteristics_ch1", "individual:"))
allowed_tissue <- c("tissue: liver non-tumor tissue", "tissue: liver tumor tissue")
stopifnot(
  all(tissue %in% allowed_tissue),
  all(disease == "disease state: hepatocellular carcinoma (hcc)"),
  all(grepl("^(LCS-[0-9]+|X[0-9]+-[0-9]+)[AB]$", source_id)),
  !anyDuplicated(geo_ids), !anyDuplicated(source_id)
)
geo_tumor <- as.integer(tissue == allowed_tissue[2L])
stopifnot(
  all(endsWith(source_id, ifelse(geo_tumor == 1L, "A", "B"))),
  all(endsWith(sample_title, paste0(" ", source_id))),
  all(startsWith(
    sample_title, ifelse(geo_tumor == 1L, "Liver Tumor Tissue ", "Liver Non-Tumor Tissue ")
  ))
)
geo_meta <- data.frame(
  sample = geo_ids, patient = sub("[AB]$", "", source_id), tumor = geo_tumor,
  source_id = source_id, individual = individual
)
stopifnot(!anyDuplicated(geo_meta[c("patient", "tumor")]))
geo <- expression_matrix(read.delim(
  text = geo_lines[geo_bounds], check.names = FALSE, comment.char = ""
))
stopifnot(identical(colnames(geo), geo_meta$sample))

# Use the official, fixed probe annotation. Exclude probes with ambiguous symbols.
# Collapse remaining probes within each sample by their median, independent of outcome.
annotation_lines <- read_gzip_lines(input[3L])
annotation_bounds <- table_bounds(annotation_lines, c("!platform_table_begin", "!platform_table_end"))
annotation <- read.delim(
  text = annotation_lines[annotation_bounds], check.names = FALSE,
  comment.char = "", quote = "\"", colClasses = "character"
)
mapping <- data.frame(
  probe = trimws(annotation$ID), gene = trimws(annotation[["Gene symbol"]])
)
mapping <- mapping[
  !is.na(mapping$gene) & grepl("^[A-Za-z][A-Za-z0-9_.-]*$", mapping$gene) &
    mapping$probe %in% rownames(geo) & mapping$gene %in% rownames(tcga),
]
stopifnot(nrow(mapping) > 0L, !anyDuplicated(mapping$probe))
mapping <- mapping[order(mapping$gene, mapping$probe, method = "radix"), ]
probe_groups <- split(match(mapping$probe, rownames(geo)), mapping$gene)
geo_by_gene <- t(vapply(probe_groups, function(rows) {
  if (length(rows) == 1L) {
    geo[rows, ]
  } else {
    apply(geo[rows, , drop = FALSE], 2L, median)
  }
}, FUN.VALUE = numeric(ncol(geo))))
colnames(geo_by_gene) <- colnames(geo)
