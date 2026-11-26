# GSE297004 donor-level neuronal program contrasts
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

## Frozen 13H32 gene universe; GSE297004 MTG nuclei; donor is the inference unit.
local({
  for (pkg in c("Matrix", "data.table", "sandwich"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)

  frozen <- file.path(REPO_ROOT, "definitions/mtg_gene_programs.csv")
  out <- file.path(study_path("project", "results"),
    paste0("13H34_GSE297004_MTG_neuron_bridge_", "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  status_file <- file.path(out, "00_status.txt")
  status <- function(label, note = NULL) {
    lines <- c(paste0("STATUS=", label), if (!is.null(note)) paste0("NOTE=", note),
      paste0("OUTPUT_DIR=", out))
    writeLines(lines, status_file)
    cat(paste(lines, collapse = "\n"), "\n")
  }
  put <- function(x, name) data.table::fwrite(x, file.path(out, name), na = "NA", bom = TRUE)
  status("RUNNING")

  meta_file <- study_path("mtg_metadata")
  count_file <- study_path("mtg_counts")
  if (!all(file.exists(c(meta_file, count_file, frozen)))) stop("Missing downloaded RDS or frozen 13H32 gene table")
  inputs <- c(metadata = meta_file, counts = count_file, frozen_program = frozen)
  put(data.frame(role = names(inputs), path = inputs,
    bytes = as.numeric(file.info(inputs)$size)), "00_input_provenance.csv")
  meta <- readRDS(meta_file)
  if (!is.data.frame(meta)) meta <- as.data.frame(meta)
  put(data.frame(column = names(meta), class = vapply(meta, function(x) class(x)[1], ""),
    example = vapply(meta, function(x) paste(head(unique(as.character(x[!is.na(x)])), 3), collapse = " | "), "")),
    "01_metadata_column_audit.csv")

  pick <- function(options) {
    canonical <- function(x) tolower(gsub("[^a-z0-9]", "", x))
    idx <- match(canonical(options), canonical(names(meta)))
    idx <- idx[!is.na(idx)]
    if (length(idx)) names(meta)[idx[1]] else NA_character_
  }
  donor_col <- "individual"
  diagnosis_col <- "diagnosis"
  cell_col <- "anno_lvl_2"
  if (!all(c(donor_col, diagnosis_col, cell_col, "barcode") %in% names(meta)))
    stop("GSE297004 metadata schema differs from 01_metadata_column_audit.csv")
  put(data.frame(role = c("donor", "diagnosis", "cell_subtype"),
    column = c(donor_col, diagnosis_col, cell_col)), "02_selected_metadata_columns.csv")
  cells <- data.frame(cell_id = as.character(meta$barcode), donor = as.character(meta[[donor_col]]),
    diagnosis_raw = as.character(meta[[diagnosis_col]]), subtype = as.character(meta[[cell_col]]))
  levels <- data.table::as.data.table(cells)[, .(n_nuclei = .N, n_donors = data.table::uniqueN(donor)), by = subtype]
  put(levels[order(-n_nuclei)], "03_author_subtype_levels.csv")
  is_vulnerable <- grepl("RORB", cells$subtype, ignore.case = TRUE) & grepl("LCN15", cells$subtype, ignore.case = TRUE)
  is_resilient <- grepl("RORB", cells$subtype, ignore.case = TRUE) & grepl("LINC01202", cells$subtype, ignore.case = TRUE)
  if (length(unique(cells$subtype[is_vulnerable])) != 1L ||
      length(unique(cells$subtype[is_resilient])) != 1L) {
    status("SUBTYPE_MAPPING_REQUIRED", "Inspect 03_author_subtype_levels.csv and the paper's RNA-to-IMC cluster mapping")
    stop("Published vulnerable/resilient RNA cluster names could not be matched uniquely. Output: ", out)
  }
  cells$target <- ifelse(is_vulnerable, "RORB_LCN15_vulnerable",
    ifelse(is_resilient, "RORB_LINC01202_resilient", NA_character_))
  dx <- tolower(trimws(cells$diagnosis_raw))
  cells$diagnosis <- ifelse(grepl("cntrl|ctrl|control|healthy|non.ad", dx), "control",
    ifelse(grepl("alz|(^|[^a-z])ad($|[^a-z])", dx), "AD", NA_character_))
  if (any(is.na(cells$diagnosis[!is.na(cells$target)]))) {
    put(unique(cells[, c("diagnosis_raw", "diagnosis")]), "04_diagnosis_mapping_review.csv")
    status("DIAGNOSIS_MAPPING_REQUIRED", "Review diagnosis labels before treating any donor as AD or control")
    stop("Diagnosis could not be mapped unambiguously. Output: ", out)
  }
  if (anyNA(cells$donor) || any(!nzchar(cells$donor))) stop("Missing donor ID")
  check <- unique(cells[, c("donor", "diagnosis")])
  if (anyDuplicated(check$donor)) stop("One donor has conflicting diagnosis labels")
  put(data.table::as.data.table(cells)[, .(n_nuclei = .N), by = .(donor, diagnosis, subtype, target)],
    "04_donor_subtype_capture_audit.csv")

  age_col <- pick(c("age", "age_at_death", "age_death"))
  sex_col <- pick(c("sex", "gender"))
  trem_col <- pick(c("trem2_genotype", "trem2_status", "trem2", "genotype"))
  trem_value <- if (!is.na(trem_col)) as.character(meta[[trem_col]]) else {
    x <- as.character(meta[[diagnosis_col]])
    ifelse(grepl("TREM2|R47H|R62H", x, ignore.case = TRUE), "risk",
      ifelse(grepl("CV", x, ignore.case = TRUE), "common", NA_character_))
  }
  donor_meta <- unique(data.frame(donor = cells$donor,
    age = if (!is.na(age_col)) suppressWarnings(as.numeric(as.character(meta[[age_col]]))) else NA_real_,
    sex = if (!is.na(sex_col)) as.character(meta[[sex_col]]) else NA_character_,
    trem2 = trem_value))
  if (anyDuplicated(donor_meta$donor)) {
    put(donor_meta, "05_donor_covariate_conflicts.csv")
    status("DONOR_COVARIATE_REVIEW_REQUIRED", "Age, sex or TREM2 varies within donor")
    stop("Donor covariates are inconsistent. Output: ", out)
  }
  donor_meta <- merge(donor_meta, check, by = "donor")
  put(donor_meta, "05_donor_metadata_audit.csv")

  cat("READ_COUNTS=", count_file, "\n", sep = "")
  counts <- readRDS(count_file)
  if (is.list(counts) && !inherits(counts, "Matrix") && !is.matrix(counts)) {
    if (!is.null(counts$counts)) counts <- counts$counts
  }
  if (!inherits(counts, "Matrix") && !is.matrix(counts)) stop("Counts RDS is not a gene-by-nucleus matrix")
  if (is.null(rownames(counts)) || is.null(colnames(counts))) stop("Count matrix has no gene or nucleus identifiers")
  if (nrow(counts) < 10000 || ncol(counts) < 1000) stop("Unexpected count matrix dimensions")
  if (anyNA(cells$cell_id) || anyDuplicated(cells$cell_id) || anyDuplicated(colnames(counts)) ||
      !setequal(cells$cell_id, colnames(counts))) {
    status("NUCLEUS_ID_LINK_UNRESOLVED", "Metadata IDs do not form a one-to-one match with count columns")
    stop("Metadata/count nucleus link failed. Output: ", out)
  }
  cells <- cells[match(colnames(counts), cells$cell_id), ]
  if (inherits(counts, "sparseMatrix")) {
    if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-6))
      stop("Counts are not finite nonnegative integer-like values")
  }

  pg <- data.table::fread(frozen)
  pg <- pg[program %in% c("gAD_priority_top200", "gMAV_priority_top200"), .(gene_id, program)]
  pg[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
  if (anyDuplicated(pg$gene_id)) stop("Frozen gene sets overlap; direct score needs review")
  feature <- sub("\\.[0-9]+$", "", rownames(counts))
  if (!all(grepl("^ENSG", feature))) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE))
      stop("Counts use gene symbols. Install Bioconductor packages org.Hs.eg.db and AnnotationDbi")
    symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
      keys = pg$gene_id, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
    pg[, feature_key := as.character(symbols[gene_id])]
  } else pg[, feature_key := gene_id]
  pg[, count_row := match(feature_key, feature)]
  pg[, present := !is.na(count_row)]
  put(pg, "06_frozen_gene_feature_coverage.csv")
  coverage <- pg[, .(frozen_genes = .N, present_genes = sum(present),
    coverage = mean(present)), by = program]
  put(coverage, "07_program_coverage_summary.csv")
  if (any(coverage$coverage < .8)) {
    status("GENE_COVERAGE_INSUFFICIENT", "At least one frozen program has less than 80% feature coverage")
    stop("Frozen program coverage failed. Output: ", out)
  }
  pg <- pg[pg$present == TRUE]
  if (anyDuplicated(pg$count_row)) stop("Multiple frozen genes map to the same count row")

  target_idx <- which(!is.na(cells$target))
  group <- paste(cells$donor[target_idx], cells$target[target_idx], sep = "||")
  group_levels <- unique(group)
  design <- Matrix::sparseMatrix(i = seq_along(group), j = match(group, group_levels), x = 1,
    dims = c(length(group), length(group_levels)))
  sub_counts <- counts[pg$count_row, target_idx, drop = FALSE]
  pb <- as.matrix(sub_counts %*% design)
  rownames(pb) <- pg$gene_id
  colnames(pb) <- group_levels
  all_library_totals <- Matrix::colSums(counts)
  depth <- as.numeric(Matrix::crossprod(design, all_library_totals[target_idx]))
  nuclei <- as.numeric(Matrix::colSums(design))
  grp <- cells[target_idx[match(group_levels, group)], c("donor", "diagnosis", "target")]
  grp$n_nuclei <- nuclei
  grp$full_UMI <- depth
  grp$captured_fraction_among_neurons <- nuclei / as.numeric(table(cells$donor)[grp$donor])
  grp <- merge(grp, donor_meta[, c("donor", "age", "sex", "trem2")], by = "donor", sort = FALSE)
  grp <- grp[match(group_levels, paste(grp$donor, grp$target, sep = "||")), ]
  if (anyNA(grp$donor)) stop("Pseudobulk metadata aggregation failed")
  put(grp, "08_donor_subtype_pseudobulk_QC.csv")
  keep <- grp$n_nuclei >= 20 & grp$full_UMI > 0
  grp <- grp[keep, ]
  pb <- pb[, keep, drop = FALSE]
  if (!ncol(pb)) stop("No donor-subtype pseudobulk passes 20 nuclei")
  logcpm <- log2(sweep(pb, 2, grp$full_UMI / 1e6, "/") + 1)
  gene_sd <- apply(logcpm, 1, sd)
  usable <- is.finite(gene_sd) & gene_sd > 0
  pg <- pg[match(rownames(logcpm), pg$gene_id)]
  pg[, usable := usable]
  put(pg[, .(gene_id, program, present, usable)], "09_scored_gene_audit.csv")
  z <- sweep(sweep(logcpm[usable, , drop = FALSE], 1, rowMeans(logcpm[usable, , drop = FALSE]), "-"),
    1, gene_sd[usable], "/")
  for (label in c("gAD_priority_top200", "gMAV_priority_top200")) {
    ix <- pg$program[usable] == label
    if (sum(ix) < 0.8 * sum(pg$program == label)) stop("Too few variable genes in ", label)
    grp[[if (label == "gAD_priority_top200") "gAD" else "gMAV"]] <- colMeans(z[ix, , drop = FALSE])
  }
  grp$direct <- grp$gAD - grp$gMAV
  put(grp, "10_donor_subtype_program_scores.csv")

  wide <- data.table::dcast(data.table::as.data.table(grp), donor + diagnosis ~ target,
    value.var = c("gAD", "gMAV", "direct"))
  wide <- merge(wide, donor_meta[, c("donor", "age", "sex", "trem2")], by = "donor", sort = FALSE)
  v <- "RORB_LCN15_vulnerable"; r <- "RORB_LINC01202_resilient"
  required <- as.vector(outer(c("gAD", "gMAV", "direct"), c(v, r), paste, sep = "_"))
  if (!all(required %in% names(wide))) {
    status("PAIRED_SUBTYPE_COVERAGE_INSUFFICIENT", "At least one published target subtype has no pseudobulk passing the nucleus filter")
    stop("Vulnerable or resilient subtype scores missing. Output: ", out)
  }
  wide <- as.data.frame(wide)
  wide$paired_direct <- wide[[paste0("direct_", v)]] - wide[[paste0("direct_", r)]]
  wide$paired_gAD <- wide[[paste0("gAD_", v)]] - wide[[paste0("gAD_", r)]]
  wide$paired_gMAV <- wide[[paste0("gMAV_", v)]] - wide[[paste0("gMAV_", r)]]
  wide <- wide[is.finite(wide$paired_direct), ]
  wide$AD01 <- as.integer(wide$diagnosis == "AD")
  put(wide, "11_paired_donor_contrasts.csv")
  counts_by_dx <- table(wide$diagnosis)
  if (!all(c("AD", "control") %in% names(counts_by_dx)) || any(counts_by_dx[c("AD", "control")] < 8)) {
    status("PAIRED_DONORS_INSUFFICIENT", "Primary comparison needs at least eight complete donors per diagnosis group")
    stop("Insufficient paired donor coverage. Output: ", out)
  }

  model <- function(y, adjusted = FALSE) {
    d <- wide
    rhs <- "AD01"
    if (adjusted) {
      d$TREM2_risk <- as.integer(grepl("R47H|R62H|risk", d$trem2, ignore.case = TRUE))
      d$TREM2_risk[is.na(d$trem2)] <- NA_integer_
      covars <- c("age", "sex", "TREM2_risk")
      covars <- covars[vapply(covars, function(a) {
        x <- d[[a]]
        sum(!is.na(x)) >= .8 * nrow(d) && length(unique(x[!is.na(x)])) > 1L
      }, logical(1))]
      if (!length(covars)) return(NULL)
      rhs <- paste(c("AD01", covars), collapse = " + ")
    }
    f <- stats::as.formula(paste(y, "~", rhs))
    fit <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    if (fit$rank != length(stats::coef(fit)) || stats::df.residual(fit) < 8) return(NULL)
    se <- sqrt(sandwich::vcovHC(fit, type = "HC3")["AD01", "AD01"])
    beta <- unname(stats::coef(fit)["AD01"])
    if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NULL)
    data.frame(outcome = y, model = if (adjusted) "covariate_sensitivity" else "primary_unadjusted",
      n_donors = stats::nobs(fit), AD = sum(stats::model.frame(fit)$AD01 == 1),
      control = sum(stats::model.frame(fit)$AD01 == 0),
      estimate_AD_minus_control = beta, se_HC3 = se,
      ci_low = beta - stats::qt(.975, stats::df.residual(fit)) * se,
      ci_high = beta + stats::qt(.975, stats::df.residual(fit)) * se,
      p = 2 * stats::pt(-abs(beta / se), stats::df.residual(fit)),
      formula = paste(deparse(f), collapse = ""))
  }
  outcomes <- c("paired_direct", "paired_gAD", "paired_gMAV",
    paste0("direct_", v), paste0("direct_", r))
  fits <- data.table::rbindlist(c(lapply(outcomes, model),
    lapply(outcomes, model, adjusted = TRUE)), fill = TRUE)
  if (!nrow(fits) || !any(fits$outcome == "paired_direct" & fits$model == "primary_unadjusted"))
    stop("Primary donor-level contrast was not estimable. Output: ", out)
  fits[, q_secondary := NA_real_]
  fits[model == "primary_unadjusted" & outcome != "paired_direct",
    q_secondary := p.adjust(p, "BH")]
  put(fits, "12_donor_level_models.csv")
  put(data.frame(metric = c("donors_total", "donors_AD", "donors_control", "gAD_gMAV_score_correlation"),
    value = c(nrow(wide), sum(wide$AD01 == 1), sum(wide$AD01 == 0),
      cor(grp$gAD, grp$gMAV))), "13_interpretation_checks.csv")

  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  status("GSE297004_FROZEN_NEURON_BRIDGE_COMPLETE",
    paste0("PAIRED_DONORS=", nrow(wide), "; PRIMARY_P=", signif(fits$p[fits$outcome == "paired_direct" & fits$model == "primary_unadjusted"], 4),
      "; IMC_SPATIAL_LINK=NOT_ANALYZED"))
})
