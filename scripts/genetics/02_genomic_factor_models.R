# Genomic structural equation models and signed loadings
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

#!/usr/bin/env Rscript

# Step 09E: compare prespecified Genomic SEM measurement structures for the
# four MAV domains, followed by a separate AD genetic bridge. This script uses
# the frozen Step 09D LDSC covariance only.

options(stringsAsFactors = FALSE, scipen = 999, width = 180, warn = 1)

PROJECT_DIR <- study_path("project")
INPUT_DIR <- file.path(PROJECT_DIR, "results", "09D_GWAS_QC_LDSC")
RESULT_DIR <- file.path(PROJECT_DIR, "results", "09E_GenomicSEM_model_comparison")
INPUT_RDS <- file.path(INPUT_DIR, "05_primary_five_trait_covstruc.rds")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(RESULT_DIR, "09E_MAV_GenomicSEM_model_comparison.log")
if (file.exists(LOG_FILE)) unlink(LOG_FILE, force = TRUE)

log_line <- function(level = "INFO", ...) {
  line <- sprintf(
    "[%s] [%s] %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level,
    paste0(..., collapse = "")
  )
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

write_csv <- function(x, filename, row.names = FALSE) {
  path <- file.path(RESULT_DIR, filename)
  utils::write.csv(x, path, row.names = row.names, na = "", fileEncoding = "UTF-8")
  log_line("INFO", "Wrote ", path, " (", nrow(x), " rows)")
}

if (!file.exists(INPUT_RDS)) {
  stop("Missing Step 09D covariance structure: ", INPUT_RDS, call. = FALSE)
}
if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  stop("GenomicSEM is not installed.", call. = FALSE)
}
if (!"usermodel" %in% getNamespaceExports("GenomicSEM")) {
  stop("The installed GenomicSEM does not export usermodel().", call. = FALSE)
}

covstruc <- readRDS(INPUT_RDS)
required_parts <- c("S", "V")
if (!all(required_parts %in% names(covstruc))) {
  stop("The Step 09D object does not contain both S and V.", call. = FALSE)
}

trait_names <- c("CMM", "INSOMNIA", "HEARING", "GRIP", "AD")
if (!identical(dim(covstruc$S), c(5L, 5L))) {
  stop("Expected a 5 x 5 genetic covariance matrix.", call. = FALSE)
}
rownames(covstruc$S) <- colnames(covstruc$S) <- trait_names

subset_covstruc <- function(x, keep_names) {
  keep <- match(keep_names, trait_names)
  if (anyNA(keep)) stop("Unknown trait in covariance subset.", call. = FALSE)
  all_pairs <- which(lower.tri(matrix(TRUE, nrow = length(trait_names),
                                      ncol = length(trait_names)), diag = TRUE),
                     arr.ind = TRUE)
  keep_pair <- which(all_pairs[, 1] %in% keep & all_pairs[, 2] %in% keep)
  out <- x
  out$S <- x$S[keep, keep, drop = FALSE]
  out$V <- x$V[keep_pair, keep_pair, drop = FALSE]
  rownames(out$S) <- colnames(out$S) <- keep_names
  out
}

covstruc_mav <- subset_covstruc(covstruc, trait_names[1:4])

if (any(!is.finite(covstruc$S)) || any(!is.finite(covstruc$V))) {
  stop("S or V contains non-finite values.", call. = FALSE)
}
if (max(abs(covstruc$S - t(covstruc$S))) > 1e-8) {
  stop("S is not symmetric.", call. = FALSE)
}

eigenvalues <- eigen(covstruc$S, symmetric = TRUE, only.values = TRUE)$values
qc <- data.frame(
  check = c("S_dimension", "V_dimension", "S_min_eigenvalue", "S_positive_definite"),
  value = c(
    paste(dim(covstruc$S), collapse = "x"),
    paste(dim(covstruc$V), collapse = "x"),
    format(min(eigenvalues), scientific = TRUE),
    as.character(min(eigenvalues) > 0)
  )
)
write_csv(qc, "01_input_covariance_QC.csv")
if (min(eigenvalues) <= 0) {
  stop("The genetic covariance matrix is not positive definite.", call. = FALSE)
}

sd_g <- sqrt(diag(covstruc$S))
rg <- covstruc$S / outer(sd_g, sd_g)
write_csv(
  data.frame(trait = trait_names, rg, check.names = FALSE),
  "02_genetic_correlation_matrix.csv"
)

# The primary factor uses four traits and excludes AD. Loadings retain their estimated signs.
model_m1 <- paste(
  "gMAV =~ CMM + INSOMNIA + HEARING + GRIP",
  sep = "\n"
)

# M3 estimates mutually adjusted domain-specific genetic paths to AD.
model_m3 <- paste(
  "AD ~ CMM + INSOMNIA + HEARING + GRIP",
  "CMM ~~ INSOMNIA + HEARING + GRIP",
  "INSOMNIA ~~ HEARING + GRIP",
  "HEARING ~~ GRIP",
  sep = "\n"
)

# M4 estimates the association of the single multisystem factor with AD.
model_m4 <- paste(
  "gMAV =~ CMM + INSOMNIA + HEARING + GRIP",
  "gMAV ~~ AD",
  sep = "\n"
)

model_specs <- list(
  M1_single_gMAV = model_m1,
  M3_domain_specific_AD = model_m3,
  M4_gMAV_AD_bridge = model_m4
)

writeLines(
  unlist(Map(function(name, syntax) c(paste0("[", name, "]"), syntax, ""),
             names(model_specs), model_specs)),
  file.path(RESULT_DIR, "03_frozen_model_syntax.txt"),
  useBytes = TRUE
)

fit_one <- function(model_name, syntax) {
  log_line("INFO", "Fitting ", model_name)
  warnings_seen <- character()
  model_covstruc <- if (model_name == "M1_single_gMAV") {
    covstruc_mav
  } else {
    covstruc
  }
  fit <- tryCatch(
    withCallingHandlers(
      GenomicSEM::usermodel(
        covstruc = model_covstruc,
        estimation = "DWLS",
        model = syntax,
        CFIcalc = TRUE,
        std.lv = TRUE,
        imp_cov = FALSE,
        fix_resid = TRUE,
        toler = NULL,
        Q_Factor = FALSE
      ),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    log_line("ERROR", model_name, ": ", conditionMessage(fit))
    return(list(
      status = data.frame(model = model_name, status = "FAILED", message = conditionMessage(fit)),
      fit = NULL, warnings = warnings_seen
    ))
  }
  if (is.null(fit)) {
    message <- "GenomicSEM::usermodel() returned NULL; model was not estimable."
    log_line("ERROR", model_name, ": ", message)
    return(list(
      status = data.frame(model = model_name, status = "FAILED", message = message),
      fit = NULL, warnings = warnings_seen
    ))
  }

  saveRDS(fit, file.path(RESULT_DIR, paste0("04_", model_name, ".rds")))
  log_line("INFO", model_name, " completed")
  list(
    status = data.frame(model = model_name, status = "SUCCESS", message = ""),
    fit = fit, warnings = warnings_seen
  )
}

fits <- Map(fit_one, names(model_specs), model_specs)
status <- do.call(rbind, lapply(fits, `[[`, "status"))
write_csv(status, "05_model_status.csv")

successful <- fits[vapply(fits, function(x) !is.null(x$fit), logical(1))]
if (!length(successful)) {
  stop("All prespecified Genomic SEM models failed. See 05_model_status.csv.", call. = FALSE)
}

extract_component <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) return(NULL)
  value <- x[[hit[[1]]]]
  if (is.null(value)) return(NULL)
  as.data.frame(value, stringsAsFactors = FALSE)
}

parameter_tables <- list()
fit_tables <- list()
warning_tables <- list()

for (model_name in names(successful)) {
  obj <- successful[[model_name]]$fit
  parameters <- extract_component(obj, c("results", "Results", "estimates", "parameters"))
  fit_stats <- extract_component(obj, c("modelfit", "model.fit", "fit", "fitmeasures"))

  if (!is.null(parameters)) {
    parameters$model <- model_name
    parameter_tables[[model_name]] <- parameters[, c("model", setdiff(names(parameters), "model"))]
  }
  if (!is.null(fit_stats)) {
    fit_stats$model <- model_name
    fit_tables[[model_name]] <- fit_stats[, c("model", setdiff(names(fit_stats), "model"))]
  }
  warnings_seen <- unique(successful[[model_name]]$warnings)
  if (length(warnings_seen)) {
    warning_tables[[model_name]] <- data.frame(model = model_name, warning = warnings_seen)
  }
}

if (length(parameter_tables)) {
  write_csv(do.call(rbind, parameter_tables), "06_all_model_parameters.csv")
}
if (length(fit_tables)) {
  write_csv(do.call(rbind, fit_tables), "07_model_fit_comparison.csv")
}
if (length(warning_tables)) {
  write_csv(do.call(rbind, warning_tables), "08_model_warnings.csv")
} else {
  write_csv(data.frame(model = character(), warning = character()), "08_model_warnings.csv")
}

component_audit <- do.call(rbind, lapply(names(successful), function(model_name) {
  data.frame(
    model = model_name,
    returned_components = paste(names(successful[[model_name]]$fit), collapse = ";")
  )
}))
write_csv(component_audit, "09_returned_component_audit.csv")

if (status$status[status$model == "M1_single_gMAV"] != "SUCCESS") {
  stop("The primary gMAV measurement model failed. See 05_model_status.csv.", call. = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(RESULT_DIR, "11_sessionInfo.txt"), useBytes = TRUE)
log_line("INFO", "Genomic SEM model estimation completed.")
