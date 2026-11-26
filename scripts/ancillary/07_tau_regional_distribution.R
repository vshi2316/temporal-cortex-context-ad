# FTP regional AD versus control distribution
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

root <- study_path("project")
map_file <- file.path(root, "results/13H1_aligned_regional_slope_map_release",
                      "01_aligned_34_region_contrasts.csv")
tau_file <- study_path("cohorts", "ADNI/data/PET/PET_Tau_Berkeley.csv")
clinical_file <- study_path("cohorts", "ADNI/data/Clinical/ADNIMERGE_03Jan2026.csv")
stopifnot(all(file.exists(c(map_file, tau_file, clinical_file))))

out <- file.path(root, "results", paste0("13H3_tau_spatial_specificity_",
                                          "release"))
dir.create(out, recursive = TRUE)

map <- fread(map_file)
stopifnot(nrow(map) == 34L, all(map$status == "SUCCESS"),
          !anyDuplicated(map$region), all(is.finite(map$estimate)),
          sum(map$temporal_target) == 6L)
setorder(map, region)
regions <- map$region
unit_weight <- function(x) {
  x <- x - mean(x)
  x / sqrt(sum(x^2))
}
map[, genetic_weight := unit_weight(-estimate)]
map[, temporal_weight := unit_weight(as.numeric(temporal_target))]
map[, specific_weight := unit_weight(genetic_weight -
                                     sum(genetic_weight * temporal_weight) * temporal_weight)]

left <- paste0("CTX_LH_", toupper(regions), "_SUVR")
right <- paste0("CTX_RH_", toupper(regions), "_SUVR")
pet <- fread(tau_file, select = c("RID", "LONIUID", "SCANDATE",
                                  "TRACER", left, right), showProgress = FALSE)
tracers <- pet[, .(scans = .N, people = uniqueN(RID)), by = TRACER]
fwrite(tracers, file.path(out, "01_tracer_counts.csv"))
total_scans <- nrow(pet)
pet[, `:=`(RID = as.character(RID), pet_date = as.IDate(SCANDATE))]
pet <- pet[toupper(trimws(TRACER)) == "FTP" & !is.na(pet_date)]
ftp_scans <- nrow(pet)
lh <- as.matrix(pet[, ..left])
rh <- as.matrix(pet[, ..right])
storage.mode(lh) <- "double"
storage.mode(rh) <- "double"
tau <- (lh + rh) / 2
valid <- rowSums(is.finite(tau) & tau > 0) == length(regions)
pet_meta <- pet[valid, .(RID, pet_image_id = as.character(LONIUID), pet_date)]
tau <- tau[valid, , drop = FALSE]
tau_wide <- as.data.table(tau)
setnames(tau_wide, regions)
pet <- cbind(pet_meta, tau_wide)
complete_scans <- nrow(pet)
setorder(pet, RID, pet_date, pet_image_id)
pet <- pet[, .SD[1L], by = RID]
tau <- as.matrix(pet[, ..regions])
pet[, `:=`(tau_global = rowMeans(tau),
           score_genetic = as.vector(tau %*% map$genetic_weight),
           score_specific = as.vector(tau %*% map$specific_weight))]

clinical <- fread(clinical_file,
                  select = c("RID", "EXAMDATE", "EXAMDATE_bl", "DX", "AGE", "PTGENDER"),
                  showProgress = FALSE)
clinical[, `:=`(RID = as.character(RID), exam_date = as.IDate(EXAMDATE),
                baseline_date = as.IDate(EXAMDATE_bl),
                DX = toupper(trimws(as.character(DX))))]
clinical <- clinical[!is.na(exam_date) & !is.na(DX) & nzchar(DX),
                     .(RID, exam_date, baseline_date, DX, AGE, PTGENDER)]
linked <- merge(pet, clinical, by = "RID", allow.cartesian = TRUE)
linked[, gap_days := abs(as.integer(exam_date - pet_date))]
linked <- linked[gap_days <= 180L]
setorder(linked, RID, gap_days, exam_date)
linked <- linked[, .SD[1L], by = RID]
matched_people <- nrow(linked)
linked[, diagnosis := fifelse(DX == "CN", "CN",
                              fifelse(DX %in% c("AD", "DEMENTIA"), "AD", "OTHER"))]
d <- linked[diagnosis %in% c("CN", "AD")]
d[, `:=`(case = as.integer(diagnosis == "AD"),
         age = as.numeric(AGE) + as.numeric(pet_date - baseline_date) / 365.25,
         sex = factor(PTGENDER))]
d <- droplevels(d[is.finite(age) & !is.na(sex)])

flow <- data.table(
  stage = c("source_scans", "FTP_scans", "FTP_complete_34_scans",
            "first_FTP_per_person", "clinical_visit_within_180d",
            "CN_complete", "AD_complete"),
  n = c(total_scans, ftp_scans, complete_scans, nrow(pet), matched_people,
        sum(d$case == 0L), sum(d$case == 1L))
)
fwrite(flow, file.path(out, "02_sample_flow.csv"))
if (length(unique(d$case)) != 2L || min(table(d$case)) < 10L)
  stop("Fewer than 10 complete FTP participants in CN or AD. See 02_sample_flow.csv")

case_effect <- function(response, covariates, analysis) {
  fit <- lm(reformulate(c("case", covariates), response), data = d)
  b <- coef(summary(fit))["case", ]
  ci <- confint(fit, "case")
  data.table(analysis, estimate_AD_minus_CN = unname(b["Estimate"]),
             lower_95 = ci[1L], upper_95 = ci[2L],
             p_value = unname(b["Pr(>|t|)"]), n = nobs(fit))
}
regional <- rbindlist(lapply(regions, function(region) {
  fit <- lm(reformulate(c("case", "age", "sex"), region), data = d)
  b <- coef(summary(fit))["case", ]
  data.table(region, tau_AD_minus_CN = unname(b["Estimate"]),
             tau_p_value = unname(b["Pr(>|t|)"]))
}))
regional[, tau_q_value := p.adjust(tau_p_value, method = "BH")]
regional <- merge(map[, .(region, estimate, q_value, temporal_target,
                          genetic_weight, specific_weight)], regional, by = "region")
fwrite(regional, file.path(out, "04_regional_tau_and_frozen_genetic_map.csv"))
fwrite(d[, .(RID, pet_image_id, pet_date, exam_date, gap_days, diagnosis,
             age, sex, tau_global, score_genetic, score_specific)],
       file.path(out, "05_subject_scores.csv"))

writeLines(capture.output(sessionInfo()),file.path(out,"99_sessionInfo.txt"))
