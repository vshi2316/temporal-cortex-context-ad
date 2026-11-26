# Three-scan current-state and prior-loss event models
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

local({
  options(stringsAsFactors = FALSE, scipen = 999)
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })

  root <- study_path("project", "results")
  input <- file.path(root,
    "13H55_transition_ordered_structural_chain_release")
  out <- file.path(root, paste0("13H56_event_proximity_state_vs_rate_",
                                 "release"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, name) fwrite(x, file.path(out, name), bom = TRUE)
  z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

  d <- fread(file.path(input, "06_first_triple_event_frame.csv"))
  frozen <- fread(file.path(input, "01_frozen_frame.csv"),
                  select = c("RID", "mri_date", "global_z"))
  g1 <- copy(frozen)
  g2 <- copy(frozen)
  setnames(g1, c("mri_date", "global_z"), c("t1", "global1"))
  setnames(g2, c("mri_date", "global_z"), c("t2", "global2"))
  d <- merge(d, g1, by = c("RID", "t1"), all.x = TRUE, sort = FALSE)
  d <- merge(d, g2, by = c("RID", "t2"), all.x = TRUE, sort = FALSE)
  if (anyDuplicated(d$RID)) stop("The event frame is no longer one row per participant")

  d[, `:=`(loss01_1sd = z(loss01), loss12_1sd = z(loss12),
            cog12_1sd = z(cog12),
            global_loss12 = (global1 - global2) / dt12,
            stage_t2 = factor(stage_t2), site = factor(site), sex = factor(sex))]
  d[, global_loss12_1sd := z(global_loss12)]
  needed <- c("RID", "followup_years", "event", "thickness2", "global2",
              "loss01_1sd", "loss12_1sd", "global_loss12_1sd",
              "cog12_1sd", "cognition2", "age2", "sex", "education",
              "gAD", "gMAV", "dt01", "dt12", "dtc12", "stage_t2", "site")
  d <- droplevels(d[complete.cases(d[, ..needed]) & followup_years > 0])
  flow <- data.table(people = nrow(d), events = sum(d$event),
                     stage_CN = sum(d$stage_t2 == "CN"),
                     stage_MCI = sum(d$stage_t2 == "MCI"))
  write_csv(flow, "00_common_sample.csv")
  if (nrow(d) < 100L || sum(d$event) < 20L)
    stop("The common sample cannot support the event models")

  core <- paste(
    "thickness2 + loss01_1sd + loss12_1sd + cognition2 +",
    "age2 + sex + education + gAD + gMAV + dt01 + dt12 +",
    "strata(stage_t2) + cluster(site)"
  )
  rhs <- list(
    structural_state_and_two_intervals = core,
    plus_recent_cognition = paste(core, "+ cog12_1sd + dtc12"),
    plus_global_cortex = paste(core,
      "+ cog12_1sd + dtc12 + global2 + global_loss12_1sd")
  )
  fits <- lapply(rhs, function(x) {
    coxph(as.formula(paste("Surv(followup_years, event) ~", x)),
          data = d, ties = "efron")
  })

  coefficients <- rbindlist(lapply(names(fits), function(model) {
    fit <- fits[[model]]
    b <- coef(fit)
    se <- sqrt(diag(vcov(fit)))
    data.table(model = model, term = names(b), log_HR = as.numeric(b),
               HR = exp(as.numeric(b)), lower = exp(as.numeric(b) - 1.96 * se),
               upper = exp(as.numeric(b) + 1.96 * se),
               p = 2 * pnorm(-abs(b / se)), people = nrow(d),
               events = sum(d$event))
  }), fill = TRUE)
  write_csv(coefficients, "01_event_cox_tiers.csv")

  contrasts <- rbindlist(lapply(names(fits), function(model) {
    fit <- fits[[model]]
    b <- coef(fit)
    L <- setNames(rep(0, length(b)), names(b))
    L["loss12_1sd"] <- 1
    L["loss01_1sd"] <- -1
    est <- sum(L * b)
    se <- sqrt(as.numeric(t(L) %*% vcov(fit) %*% L))
    data.table(model = model, comparison = "recent_minus_earlier_loss",
               log_HR_difference = est, se = se,
               lower = est - 1.96 * se, upper = est + 1.96 * se,
               p = 2 * pnorm(-abs(est / se)))
  }))
  write_csv(contrasts, "02_recency_contrast.csv")
  writeLines(c(
    "Same first-triple risk set in every Cox tier.",
    "Loss01 is MRI t0-t1; loss12 is MRI t1-t2; both are standardized within this event sample.",
    "All outcomes occur after MRI t2 and its matched cognitive examination.",
    "A recency contrast tests the recent minus earlier loss coefficient directly.",
    "The global-cortex tier tests whether six-region state and recent loss retain information beyond overall cortical thickness."
  ), file.path(out, "03_readme.txt"))
  capture.output(sessionInfo(), file = file.path(out, "99_sessionInfo.txt"))
  cat("13H56 complete. Results: ", out, "\n", sep = "")
})
