###############################################################################
## generate_pseudo_data
##
## Generates fully synthetic, internally-consistent pseudo data that mimics the
## real BD-Tau TBI / SCI / Control dataset, for public release on GitHub.
## No real patient values are used: every value is invented from summary-level
## distributions.
##
## HOW THIS FILE IS ORGANISED -------------------------------------------------
##   PART 1  Inspect the REAL objects (only runs if they are loaded now).
##   SNAPSHOT  Copy the reals into a side env BEFORE the generator overwrites
##             the global names -- this is what makes verification meaningful.
##   PART 2  The generator (overwrites the global names with synthetics).
##   PART 3  Verify synthetics against the snapshot (names + types).
## Run standalone (no reals in session) and PART 1 / PART 3 simply no-op, so the
## committed repo script still works for anyone who clones it.
##
## DESIGN PRINCIPLE -----------------------------------------------------------
## Most columns are DETERMINISTIC functions of a few primitives (BDTau,
## time-of-draw, Age, Sex, Cohort). We generate only the primitives per
## subject/visit and DERIVE the rest, building *_demo / *_baseline / *_match by
## JOINING off one Subject.ID-keyed registry. This preserves:
##   *_log2_BDTau == log2(*_BDTau);  mean_*_log2_BDTau == log2(mean_*_BDTau)
##   *_Raw_12Hr_Value == 2*Raw_Day;  *_First_Digit_* == floor(...)
##   mean_BDTau_Ratio == mean_serum_BDTau / mean_csf_BDTau
##   Age_40/65_Recode, baseline_*, and cross-table subject keys
## so Kruskal-Wallis / Dunn / median fold-change run and stay cohort-separated.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr)
})

## ============================================================================
## PART 1 -- INSPECT REAL OBJECTS  (guarded: only if they exist in this session)
## ============================================================================

real_obj_names <- c(
  "cntr",
  "sci_csf", "sci_csf_demo", "sci_serum", "sci_serum_demo", "sci_demo",
  "tbi_csf", "tbi_csf_demo", "tbi_serum", "tbi_serum_demo", "tbi_demo",
  ls()[grepl("_cleaned", ls())]
)
# keep only those actually present now, so the script is standalone-safe
real_obj_names <- real_obj_names[vapply(real_obj_names, exists, logical(1))]

if (length(real_obj_names)) {
  cat("\n--- inspecting", length(real_obj_names), "real objects ---\n")
  overview <- purrr::map_dfr(real_obj_names, function(nm) {
    x <- get(nm)
    tibble::tibble(
      object = nm,
      class  = paste(class(x), collapse = "/"),
      nrow   = if (is.data.frame(x)) nrow(x) else length(x),
      ncol   = if (is.data.frame(x)) ncol(x) else NA_integer_
    )
  })
  print(overview)
  for (nm in real_obj_names) {
    cat("\n=====", nm, "=====\n")
    str(get(nm), max.level = 1)
    if (is.data.frame(get(nm))) print(utils::head(get(nm), 3))
  }
}

## ---- SNAPSHOT reals into a side env BEFORE the generator overwrites them ----

.real_env <- new.env()
for (nm in real_obj_names) assign(nm, get(nm), envir = .real_env)
if (length(real_obj_names))
  cat("\nSnapshotted", length(real_obj_names), "real objects into .real_env\n")

set.seed(42)

## ============================================================================
## PART 2 -- GENERATOR
## ============================================================================

## 0. CONFIG  (tune subject counts / visit counts / magnitudes here) ----------
cfg <- list(
  n_control = 22,
  n_sci = 57,
  n_tbi = 298,
  tbi_severe_frac = 0.72,
  
  visits_sci_csf = function(n) pmax(1L, rpois(n, 4)),
  visits_sci_serum = function(n) pmax(1L, rpois(n, 7)),
  visits_tbi_csf = function(n) pmax(1L, rpois(n, 12)),
  visits_tbi_serum = function(n) pmax(1L, rpois(n, 6)),
  
  bdtau = list(
    Control = list(serum = c(med =   6, sd = 0.5)),
    SCI = list(serum = c(med =  40, sd = 0.8), csf = c(med = 1100, sd = 0.9)),
    `Moderate TBI` = list(serum = c(med =  70, sd = 0.8), csf = c(med = 3000, sd = 0.9)),
    `Severe TBI` = list(serum = c(med = 120, sd = 0.8), csf = c(med = 8000, sd = 0.9))
  ),
  decay_per_day = 0.06
)

processed_data_path <- "~/Desktop/Projects/manuscript writing/TBI BD-tau/table"

## 1. DERIVATION HELPERS ------------------------------------------------------

recode_age <- function(age, cut) ifelse(age <= cut, paste0("<=", cut), paste0(">", cut))

add_derived <- function(df, pfx) {
  raw_day <- paste0(pfx, "_Raw_Day_Value"); raw_12 <- paste0(pfx, "_Raw_12Hr_Value")
  fd24 <- paste0(pfx, "_First_Digit_24Hr"); fd12 <- paste0(pfx, "_First_Digit_12Hr")
  bd <- paste0(pfx, "_BDTau"); l2 <- paste0(pfx, "_log2_BDTau")
  m_raw <- paste0("mean_", pfx, "_Raw_Day_Value")
  m_bd <- paste0("mean_", pfx, "_BDTau"); m_l2 <- paste0("mean_", pfx, "_log2_BDTau")
  
  df[[raw_12]] <- 2 * df[[raw_day]]
  df[[fd24]] <- as.character(floor(df[[raw_day]]))
  df[[fd12]] <- as.character(floor(df[[raw_12]]))
  df[[l2]] <- log2(df[[bd]])
  
  df %>%
    group_by(Subject.ID, day_idx) %>%
    mutate(!!m_raw := round(mean(.data[[raw_day]]), 2),
           !!m_bd := mean(.data[[bd]]),
           !!m_l2 := log2(mean(.data[[bd]]))) %>%
    ungroup()
}

draw_traj <- function(n_visits, med, sd, decay) {
  days <- sort(sample(0:max(0, n_visits + 2), n_visits, replace = FALSE))
  base <- rlnorm(1, meanlog = log(med), sdlog = sd)
  vals <- base * exp(-decay * days) * rlnorm(n_visits, 0, 0.25)
  tibble(day_idx = days, raw_day = days + runif(n_visits, 0.35, 0.99), bdtau = vals)
}

## 2. SUBJECT REGISTRIES ------------------------------------------------------

sci_loc <- c("C4","C5","C6","T4","T6","T10","L1", NA)
sci_locgrp <- c(C4="Cervical",C5="Cervical",C6="Cervical",
                T4="Thoracic",T6="Thoracic",T10="Thoracic",L1="Lumbar")
asia_levels <- c("A","B","C/D", NA)

sci_subjects <- tibble(
  Subject.ID = sprintf("CP%03d", 4:(4 + cfg$n_sci - 1)),
  Study = "CP", Group = "SCI", Cohort = "SCI",
  Age = round(rnorm(cfg$n_sci, 50, 18)) %>% pmin(90) %>% pmax(18),
  Sex = sample(c("Male","Female"), cfg$n_sci, TRUE, c(.65,.35)),
  ASIA = sample(asia_levels, cfg$n_sci, TRUE),
  Location = sample(sci_loc, cfg$n_sci, TRUE),
  Coenrolled = ifelse(runif(cfg$n_sci) < .2, sprintf("ST-%04d", sample(1000:1999, cfg$n_sci, TRUE)), NA),
  In_Hospital_Dead = sample(c("Survived","Dead"), cfg$n_sci, TRUE, c(.85,.15)),
  Death_06 = sample(c("Survived","Dead", NA), cfg$n_sci, TRUE, c(.6,.15,.25))
) %>%
  mutate(Location_Group = unname(sci_locgrp[Location]),
         Age_40_Recode = recode_age(Age, 40),
         Age_65_Recode = recode_age(Age, 65))

gcs_group <- function(g) cut(g, c(0,2,4,8,12,15),
                             labels = c("GCS 1-2","GCS 3-4","GCS 5-8","GCS 9-12","GCS 13-15"))
n_sev <- round(cfg$n_tbi * cfg$tbi_severe_frac)
tbi_subjects <- tibble(
  Subject.ID = as.character(1000:(1000 + cfg$n_tbi - 1)),
  Study = "BTRC",
  Group  = c(rep("Severe", n_sev), rep("Moderate", cfg$n_tbi - n_sev)),
  Cohort = c(rep("Severe TBI", n_sev), rep("Moderate TBI", cfg$n_tbi - n_sev)),
  Age = round(rnorm(cfg$n_tbi, 42, 18)) %>% pmin(90) %>% pmax(16),
  Sex = sample(c("Male","Female"), cfg$n_tbi, TRUE, c(.72,.28)),
  GCS = c(sample(3:8, n_sev, TRUE), sample(9:12, cfg$n_tbi - n_sev, TRUE)),
  GOSE_01 = NA, Note_from_Email = NA
) %>%
  mutate(GCS_Group = as.character(gcs_group(GCS)),
         Age_40_Recode = recode_age(Age, 40),
         Age_65_Recode = recode_age(Age, 65),
         .dead = sample(c(TRUE, FALSE), cfg$n_tbi, TRUE, c(.3,.7)))

for (col in c("GOS_03","GOS_06","GOS_12","GOS_24","Death_03","Death_06","Death_12","Death_24"))
  tbi_subjects[[col]] <- ifelse(tbi_subjects$.dead, "Dead",
                                sample(c("Survived", NA), cfg$n_tbi, TRUE, c(.7,.3)))
tbi_subjects$.dead <- NULL

ctrl_subjects <- tibble(
  Study = "CTE Internal", Group = "Controls", Cohort = "Control",
  Age = round(rnorm(cfg$n_control, 35, 12)) %>% pmin(70) %>% pmax(18),
  Sex = sample(c("Male","Female"), cfg$n_control, TRUE),
  Subject.ID = sprintf("CTE-%03d", sample(1:99, cfg$n_control))
) %>%
  mutate(Age_40_Recode = recode_age(Age, 40), Age_65_Recode = recode_age(Age, 65))

## 3. LONGITUDINAL MEASUREMENTS -----------------------------------------------

make_long <- function(subjects, pfx, n_visits_fun, visit_fmt) {
  par <- cfg$bdtau
  rows <- pmap_dfr(subjects, function(...) {
    s <- tibble(...); p <- par[[s$Cohort]][[pfx]]
    if (is.null(p)) return(NULL)
    nv <- n_visits_fun(1); tr <- draw_traj(nv, p["med"], p["sd"], cfg$decay_per_day)
    tibble(Subject.ID = s$Subject.ID,
           Subject.Visit = visit_fmt(s$Subject.ID, seq_len(nv)),
           day_idx = tr$day_idx,
           !!paste0(pfx, "_Raw_Day_Value") := tr$raw_day,
           !!paste0(pfx, "_BDTau")         := tr$bdtau,
           Cohort = s$Cohort)
  })
  add_derived(rows, pfx)
}

fmt_sci_csf <- function(id, k) sprintf("%s.C%02d%s", id, k, sample(LETTERS, length(k), TRUE))
fmt_sci_serum <- function(id, k) sprintf("%s.S%02d%s", id, k, sample(LETTERS, length(k), TRUE))
fmt_tbi <- function(id, k) as.numeric(sprintf("%s.%02d", id, k))

sci_csf_long <- make_long(sci_subjects, "csf", cfg$visits_sci_csf, fmt_sci_csf)
sci_serum_long <- make_long(sci_subjects, "serum", cfg$visits_sci_serum, fmt_sci_serum)
tbi_csf_long <- make_long(tbi_subjects, "csf", cfg$visits_tbi_csf, fmt_tbi)
tbi_serum_long <- make_long(tbi_subjects, "serum", cfg$visits_tbi_serum, fmt_tbi)

## 4. ASSEMBLE NAMED OBJECTS --------------------------------------------------

cntr <- ctrl_subjects %>%
  mutate(serum_BDTau = rlnorm(n(), log(cfg$bdtau$Control$serum["med"]), cfg$bdtau$Control$serum["sd"]),
         serum_First_Digit_24Hr = 180,
         serum_log2_BDTau = log2(serum_BDTau),
         Subject.Visit = sprintf("%s-%02dS", Subject.ID, sample(1:12, n(), TRUE)),
         serum_Raw_Day_Value = NA, serum_Raw_12Hr_Value = NA,
         serum_First_Digit_12Hr = 360, mean_serum_Raw_Day_Value = NA,
         mean_serum_BDTau = serum_BDTau, mean_serum_log2_BDTau = serum_log2_BDTau,
         Cohort = "Control") %>%
  select(Study, Group, Age, Sex, Subject.ID, serum_BDTau, serum_First_Digit_24Hr,
         serum_log2_BDTau, Subject.Visit, Age_40_Recode, Age_65_Recode,
         serum_Raw_Day_Value, serum_Raw_12Hr_Value, serum_First_Digit_12Hr,
         mean_serum_Raw_Day_Value, mean_serum_BDTau, mean_serum_log2_BDTau, Cohort)

sci_csf <- sci_csf_long
sci_serum <- sci_serum_long
tbi_csf <- tbi_csf_long %>% mutate(csf_Comment = NA)
tbi_serum <- tbi_serum_long %>% mutate(serum_Comment = NA)

sci_demo <- sci_subjects %>%
  select(Subject.ID, Study, Group, Cohort, Age, Age_40_Recode, Age_65_Recode,
         Coenrolled, Sex, Location, Location_Group, ASIA, In_Hospital_Dead, Death_06)
tbi_demo <- tbi_subjects %>%
  select(Subject.ID, Study, GCS, GCS_Group, Group, Cohort, Age, Age_40_Recode,
         Age_65_Recode, Sex, GOS_06, GOS_03, GOS_12, GOS_24, GOSE_01,
         Death_03, Death_06, Death_12, Death_24, Note_from_Email)

sci_csf_demo <- sci_csf_long %>%
  left_join(sci_subjects, by = "Subject.ID", suffix = c("", ".y")) %>%
  group_by(Subject.ID) %>%
  mutate(baseline_time = first(csf_Raw_Day_Value[order(day_idx)]),
         baseline_csf = first(csf_BDTau[order(day_idx)])) %>% ungroup()
sci_serum_demo <- sci_serum_long %>%
  left_join(sci_subjects, by = "Subject.ID", suffix = c("", ".y")) %>%
  group_by(Subject.ID) %>%
  mutate(baseline_time = first(serum_Raw_Day_Value[order(day_idx)]),
         baseline_serum = first(serum_BDTau[order(day_idx)])) %>% ungroup()
tbi_csf_demo <- tbi_csf_long %>% mutate(csf_Comment = NA) %>%
  left_join(tbi_subjects, by = "Subject.ID", suffix = c("", ".y")) %>%
  group_by(Subject.ID) %>%
  mutate(baseline_time = first(csf_Raw_Day_Value[order(day_idx)]),
         baseline_csf = first(csf_BDTau[order(day_idx)])) %>% ungroup()
tbi_serum_demo <- tbi_serum_long %>% mutate(serum_Comment = NA) %>%
  left_join(tbi_subjects, by = "Subject.ID", suffix = c("", ".y")) %>%
  group_by(Subject.ID) %>%
  mutate(baseline_time = first(serum_Raw_Day_Value[order(day_idx)]),
         baseline_serum = first(serum_BDTau[order(day_idx)])) %>% ungroup()

## 5. BASELINE + MATCHED + *_cleaned -----------------------------------------

baseline_rows <- function(df) df %>% group_by(Subject.ID) %>%
  slice_min(day_idx, with_ties = FALSE) %>% ungroup()

cntr_cleaned <- as.data.frame(cntr)

serum_base_all <- bind_rows(
  cntr %>% transmute(Study, Group, Age, Sex, Subject.ID, serum_BDTau,
                     serum_log2_BDTau, mean_serum_BDTau, mean_serum_log2_BDTau,
                     Cohort, Age_40_Recode, Age_65_Recode),
  baseline_rows(sci_serum_demo) %>% transmute(Study, Group, Age, Sex, Subject.ID,
                                              serum_BDTau, serum_log2_BDTau, mean_serum_BDTau,
                                              mean_serum_log2_BDTau, Cohort, Age_40_Recode, Age_65_Recode),
  baseline_rows(tbi_serum_demo) %>% transmute(Study, Group, Age, Sex, Subject.ID,
                                              serum_BDTau, serum_log2_BDTau, mean_serum_BDTau,
                                              mean_serum_log2_BDTau, Cohort, Age_40_Recode, Age_65_Recode))
serum_baseline_cleaned <- serum_base_all %>%
  mutate(Cohort = factor(Cohort, levels = c("Control","SCI","Moderate TBI","Severe TBI","Other")),
         Cohort2 = factor(case_when(Cohort == "Control" ~ "Control",
                                    Cohort == "SCI" ~ "SCI", TRUE ~ "TBI"),
                          levels = c("Control","SCI","TBI"))) %>% as.data.frame()

csf_baseline_cleaned <- bind_rows(
  baseline_rows(sci_csf_demo) %>% mutate(Cohort2 = "SCI", Subject.Visit = as.character(Subject.Visit)),
  baseline_rows(tbi_csf_demo) %>% mutate(Cohort2 = "TBI", Subject.Visit = as.character(Subject.Visit))
) %>%
  mutate(Cohort = factor(Cohort, levels = c("SCI","Moderate TBI","Severe TBI")),
         Cohort2 = factor(Cohort2, levels = c("SCI","TBI"))) %>% as.data.frame()

make_match <- function(csf_demo, serum_demo) {
  csf_b <- baseline_rows(csf_demo) %>% select(Subject.ID, mean_csf_BDTau,
                                                  mean_csf_log2_BDTau, mean_csf_Raw_Day_Value, Cohort, Group, Study, Age,
                                                  Age_40_Recode, Age_65_Recode, Sex)
  serum_b <- baseline_rows(serum_demo) %>% select(Subject.ID, mean_serum_BDTau,
                                                  mean_serum_log2_BDTau, mean_serum_Raw_Day_Value)
  inner_join(serum_b, csf_b, by = "Subject.ID") %>%
    mutate(mean_Day_Diff = round(abs(mean_serum_Raw_Day_Value - mean_csf_Raw_Day_Value), 2),
           mean_Hour_Diff = round(mean_Day_Diff * 24, 2),
           mean_BDTau_Ratio = mean_serum_BDTau / mean_csf_BDTau,
           mean_log2_BDTau_Ratio = log2(mean_BDTau_Ratio),
           Day_6hr_Unit_Bin = cut(mean_Hour_Diff, c(-Inf,6,12,18,24,Inf),
                                  labels = c("0-6 hr","6-12 hr","12-18 hr","18-24 hr",">24 hr")),
           within_12hr_filter = as.character(mean_Hour_Diff <= 12),
           Subject.Visit = paste0(Subject.ID, ".01"))
}
match_baseline_cleaned <- bind_rows(
  make_match(sci_csf_demo, sci_serum_demo) %>% mutate(Cohort2 = "SCI"),
  make_match(tbi_csf_demo, tbi_serum_demo) %>% mutate(Cohort2 = "TBI")) %>%
  mutate(Cohort = factor(Cohort, levels = c("SCI","Moderate TBI","Severe TBI")),
         Cohort2 = factor(Cohort2, levels = c("SCI","TBI"))) %>% as.data.frame()

sci_csf_demo_BL_cleaned <- as.data.frame(baseline_rows(sci_csf_demo))
sci_serum_demo_BL_cleaned <- as.data.frame(baseline_rows(sci_serum_demo))
tbi_csf_demo_BL_cleaned <- as.data.frame(baseline_rows(tbi_csf_demo))
tbi_serum_demo_BL_cleaned <- as.data.frame(baseline_rows(tbi_serum_demo))
sci_match_demo_BL_cleaned <- as.data.frame(make_match(sci_csf_demo, sci_serum_demo))
tbi_match_demo_BL_cleaned <- as.data.frame(make_match(tbi_csf_demo, tbi_serum_demo))

## 6. SANITY CHECKS + SAVE ----------------------------------------------------
obj_names <- c(
  "cntr","sci_csf","sci_csf_demo","sci_serum","sci_serum_demo","sci_demo",
  "tbi_csf","tbi_csf_demo","tbi_serum","tbi_serum_demo","tbi_demo",
  ls()[grepl("_cleaned", ls())])

cat("\n--- dimensions of generated (synthetic) objects ---\n")
print(purrr::map_dfr(obj_names, ~tibble(object = .x, nrow = nrow(get(.x)), ncol = ncol(get(.x)))))

stopifnot(
  with(sci_serum, all.equal(serum_log2_BDTau, log2(serum_BDTau))),
  with(sci_csf,   all.equal(mean_csf_log2_BDTau, log2(mean_csf_BDTau))),
  with(match_baseline_cleaned, all.equal(mean_BDTau_Ratio, mean_serum_BDTau / mean_csf_BDTau)))
cat("\nAll internal-consistency invariants hold.\n")

save(list = obj_names, file = file.path(processed_data_path, "tbi_sci_cntr_SYNTHETIC.RData"))
cat("\nSaved -> ", file.path(processed_data_path, "tbi_sci_cntr_SYNTHETIC.RData"), "\n")

## ============================================================================
## PART 3 -- VERIFY synthetics vs the SNAPSHOT (names + types). Safe no-op when
## no reals were loaded. Compares the .real_env snapshot to the current globals.
## ============================================================================

verify_against_originals <- function(real_env = .real_env) {
  nms <- ls(real_env)
  if (!length(nms)) { cat("\n[verify] no real objects were snapshotted - skipping.\n"); return(invisible()) }
  cat("\n[verify] comparing", length(nms), "synthetic tables against the real snapshot\n")
  clean <- TRUE
  for (nm in nms) {
    real <- get(nm, envir = real_env)
    if (!is.data.frame(real) || !exists(nm)) next
    syn <- get(nm)
    miss <- setdiff(names(real), names(syn))
    extra <- setdiff(names(syn), names(real))
    shared <- intersect(names(real), names(syn))
    type_diff <- shared[vapply(shared, function(c) class(real[[c]])[1] != class(syn[[c]])[1], logical(1))]
    if (length(miss) || length(extra) || length(type_diff)) {
      clean <- FALSE
      cat("\n[", nm, "]\n")
      if (length(miss))  cat("  MISSING:", paste(miss,  collapse = ", "), "\n")
      if (length(extra)) cat("  EXTRA  :", paste(extra, collapse = ", "), "\n")
      for (c in type_diff)
        cat("  TYPE   :", c, " real=", class(real[[c]])[1], " synth=", class(syn[[c]])[1], "\n")
    }
  }
  if (clean) cat("\n[verify] all snapshotted tables match on column names and types.\n")
  invisible()
}

verify_against_originals()
