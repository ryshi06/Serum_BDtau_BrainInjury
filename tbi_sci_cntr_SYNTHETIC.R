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
  n_mild_tbi = 40,
  tbi_severe_frac = 0.72,
  
  visits_sci_csf = function(n) pmax(1L, rpois(n, 4)),
  visits_sci_serum = function(n) pmax(1L, rpois(n, 7)),
  visits_tbi_csf = function(n) pmax(1L, rpois(n, 12)),
  visits_tbi_serum = function(n) pmax(1L, rpois(n, 6)),
  
  bdtau = list(
    Control = list(serum = c(med =   6, sd = 0.5)),
    SCI = list(serum = c(med =  40, sd = 0.8), csf = c(med = 1100, sd = 0.9)),
    `Mild TBI` = list(serum = c(med =  35, sd = 0.7)),
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
asia_levels <- c("A","B","C","D", NA)

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
n_mild <- cfg$n_mild_tbi
n_non_mild <- cfg$n_tbi - n_mild
n_sev <- round(n_non_mild * cfg$tbi_severe_frac)
n_mod <- n_non_mild - n_sev
tbi_subjects <- tibble(
  Subject.ID = as.character(1000:(1000 + cfg$n_tbi - 1)),
  Study = "BTRC",
  Group  = c(rep("Severe", n_sev), rep("Moderate", n_mod), rep("Mild", n_mild)),
  Cohort = c(rep("Severe TBI", n_sev), rep("Moderate TBI", n_mod), rep("Mild TBI", n_mild)),
  Age = round(rnorm(cfg$n_tbi, 42, 18)) %>% pmin(90) %>% pmax(16),
  Sex = sample(c("Male","Female"), cfg$n_tbi, TRUE, c(.72,.28)),
  GCS = c(sample(3:8, n_sev, TRUE), sample(9:12, n_mod, TRUE), sample(13:15, n_mild, TRUE)),
  GOSE_01 = NA, Note_from_Email = NA
) %>%
  mutate(GCS_Group = as.character(gcs_group(GCS)),
         Age_40_Recode = recode_age(Age, 40),
         Age_65_Recode = recode_age(Age, 65),
         .unfavorable = sample(c(TRUE, FALSE), cfg$n_tbi, TRUE, c(.35,.65)))

for (col in c("GOS_03","GOS_06","GOS_12","GOS_24"))
  tbi_subjects[[col]] <- ifelse(tbi_subjects$.unfavorable, "Unfavorable",
                                sample(c("Favorable", NA), cfg$n_tbi, TRUE, c(.85,.15)))
for (col in c("Death_03","Death_06","Death_12","Death_24"))
  tbi_subjects[[col]] <- ifelse(tbi_subjects$.unfavorable & runif(cfg$n_tbi) < .35, "Dead",
                                sample(c("Survived", NA), cfg$n_tbi, TRUE, c(.85,.15)))
tbi_subjects$.unfavorable <- NULL

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
fmt_tbi <- function(id, k) sprintf("%s.%02d", id, k)

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
  mutate(Cohort = factor(Cohort, levels = c("Control","SCI","Mild TBI","Moderate TBI","Severe TBI","Other")),
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

## 6. SYNTHETIC N4PD COMPARISON DATA ------------------------------------------

make_n4pd_rows <- function(df, id_col = "Subject.ID", specimen_col = "Subject.Visit") {
  cohort_multiplier <- case_when(
    df$Cohort == "Control" ~ 0.55,
    df$Cohort == "SCI" ~ 0.85,
    df$Cohort == "Mild TBI" ~ 0.95,
    df$Cohort == "Moderate TBI" ~ 1.15,
    df$Cohort == "Severe TBI" ~ 1.45,
    TRUE ~ 1
  )
  tibble(
    !!id_col := as.character(df$Subject.ID),
    !!specimen_col := df$Subject.Visit,
    `Time From Injury` = suppressWarnings(as.numeric(df$serum_Raw_Day_Value * 24)),
    `BD-tau (pg/ml)` = round(df$serum_BDTau * rlnorm(nrow(df), 0, 0.08), 3),
    `BD-tau_N4PD` = round(df$serum_BDTau * rlnorm(nrow(df), log(1.05), 0.12), 3),
    GFAP_N4PD = round(40 * cohort_multiplier * rlnorm(nrow(df), 0, 0.7), 3),
    NfL_N4PD = round(18 * cohort_multiplier * rlnorm(nrow(df), 0, 0.65), 3),
    `UCH-L1_N4PD` = round(28 * cohort_multiplier * rlnorm(nrow(df), 0, 0.75), 3)
  )
}

write_synthetic_n4pd <- function(path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    warning("Package 'openxlsx' is not installed; skipping synthetic N4PD workbook.")
    return(invisible(FALSE))
  }
  
  severe_sheet <- tbi_serum_demo %>%
    filter(Cohort == "Severe TBI") %>%
    make_n4pd_rows(id_col = "Subject ID", specimen_col = "Specimen ID")
  
  sci_sheet <- sci_serum_demo %>%
    make_n4pd_rows(id_col = "Participant ID", specimen_col = "Specimen ID")
  
  mild_control_sheet <- bind_rows(
    tbi_serum_demo %>%
      filter(Cohort == "Mild TBI") %>%
      transmute(Subject.ID, Subject.Visit = as.character(Subject.Visit), Cohort,
                serum_Raw_Day_Value, serum_BDTau),
    cntr_cleaned %>%
      transmute(Subject.ID, Subject.Visit = as.character(Subject.Visit), Cohort,
                serum_Raw_Day_Value, serum_BDTau)
  ) %>%
    make_n4pd_rows(id_col = "Subject ID", specimen_col = "Current Label") %>%
    mutate(
      Group = case_when(
        `Subject ID` %in% cntr_cleaned$Subject.ID ~ "Controls",
        TRUE ~ "Mild"
      ),
      `BD-tau_June2025` = `BD-tau (pg/ml)`
    ) %>%
    select(`Subject ID`, Group, `Current Label`, `BD-tau_June2025`,
           `BD-tau_N4PD`, GFAP_N4PD, NfL_N4PD, `UCH-L1_N4PD`)
  
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(
    list(
      Severe_TBI = severe_sheet,
      SCI = sci_sheet,
      Mild_TBI_Control = mild_control_sheet
    ),
    file = path,
    overwrite = TRUE
  )
  invisible(TRUE)
}

synthetic_n4pd_path <- file.path("~/Desktop/Projects/2025_05_19_TBI_SCI_BDTau/Raw_Data/N4PD/TBI SCI_N4PD_Results_11212025.xlsx")

make_n4pd_direct <- function(df) {
  make_n4pd_rows(df) %>%
    transmute(
      Subject.ID,
      Subject.Visit = as.character(Subject.Visit),
      `BD-tau (pg/ml)`,
      n4pd_BDtau_conc = `BD-tau_N4PD`,
      GFAP_N4PD,
      NfL_N4PD,
      `UCH-L1_N4PD`
    )
}

commercial_bdtau <- bind_rows(
  tbi_serum_demo %>%
    dplyr::select(Subject.ID, Subject.Visit, Cohort, Age, Sex, serum_BDTau, serum_Raw_Day_Value),
  sci_serum_demo %>%
    dplyr::select(Subject.ID, Subject.Visit, Cohort, Age, Sex, serum_BDTau, serum_Raw_Day_Value),
  cntr_cleaned %>%
    dplyr::select(Subject.ID, Subject.Visit, Cohort, Age, Sex, serum_BDTau, serum_Raw_Day_Value)
) %>%
  mutate(
    Subject.ID = as.character(Subject.ID),
    Subject.Visit = as.character(Subject.Visit)
  ) %>%
  filter(!is.na(serum_BDTau)) %>%
  group_by(Subject.ID, Subject.Visit) %>%
  arrange(desc(!is.na(serum_Raw_Day_Value)), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  distinct()

n4pd_bdtau <- bind_rows(
  tbi_serum_demo %>%
    filter(Cohort %in% c("Mild TBI", "Moderate TBI", "Severe TBI")) %>%
    make_n4pd_direct(),
  sci_serum_demo %>%
    make_n4pd_direct(),
  cntr_cleaned %>%
    make_n4pd_direct()
) %>%
  distinct()

bdtau_conc <- commercial_bdtau %>%
  filter(!is.na(serum_BDTau)) %>%
  dplyr::select(Subject.ID, Subject.Visit, Cohort, Sex, Age, serum_BDTau) %>%
  inner_join(n4pd_bdtau, by = c("Subject.ID", "Subject.Visit")) %>%
  distinct()

bdtau_demo <- bdtau_conc %>%
  group_by(Subject.ID) %>%
  summarise(
    Cohort = first(na.omit(Cohort)),
    Age = first(na.omit(Age)),
    Sex = first(na.omit(Sex)),
    .groups = "drop"
  ) %>%
  distinct() %>%
  filter(!is.na(Sex))

bdtau <- bdtau_conc %>%
  dplyr::select(-Cohort, -Age, -Sex) %>%
  left_join(bdtau_demo, by = "Subject.ID") %>%
  filter(!is.na(Sex)) %>%
  mutate(commercial_BDtau_conc = serum_BDTau) %>%
  dplyr::select(Subject.ID, Subject.Visit, Cohort, Sex, Age,
                commercial_BDtau_conc, n4pd_BDtau_conc, GFAP_N4PD, NfL_N4PD, `UCH-L1_N4PD`) %>%
  mutate(
    across(where(is.character), stringr::str_squish),
    across(where(is.numeric), ~ round(.x, 6))
  ) %>%
  distinct()

dup_ids <- bdtau$Subject.Visit[duplicated(bdtau$Subject.Visit)]

bdtau_clean <- bdtau %>%
  filter(!(Subject.Visit %in% dup_ids)) %>%
  filter(!is.na(n4pd_BDtau_conc) & !is.na(`UCH-L1_N4PD`))

n_visits <- bdtau_clean %>%
  group_by(Subject.ID) %>%
  summarise(n_visits = n_distinct(Subject.Visit), .groups = "drop")

asia_a_ids <- sci_serum_demo_BL_cleaned[
  !is.na(sci_serum_demo_BL_cleaned$ASIA) & sci_serum_demo_BL_cleaned$ASIA == "A", ]$Subject.ID
asia_b_ids <- sci_serum_demo_BL_cleaned[
  !is.na(sci_serum_demo_BL_cleaned$ASIA) & sci_serum_demo_BL_cleaned$ASIA == "B", ]$Subject.ID
asia_cd_ids <- sci_serum_demo_BL_cleaned[
  !is.na(sci_serum_demo_BL_cleaned$ASIA) & sci_serum_demo_BL_cleaned$ASIA %in% c("C","D"), ]$Subject.ID

bdtau_conc_baseline <- bdtau_clean %>%
  filter(!is.na(commercial_BDtau_conc), !is.na(n4pd_BDtau_conc)) %>%
  group_by(Subject.ID) %>%
  arrange(stringr::str_order(Subject.Visit, numeric = TRUE), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(n_visits, by = "Subject.ID") %>%
  mutate(
    Cohort2 = factor(case_when(
      Cohort %in% c("Mild TBI", "Moderate TBI", "Severe TBI", "TBI") ~ "TBI",
      Cohort == "SCI" ~ "SCI",
      Cohort == "Control" ~ "Control",
      TRUE ~ NA_character_
    ), levels = c("Control", "SCI", "TBI")),
    n4pd_BDtau_conc = suppressWarnings(as.numeric(n4pd_BDtau_conc)),
    GFAP_N4PD = suppressWarnings(as.numeric(GFAP_N4PD)),
    NfL_N4PD = suppressWarnings(as.numeric(NfL_N4PD)),
    `UCH-L1_N4PD` = suppressWarnings(as.numeric(`UCH-L1_N4PD`)),
    Cohort3 = case_when(
      Subject.ID %in% asia_a_ids ~ "A SCI",
      Subject.ID %in% asia_b_ids ~ "B SCI",
      Subject.ID %in% asia_cd_ids ~ "C/D SCI",
      TRUE ~ Cohort
    )
  )

bdtau_conc_baseline_roc <- bdtau_conc_baseline %>%
  rename(UCH_L1_N4PD = `UCH-L1_N4PD`)

## 7. SANITY CHECKS + SAVE ----------------------------------------------------
obj_names <- c(
  "cntr","sci_csf","sci_csf_demo","sci_serum","sci_serum_demo","sci_demo",
  "tbi_csf","tbi_csf_demo","tbi_serum","tbi_serum_demo","tbi_demo",
  "commercial_bdtau", "n4pd_bdtau", "bdtau_conc", "bdtau_demo", "bdtau",
  "bdtau_clean", "n_visits", "bdtau_conc_baseline", "bdtau_conc_baseline_roc",
  ls()[grepl("_cleaned", ls())])

cat("\n--- dimensions of generated (synthetic) objects ---\n")
print(purrr::map_dfr(obj_names, ~tibble(object = .x, nrow = nrow(get(.x)), ncol = ncol(get(.x)))))

stopifnot(
  with(sci_serum, all.equal(serum_log2_BDTau, log2(serum_BDTau))),
  with(sci_csf,   all.equal(mean_csf_log2_BDTau, log2(mean_csf_BDTau))),
  with(match_baseline_cleaned, all.equal(mean_BDTau_Ratio, mean_serum_BDTau / mean_csf_BDTau)),
  all(c("Control", "SCI", "Mild TBI", "Moderate TBI", "Severe TBI") %in%
        unique(as.character(bdtau_conc_baseline$Cohort))))
cat("\nAll internal-consistency invariants hold.\n")

save(list = obj_names, file = file.path(processed_data_path, "tbi_sci_cntr_SYNTHETIC.RData"))
cat("\nSaved -> ", file.path(processed_data_path, "tbi_sci_cntr_SYNTHETIC.RData"), "\n")

if (write_synthetic_n4pd(synthetic_n4pd_path))
  cat("Saved -> ", synthetic_n4pd_path, "\n")

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
