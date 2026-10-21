################################################################################
##### 2000HIV: Cytokine mediation analysis
### Author(s):Haoran Peng
### Last updated: 19th October, 2026
################################################################################

## COVID group -> cytokine -> phage presence, at 24 h and 7 days.
## lm mediator + binomial glm outcome; quasi-Bayesian mediation estimates.
## Keeps the original candidate filters, covariates, seeds and FDR families.
## Adds confidence intervals and an explicit record of every failed attempt.
## Participant data and upstream tables are on the local computer, not here.
## The data-driven candidate screen remains exploratory; it is not independent
## validation or proof of a causal mechanism.

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(stringr)
library(ggplot2)
library(mediation)

root <- Sys.getenv("HIV_PROJECT_DIR", "/Users/h.peng/Desktop/Fu Group/1--Project1/9--HIV-Virus")
after_pandemic <- TRUE
res_dir <- file.path(root, "results_clean", "04_Mediation",
                     if (after_pandemic) "after_pandemic" else "all_samples")
sims <- 1000L                       # Historical setting; increase for final inference.
make_plots <- TRUE
covar <- c("AGE", "SEX_BIRTH", "BMI_BASELINE", "CENTER", "ETHNICITY", "reads_depth")
timepoints <- tibble(
  timepoint = c("24h", "7d"),
  analysis = c("result/OUTPUT_update1", "result/OUTPUT_update1/cytokine_7D_analysis"),
  cytokines = c(
    "data/Cytokine production data/cytokine_production_24H/2000HIV_EX_VIVO_24H_52meas_1760samples_afterQC_TRANSFORMED.xlsx",
    "data/Cytokine production data/cytokine_production_7D/2000HIV_EX_VIVO_7DAY_38meas_1754samples_afterQC_TRANSFORMED.xlsx"
  )
)

# --- 4.1 Candidate paths ---------------------------------------------------
read_candidates <- function(root, analysis_dir, after_pandemic) {
  cyto <- read_tsv(file.path(root, analysis_dir, "maaslin3_phage_cytokine_significant_summary.tsv"),
                    show_col_types = FALSE)
  pheno <- read_tsv(file.path(root, "result/OUTPUT_update/maaslin3_phage_phenotype_significant_summary.tsv"),
                     show_col_types = FALSE)
  cyto_pheno <- read_tsv(file.path(root, analysis_dir, "phenotype_cytokine_association/significant_results.tsv"),
                          show_col_types = FALSE)
  phage_glm <- read_tsv(file.path(root, analysis_dir, "phage_presence_cytokine_phenotype_association/all_results.tsv"),
                         show_col_types = FALSE)
  label_term <- function(df) df %>%
    mutate(model = str_to_lower(model), term = ifelse(!is.na(name) & name != "", name, metadata))
  candidates <- label_term(cyto) %>%
    filter(model == "prevalence", is.finite(qval_individual), qval_individual < 0.1) %>%
    dplyr::select(phage = feature, cytokine = term, beta_cytokine = coef, qval_cytokine = qval_individual) %>%
    inner_join(label_term(pheno) %>%
      filter(model == "prevalence", is.finite(qval_individual), qval_individual < 0.1) %>%
      dplyr::select(phage = feature, phenotype = term, beta_phenotype = coef, qval_phenotype = qval_individual),
      by = "phage", relationship = "many-to-many") %>%
    inner_join(cyto_pheno %>% rename(beta_cytokine_phenotype = beta,
      pval_cytokine_phenotype = pval, fdr_cytokine_phenotype = qval),
      by = c("cytokine", "phenotype")) %>%
    inner_join(phage_glm %>% filter(is.finite(cyto_qval), cyto_qval < 0.1) %>%
      dplyr::select(phage, cytokine, phenotype, n_phage_model = n,
             beta_phage_model_cytokine = cyto_beta, pval_phage_model_cytokine = cyto_pval,
             qval_phage_model_cytokine = cyto_qval, qval_phage_model_phenotype = pheno_qval),
      by = c("phage", "cytokine", "phenotype")) %>%
    filter(is.finite(fdr_cytokine_phenotype), fdr_cytokine_phenotype < 0.1) %>% distinct()
  if (after_pandemic) candidates <- filter(candidates, str_detect(phenotype, "^COVID_GROUP[0-9]+$"))
  candidates
}

# --- 4.2 Shared participant, phage and bacterial inputs --------------------
host_label <- function(x) {
  x[is.na(x) | x == ""] <- "Unclassified"
  str_replace_all(paste0("g__", x), "[^A-Za-z0-9_]", "_")
}

read_shared_data <- function(root, candidates, after_pandemic) {
  pheno_env <- new.env()
  load(file.path(root, "data/2000hivviromerelevantinformationdata/2000HIV_updated_metadata.RData"),
       envir = pheno_env)
  metadata <- as_tibble(pheno_env$pheno) %>% mutate(Sample_ID = as.character(SampleID))
  names_map <- as_tibble(readRDS(file.path(root, "data/name_trans.rds"))) %>%
    transmute(reads_sample = as.character(sample_name2), Sample_ID = as.character(sample_name))
  reads <- read_tsv(file.path(root, "data/2000hivviromerelevantinformationdata/2000HIV_reads_final.tsv"),
                     show_col_types = FALSE) %>%
    rename(reads_sample = sample) %>% mutate(reads_sample = as.character(reads_sample)) %>%
    left_join(names_map, by = "reads_sample") %>%
    transmute(Sample_ID = coalesce(Sample_ID, reads_sample), reads_depth = as.numeric(reads_depth))
  exclusion <- read_excel(file.path(root, "data/2000hivviromerelevantinformationdata/20220322_Exclusion_2000HIV.xlsx"),
                           sheet = "Sum")[[1]]
  stopifnot(!anyDuplicated(metadata$Sample_ID), !anyDuplicated(reads$Sample_ID))

  path <- file.path(root, "result/HIV_Abundance_Table_01052025.txt")
  ids <- str_split(read_lines(path, n_max = 1), "\t")[[1]]
  raw <- read_tsv(path, skip = 1, col_names = c("Virus_ID", ids),
                   col_types = cols(Virus_ID = col_character(), .default = col_double()))
  phage_map <- candidates %>% distinct(phage) %>%
    mutate(Virus_ID = str_split_fixed(phage, fixed(" | "), 3)[, 1])
  stopifnot(!anyDuplicated(raw$Virus_ID), all(phage_map$Virus_ID %in% raw$Virus_ID))
  presence <- t(as.matrix(raw[match(phage_map$Virus_ID, raw$Virus_ID), -1] > 0))
  colnames(presence) <- phage_map$phage
  storage.mode(presence) <- "integer"
  metadata <- metadata %>% filter(Sample_ID %in% ids, !Sample_ID %in% exclusion,
                                  is.na(AB_1M) | AB_1M == 0) %>%
    dplyr::select(-any_of("reads_depth")) %>% left_join(reads, by = "Sample_ID") %>%
    mutate(AGE = as.numeric(as.character(AGE)), BMI_BASELINE = as.numeric(as.character(BMI_BASELINE)),
           SEX_BIRTH = factor(SEX_BIRTH), CENTER = factor(CENTER), ETHNICITY = factor(ETHNICITY),
           PANDEMIC_BEFOREAFTER = as.numeric(as.character(PANDEMIC_BEFOREAFTER)))
  if (after_pandemic) metadata <- filter(metadata, PANDEMIC_BEFOREAFTER == 1)

  species <- as.matrix(readRDS(file.path(root, "data/sp_re_metaphlan.rds")))
  taxonomy <- read_tsv(file.path(root, "data/mpa_vOct22_CHOCOPhlAnSGB_202212_SGB2GTDB.tsv"),
                        col_names = c("SGB", "GTDB"), show_col_types = FALSE)
  stopifnot(!anyDuplicated(taxonomy$SGB), !anyDuplicated(rownames(species)))
  genus <- str_extract(taxonomy$GTDB, "(?<=g__)[^;]*")
  genus <- str_squish(str_replace_all(genus, "_", " "))
  genus <- host_label(genus[match(sub("^t__", "", colnames(species)), taxonomy$SGB)])
  host <- t(rowsum(t(species), group = genus, na.rm = TRUE))
  annotation <- read_tsv(file.path(root, "result/HIV_Viral_Metadata_06052025.txt"),
                          col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
    transmute(Virus_ID, host_column = host_label(Bacterial_genus_host))
  stopifnot(!anyDuplicated(annotation$Virus_ID))
  list(metadata = metadata, presence = presence, host = host, annotation = annotation)
}

# --- 4.3 Phenotype coding --------------------------------------------------
phenotype_series <- function(metadata, name, after_pandemic) {
  if (name %in% names(metadata)) return(list(values = metadata[[name]], base = name))
  choices <- names(metadata)[startsWith(str_to_lower(name), str_to_lower(names(metadata)))]
  if (!length(choices)) stop("Unknown phenotype: ", name)
  base <- choices[which.max(nchar(choices))]
  level <- str_remove(substr(name, nchar(base) + 1, nchar(name)), "^_+")
  values <- metadata[[base]]
  indicator <- as.character(values) == level
  if (!sum(indicator, na.rm = TRUE) && !is.na(suppressWarnings(as.numeric(level)))) {
    indicator <- as.character(values) == as.character(as.numeric(level))
  }
  values <- ifelse(is.na(values), NA_integer_, as.integer(indicator))
  if (after_pandemic && str_detect(name, "^COVID_GROUP[0-9]+$")) {
    group <- as.character(metadata$COVID_GROUP)
    values <- ifelse(group %in% c("1", level), as.integer(group == level), NA_integer_)
  }
  list(values = values, base = base)
}

estimate_mediation <- function(mediator_fit, outcome_fit, sims) {
  mediation::mediate(mediator_fit, outcome_fit, treat = "phenotype", mediator = "cytokine",
                     sims = sims, boot = FALSE)
}

# --- 4.4 One path ----------------------------------------------------------
fit_mediation <- function(candidate, data, cytokines, add_host, covar, sims, after_pandemic) {
  phage <- candidate$phage[[1]]
  cytokine <- candidate$cytokine[[1]]
  phenotype <- candidate$phenotype[[1]]
  row <- tibble(phage = phage, cytokine = cytokine, phenotype = phenotype,
                model_set = if (add_host) "host_adjusted" else "standard",
                n = 0L, n_present = NA_integer_, n_absent = NA_integer_,
                status = "failed", message = NA_character_, warnings = NA_character_,
                host_mapping = if (add_host) "unmapped" else "not_requested",
                mediator_rank = NA_integer_, outcome_rank = NA_integer_,
                glm_converged = NA, acme_avg = NA_real_, acme_p = NA_real_,
                acme_low = NA_real_, acme_high = NA_real_, ade_avg = NA_real_, ade_p = NA_real_,
                ade_low = NA_real_, ade_high = NA_real_, total_effect = NA_real_, total_p = NA_real_,
                total_low = NA_real_, total_high = NA_real_, prop_med = NA_real_, prop_med_p = NA_real_,
                prop_med_low = NA_real_, prop_med_high = NA_real_)
  warning_text <- character()
  tryCatch(withCallingHandlers({
    if (!phage %in% colnames(data$presence) || !cytokine %in% names(cytokines)) stop("Missing phage or cytokine")
    info <- phenotype_series(data$metadata, phenotype, after_pandemic)
    covariates <- if (startsWith(info$base, "MSM_derived")) setdiff(covar, "SEX_BIRTH") else covar
    df <- as.data.frame(data$metadata[, c("Sample_ID", covariates)])
    df$phenotype <- info$values
    df$cytokine <- cytokines[[cytokine]][match(df$Sample_ID, cytokines$Sample_ID)]
    df$phage <- data$presence[match(df$Sample_ID, rownames(data$presence)), phage]
    if (add_host) {
      virus <- str_split_fixed(phage, fixed(" | "), 3)[, 1]
      column <- data$annotation$host_column[match(virus, data$annotation$Virus_ID)]
      mapped <- !is.na(column) && column %in% colnames(data$host)
      row$host_mapping <- if (mapped) "mapped" else "unmapped"
      # Historical behaviour retained: unknown host/missing bacterial profile -> 0.
      # Such cases can be singular; their failure is now recorded instead of lost.
      df$host_genus_abundance <- if (mapped)
        replace_na(data$host[match(df$Sample_ID, rownames(data$host)), column], 0) else 0
      covariates <- c(covariates, "host_genus_abundance")
    }
    df <- df[complete.cases(df[, c("phage", "cytokine", "phenotype", covariates)]), , drop = FALSE]
    row$n <- nrow(df)
    row$n_present <- sum(df$phage == 1)
    row$n_absent <- sum(df$phage == 0)
    if (nrow(df) < 10 || any(vapply(df[c("phage", "cytokine", "phenotype")], n_distinct, integer(1)) < 2)) {
      stop("Too few samples or a constant treatment/mediator/outcome")
    }
    med_fit <- lm(reformulate(c("phenotype", covariates), response = "cytokine"), data = df)
    out_fit <- glm(reformulate(c("phenotype", "cytokine", covariates), response = "phage"),
                    data = df, family = binomial())
    row$mediator_rank <- med_fit$rank
    row$outcome_rank <- out_fit$rank
    row$glm_converged <- out_fit$converged
    effect <- estimate_mediation(med_fit, out_fit, sims)
    row$acme_avg <- effect$d.avg
    row$acme_p <- effect$d.avg.p
    row$ade_avg <- effect$z.avg
    row$ade_p <- effect$z.avg.p
    row$total_effect <- effect$tau.coef
    row$total_p <- effect$tau.p
    row$prop_med <- effect$n.avg
    row$prop_med_p <- effect$n.avg.p
    row$acme_low <- effect$d.avg.ci[1]
    row$acme_high <- effect$d.avg.ci[2]
    row$ade_low <- effect$z.avg.ci[1]
    row$ade_high <- effect$z.avg.ci[2]
    row$total_low <- effect$tau.ci[1]
    row$total_high <- effect$tau.ci[2]
    row$prop_med_low <- effect$n.avg.ci[1]
    row$prop_med_high <- effect$n.avg.ci[2]
    row$status <- "ok"
    row$message <- NA_character_
  }, warning = function(w) {
    warning_text <<- c(warning_text, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) { row$message <<- conditionMessage(e) })
  row$warnings <- if (length(warning_text)) paste(unique(warning_text), collapse = "; ") else NA_character_
  row
}

adjust_mediation <- function(results) {
  # Same historical family: successful paths, separately per timepoint/model set.
  # Retaining failed rows does not silently redefine the original BH denominator.
  results$acme_fdr <- results$total_fdr <- NA_real_
  for (model in unique(results$model_set)) {
    idx <- which(results$model_set == model & results$status == "ok")
    results$acme_fdr[idx] <- p.adjust(results$acme_p[idx], "BH")
    results$total_fdr[idx] <- p.adjust(results$total_p[idx], "BH")
  }
  results
}

# --- 4.5 Run 24 h and 7 days ----------------------------------------------
candidate_tables <- lapply(seq_len(nrow(timepoints)), function(i) {
  read_candidates(root, timepoints$analysis[i], after_pandemic)
})
if (!sum(vapply(candidate_tables, nrow, integer(1)))) stop("No candidate paths passed the upstream filters")
shared <- read_shared_data(root, bind_rows(candidate_tables), after_pandemic)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
all_results <- list()
for (i in seq_len(nrow(timepoints))) {
  timepoint <- timepoints$timepoint[i]
  candidates <- candidate_tables[[i]]
  out <- file.path(res_dir, timepoint)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_tsv(candidates, file.path(out, "candidate_paths.tsv"))
  if (!nrow(candidates)) next
  cytokines <- read_xlsx(file.path(root, timepoints$cytokines[i])) %>%
    rename(Sample_ID = ID) %>% mutate(Sample_ID = as.character(Sample_ID))
  stopifnot(!anyDuplicated(cytokines$Sample_ID))
  sets <- lapply(c(FALSE, TRUE), function(add_host) {
    set.seed(if (add_host) 20260911L else 20260910L)
    bind_rows(lapply(seq_len(nrow(candidates)), function(j) {
      fit_mediation(candidates[j, ], shared, cytokines, add_host, covar, sims, after_pandemic)
    }))
  })
  results <- adjust_mediation(bind_rows(sets)) %>% mutate(timepoint = timepoint, sims = sims)
  all_results[[timepoint]] <- results
  write_tsv(results, file.path(out, "all_attempts.tsv"))
  write_tsv(filter(results, status != "ok"), file.path(out, "failed_paths.tsv"))
  standard <- filter(results, model_set == "standard", status == "ok")
  host <- filter(results, model_set == "host_adjusted", status == "ok")
  write_tsv(standard, file.path(out, "mediation_results.tsv"))
  # This legacy OR table may include total-effect-only paths; it is not an ACME-only table.
  write_tsv(filter(standard, acme_fdr < 0.1 | total_fdr < 0.1),
            file.path(out, "mediation_results_fdr_qfilt.tsv"))
  write_tsv(host, file.path(out, "mediation_results_host_genus_cov_gtdb.tsv"))
  write_tsv(filter(host, acme_fdr < 0.05),
            file.path(out, "mediation_results_host_genus_cov_gtdb_acme_fdr_lt_0.05.tsv"))
}
results <- bind_rows(all_results)
write_tsv(results, file.path(res_dir, "all_timepoints.tsv"))

# --- 4.6 Simple overview plots --------------------------------------------
# These replace the result-specific 9-path layout with a data-driven overview.
# Separate single-mediator effects and proportions must not be added together.
plot_df <- results %>%
  filter(model_set == "host_adjusted", status == "ok", acme_fdr < 0.05,
         !str_detect(cytokine, regex("RPMI", ignore_case = TRUE))) %>%
  mutate(path = paste(phenotype, cytokine, phage, sep = " / "),
         Cytokine = paste(cytokine, timepoint, sep = "_"))
if (make_plots && nrow(plot_df)) {
  p_effect <- ggplot(plot_df, aes(x = acme_avg, y = reorder(path, acme_avg))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_segment(aes(x = acme_low, xend = acme_high, yend = reorder(path, acme_avg))) +
    geom_point(size = 2) + facet_wrap(~ timepoint, scales = "free_y") +
    labs(x = "ACME (difference in predicted phage detection probability)", y = NULL) + theme_bw()
  ggsave(file.path(res_dir, "mediation_effects.pdf"), p_effect, width = 12,
          height = max(4, nrow(plot_df) * 0.3))
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p_flow <- ggplot(plot_df, aes(axis1 = phenotype, axis2 = Cytokine, axis3 = phage,
                                 y = abs(acme_avg))) +
      ggalluvial::geom_alluvium(aes(fill = phage), alpha = 0.5) +
      ggalluvial::geom_stratum(fill = "grey95", colour = "white") +
      ggalluvial::stat_stratum(aes(label = after_stat(stratum)), geom = "text", size = 3) +
      scale_x_discrete(limits = c("COVID group", "Cytokine", "Phage")) +
      labs(x = NULL, y = "Sum of |ACME|", caption = "Separate single-mediator models; flow widths are not additive causal contributions.") +
      theme_bw() + theme(legend.position = "none", panel.grid = element_blank())
    ggsave(file.path(res_dir, "mediation_sankey.pdf"), p_flow, width = 12, height = 8)
  }
}
capture.output(sessionInfo(), file = file.path(res_dir, "sessionInfo.txt"))
