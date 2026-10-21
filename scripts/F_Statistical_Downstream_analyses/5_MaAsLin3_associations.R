################################################################################
##### 2000HIV: MaAsLin3 associations
### Author(s):Haoran Peng
### Last updated: 19th October, 2026
################################################################################

## Reuses each version's prepared abundance matrices and model-specific cases.
## Original/basic6: NONE + LOG; primary target-global BH; host SGB CLR terms.
## Clinical: TSS + LOG; primary native q; one host-genus CLR term.
## Host adjustment is a sensitivity analysis, not mediation.
## Historical clinical reference levels and redundant sex/MSM coding are kept.

root <- Sys.getenv("HIV_PROJECT_DIR", "/scratch/p312334/project/16--HIV_Immun")
variant <- "MaAsLin3_clinical/base7"  # Or MaAsLin3, MaAsLin3_basic6,
                                    # MaAsLin3_clinical/base7_plus_partners.
input_dir <- file.path(root, "Result_V1.0", variant)
res_dir <- file.path(root, "results_clean", "03_MaAsLin3", variant)
refit <- FALSE                      # TRUE reruns full baseline and host models.
include_host <- TRUE
q_cutoff <- 0.1
clinical <- startsWith(variant, "MaAsLin3_clinical/")
normalization <- if (clinical) "TSS" else "NONE"
seed <- if (clinical) 1L else 20260908L
.libPaths(c(file.path(root, "Result_V1.0", "MaAsLin3", "R_library"), .libPaths()))

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(maaslin3)

# --- 3.1 Prepared inputs ---------------------------------------------------
# Reading the selected version's own panel avoids the old 36/37-target drift.
panel <- read_tsv(file.path(input_dir, "tables", "phenotype_panel.tsv"), show_col_types = FALSE)
inputs <- readRDS(file.path(input_dir, "cache", "model_inputs.rds"))
communities <- c("bacteria_species", "phage", "phage_host_genus")
suffix <- if (clinical) "_input_abundance.rds" else "_relative_abundance.rds"
abundance <- setNames(lapply(communities, function(community) {
  readRDS(file.path(input_dir, "cache", paste0(community, suffix)))
}), communities)
host_map <- read_tsv(file.path(input_dir, "tables", "votu_to_host_genus.tsv"),
                      show_col_types = FALSE) %>% dplyr::select(feature, host_genus)
stopifnot(!anyDuplicated(panel$phenotype), !anyDuplicated(host_map$feature))
get_metadata <- function(community, phenotype) {
  if (clinical) inputs[[community]][[phenotype]]$metadata else inputs[[phenotype]]$metadata
}

# --- 3.2 Fit one complete response family ----------------------------------
# Never fit only significant features: this would change the abundance median null.
fit_maaslin <- function(x, md, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  maaslin3::maaslin3(
    input_data = as.data.frame(x, check.names = FALSE), input_metadata = md,
    output = out, fixed_effects = names(md), normalization = normalization,
    transform = "LOG", min_prevalence = 0, min_abundance = 0, min_variance = 0,
    standardize = TRUE, augment = TRUE, median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE, warn_prevalence = TRUE,
    max_significance = q_cutoff, correction = "BH", cores = 1,
    plot_summary_plot = FALSE, plot_associations = FALSE,
    save_models = FALSE, verbosity = "WARN"
  )
  read_tsv(file.path(out, "all_results.tsv"),
           col_types = cols(value = col_character(), error = col_character()), show_col_types = FALSE)
}

target_support <- function(z, x, md, phenotype) {
  z <- z %>% filter(metadata == phenotype) %>% mutate(value = as.character(value))
  z$n_complete <- nrow(md)
  z$reference <- if (is.factor(md[[phenotype]])) levels(md[[phenotype]])[1] else NA_character_
  z$n_reference_present <- z$n_reference_absent <- NA_integer_
  z$n_level_present <- z$n_level_absent <- NA_integer_
  if (is.factor(md[[phenotype]])) {
    for (level in unique(z$value)) {
      idx <- which(z$value == level)
      at_ref <- md[[phenotype]] == z$reference[1]
      at_level <- md[[phenotype]] == level
      nr <- colSums(x[at_ref, z$feature[idx], drop = FALSE] > 0)
      nl <- colSums(x[at_level, z$feature[idx], drop = FALSE] > 0)
      z$n_reference_present[idx] <- nr
      z$n_reference_absent[idx] <- sum(at_ref) - nr
      z$n_level_present[idx] <- nl
      z$n_level_absent[idx] <- sum(at_level) - nl
    }
    z$support_ok <- with(z, n_reference_present >= 10 & n_level_present >= 10 &
      (model == "abundance" | (n_reference_absent >= 10 & n_level_absent >= 10)))
  } else {
    z$support_ok <- with(z, N_not_zero >= 10 & (model == "abundance" | N - N_not_zero >= 10))
  }
  z$observations_per_design_term <- ifelse(z$model == "abundance", z$N_not_zero, z$N) /
    ncol(model.matrix(reformulate(names(md)), md))
  z
}

adjust_targets <- function(z, clinical, cutoff) {
  z %>%
    mutate(valid_test = is.na(error) & is.finite(pval_individual)) %>%
    group_by(composition) %>%
    mutate(q_target_global = {
      q <- rep(NA_real_, n())
      q[valid_test] <- p.adjust(pval_individual[valid_test], "BH")
      q
    }) %>%
    ungroup() %>%
    mutate(q_primary = if (clinical) qval_individual else q_target_global,
           effect_vs_null = coef - null_hypothesis,
           significant_q10 = valid_test & !is.na(q_primary) & q_primary < cutoff,
           qc_supported = valid_test & support_ok & abs(coef) <= 15,
           supported_q10 = significant_q10 & qc_supported)
}

# --- 3.3 Baseline associations --------------------------------------------
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
tables <- designs <- samples <- list()
for (community in communities) {
  for (phenotype in panel$phenotype) {
    message(community, " / ", phenotype)
    md <- get_metadata(community, phenotype)
    stopifnot(!is.null(md), !anyNA(md), !anyDuplicated(rownames(md)))
    x <- abundance[[community]][rownames(md), , drop = FALSE]
    if (refit) {
      out <- file.path(res_dir, "models", community, phenotype)
      z <- target_support(fit_maaslin(x, md, out), x, md, phenotype)
      write_tsv(z, file.path(out, "target_results.tsv"))
    } else {
      z <- read_tsv(file.path(input_dir, "models", community, phenotype, "target_results.tsv"),
                     col_types = cols(value = col_character(), reference = col_character(),
                                      error = col_character()), show_col_types = FALSE)
    }
    X <- model.matrix(reformulate(names(md)), md)
    key <- paste(community, phenotype, sep = "/")
    tables[[key]] <- z %>% mutate(composition = community, phenotype = phenotype,
                                   category = panel$category[match(phenotype, panel$phenotype)])
    designs[[key]] <- tibble(composition = community, phenotype = phenotype, n = nrow(md),
                             design_columns = ncol(X), design_rank = qr(X)$rank,
                             formula = paste(deparse(reformulate(names(md))), collapse = " "))
    samples[[key]] <- tibble(composition = community, phenotype = phenotype, Sample_ID = rownames(md))
  }
}
results <- adjust_targets(bind_rows(tables), clinical, q_cutoff) %>%
  dplyr::select(-any_of("host_genus")) %>% left_join(host_map, by = "feature") %>%
  mutate(host_genus = ifelse(composition == "phage_host_genus", feature, host_genus))
write_tsv(results, file.path(res_dir, "all_target_associations.tsv"))
write_tsv(filter(results, significant_q10), file.path(res_dir, "significant_associations_q10.tsv"))
write_tsv(filter(results, supported_q10), file.path(res_dir, "supported_associations_q10.tsv"))
write_tsv(bind_rows(designs), file.path(res_dir, "model_designs.tsv"))
write_tsv(bind_rows(samples), file.path(res_dir, "samples_per_model.tsv"))

# --- 3.4 Additional host adjustment ---------------------------------------
fit_host <- function(cells, community, phenotype, genus, host, mapping) {
  md <- get_metadata(community, phenotype)
  md <- md[rownames(md) %in% rownames(host), , drop = FALSE]
  base_rank <- qr(model.matrix(reformulate(names(md)), md))$rank
  sgbs <- mapping$feature[!is.na(mapping$host_genus) & mapping$host_genus == genus]
  if (clinical) {
    md$HOST_GENUS_CLR <- host[rownames(md), genus]
  } else {
    for (k in seq_along(sgbs)) md[[paste0("HOST_SGB_", k)]] <- host[rownames(md), sgbs[k]]
  }
  X <- model.matrix(reformulate(names(md)), md)
  estimable <- all(is.finite(X)) && nrow(X) > ncol(X) + 1 &&
    if (clinical) qr(X)$rank > base_rank else qr(X)$rank == ncol(X)
  if (!estimable) stop("Host design is not estimable")
  out <- file.path(res_dir, "host_models", community, phenotype, genus)
  x <- abundance[[community]][rownames(md), , drop = FALSE]
  z <- fit_maaslin(x, md, out) %>% filter(metadata == phenotype) %>%
    transmute(feature, value = as.character(value), model,
              host_coef = coef, host_null = null_hypothesis,
              host_stderr = stderr, host_p = pval_individual,
              host_package_q = qval_individual, host_error = as.character(error),
              host_N = N, host_N_not_zero = N_not_zero)
  cells %>% left_join(z, by = c("feature", "value", "model")) %>%
    mutate(host_test_valid = is.na(host_error) & is.finite(host_p),
           fit_status = "completed", host_message = NA_character_)
}

if (include_host) {
  if (!refit) {
    host_results <- read_tsv(file.path(input_dir, "tables", "host_sgb_adjustment_results.tsv"),
                              col_types = cols(value = col_character()), show_col_types = FALSE)
  } else {
    mapping_file <- if (clinical) "bacterial_sgb_to_genus_mapping.tsv" else "host_sgb_gtdb_mapping.tsv"
    host_file <- if (clinical) "bacterial_genus_clr_covariates.rds" else "host_sgb_clr_covariates.rds"
    mapping <- read_tsv(file.path(input_dir, "tables", mapping_file), show_col_types = FALSE)
    host <- readRDS(file.path(input_dir, "cache", host_file))
    candidates <- results %>% filter(composition != "bacteria_species", significant_q10) %>%
      mutate(mapped = !is.na(host_genus) & if (clinical) host_genus %in% colnames(host) else
               host_genus %in% mapping$host_genus)
    batches <- split(candidates, interaction(candidates$composition, candidates$phenotype,
                                             candidates$host_genus, drop = TRUE))
    host_tables <- lapply(batches, function(cells) {
      if (!cells$mapped[1]) return(NULL)
      tryCatch(fit_host(cells, cells$composition[1], cells$phenotype[1], cells$host_genus[1],
                        host, mapping), error = function(e) {
        cells %>% mutate(host_test_valid = FALSE, fit_status = "failed",
                          host_message = conditionMessage(e))
      })
    })
    keys <- c("composition", "phenotype", "feature", "value", "model")
    evaluated <- bind_rows(host_tables)
    if (!nrow(evaluated)) evaluated <- candidates[0, keys]
    host_results <- candidates %>% left_join(
      evaluated %>% dplyr::select(any_of(c(keys, "host_coef", "host_null", "host_stderr",
        "host_p", "host_package_q", "host_error", "host_N", "host_N_not_zero", "host_test_valid",
        "fit_status", "host_message"))), by = keys)
    # Unmapped and failed cells remain in the table; they are not negative tests.
    for (name in c("host_p", "host_package_q", "host_coef", "host_null")) {
      if (!name %in% names(host_results)) host_results[[name]] <- NA_real_
    }
    if (!"host_test_valid" %in% names(host_results)) host_results$host_test_valid <- FALSE
    host_results <- host_results %>%
      mutate(host_test_valid = replace_na(host_test_valid, FALSE), q_host_primary = NA_real_)
    for (community in unique(host_results$composition)) {
      idx <- which(host_results$composition == community & host_results$host_test_valid)
      family_size <- sum(results$composition == community)
      host_results$q_host_primary[idx] <- if (clinical) host_results$host_package_q[idx] else
        p.adjust(host_results$host_p[idx], "BH", n = family_size)
    }
    host_results <- host_results %>%
      mutate(host_box = host_test_valid & !is.na(q_host_primary) & q_host_primary < q_cutoff,
             direction_concordant = ifelse(host_test_valid,
               sign(host_coef - host_null) == sign(coef - null_hypothesis), NA))
  }
  write_tsv(host_results, file.path(res_dir, "host_adjustment_results.tsv"))
}

# --- 3.5 Association heatmaps ---------------------------------------------
# Signed normal-equivalent Z uses the tested null, including abundance medians.
# This compact overview does not replace the original annotated ComplexHeatmaps.
plot_df <- results %>% filter(significant_q10) %>%
  mutate(contrast = paste(phenotype, value, sep = " / "),
         z = sign(effect_vs_null) * qnorm(pmax(pval_individual, .Machine$double.xmin) / 2,
                                          lower.tail = FALSE),
         z = pmax(-5, pmin(5, z)))
if (include_host) {
  plot_df <- plot_df %>% left_join(host_results %>%
    dplyr::select(composition, phenotype, feature, value, model, host_box),
    by = c("composition", "phenotype", "feature", "value", "model"))
} else {
  plot_df$host_box <- FALSE
}
if (nrow(plot_df)) {
  top_features <- plot_df %>% group_by(composition, feature) %>%
    summarise(q = min(q_primary), .groups = "drop") %>%
    group_by(composition) %>% slice_min(q, n = 35, with_ties = FALSE) %>% ungroup()
  plot_cells <- semi_join(plot_df, top_features, by = c("composition", "feature"))
  p_maaslin <- ggplot(plot_cells, aes(x = feature, y = contrast, fill = z)) +
    geom_tile(colour = "white") +
    geom_tile(data = filter(plot_cells, host_box %in% TRUE), fill = NA,
               colour = "black", linewidth = 0.6) +
    facet_wrap(vars(composition, model), scales = "free", ncol = 2) +
    scale_fill_gradient2(low = "#6B9EC1", mid = "#FAFAFA", high = "#EEAB82", limits = c(-5, 5)) +
    labs(x = NULL, y = NULL, fill = "Signed Z",
         caption = "Black outline: significant after host adjustment under the selected version's rule.") + theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(res_dir, "associations_overview.pdf"), p_maaslin, width = 16, height = 12)
}
capture.output(sessionInfo(), file = file.path(res_dir, "sessionInfo.txt"))
