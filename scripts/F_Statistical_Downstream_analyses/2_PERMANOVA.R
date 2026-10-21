################################################################################
##### 2000HIV: PERMANOVA analysis
### Author(s):Haoran Peng
### Last updated: 19th October, 2026
################################################################################

## Uses the prepared metadata, phenotype panel and distances in Result_V1.0.
## Phage: Bray-Curtis on log2(1 + TSS relative abundance).
## Bacteria: Euclidean on pseudocount CLR, as in the saved analysis.
## Marginal tests; 999 permutations within CENTER; BH within each community.
## The saved bacterial distance includes 6 archaeal and 2 eukaryotic features.
## This script organises the existing analysis without changing that distance.

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)
library(permute)

root <- Sys.getenv("HIV_PROJECT_DIR", "/scratch/p312334/project/16--HIV_Immun")
input_dir <- file.path(root, "Result_V1.0", "Permanova")
res_dir <- file.path(root, "results_clean", "02_PERMANOVA")
refit <- FALSE                       # FALSE: collect existing fits; TRUE: rerun.
n_perm <- 999L
seed <- 20260713L
covar <- c("AGE", "BMI_BASELINE", "CENTER", "ETHNICITY_COLLAPSED", "log10_counts")
communities <- c("bacteria_species", "phage")

# --- 2.1 Prepared inputs ---------------------------------------------------
metadata <- as.data.frame(readRDS(file.path(input_dir, "cache", "metadata.rds")))
panel <- read_tsv(file.path(input_dir, "tables", "phenotype_panel.tsv"),
                  show_col_types = FALSE)
display <- read_tsv(file.path(input_dir, "tables", "phenotype_display_names.tsv"),
                    show_col_types = FALSE)
stopifnot(!anyDuplicated(metadata$Sample_ID), !anyDuplicated(panel$phenotype))

# Feature filtering and distance construction are upstream steps. Reusing the
# prepared distances preserves the original cohort and empty-phage convention.
distances <- setNames(lapply(communities, function(community) {
  readRDS(file.path(input_dir, "cache", paste0(community, "_distance.rds")))
}), communities)

# --- 2.2 One phenotype and one community ----------------------------------
prepare_cases <- function(metadata, phenotype, type, covar) {
  df <- metadata[, c("Sample_ID", phenotype, covar), drop = FALSE]
  value <- as.character(df[[phenotype]])
  if (type == "continuous") {
    df[[phenotype]] <- as.numeric(value)
  } else if (type == "binary_numeric") {
    value <- trimws(value)
    value[value == ""] <- NA_character_
    value <- sub("^([01])\\.0+$", "\\1", value)
    df[[phenotype]] <- factor(value, levels = c("0", "1"))
  } else {
    df[[phenotype]] <- factor(value)
  }
  df <- df[complete.cases(df), , drop = FALSE]
  removed <- character()
  if (type == "categorical_group") {
    counts <- table(df[[phenotype]])
    removed <- names(counts)[counts < 10]
    df <- df[!df[[phenotype]] %in% removed, , drop = FALSE]
  }
  df <- droplevels(df)
  rownames(df) <- df$Sample_ID
  stopifnot(nrow(df) >= 30, n_distinct(df[[phenotype]]) >= 2,
            all(vapply(df[covar], n_distinct, integer(1)) > 1))
  list(data = df, removed = removed)
}

fit_permanova <- function(distance, metadata, phenotype, type, covar, n_perm, seed) {
  cases <- prepare_cases(metadata, phenotype, type, covar)
  df <- cases$data
  d <- as.dist(distance[df$Sample_ID, df$Sample_ID, drop = FALSE])
  formula <- reformulate(c(phenotype, covar), response = "d")
  set.seed(seed)
  permutations <- shuffleSet(nrow(df), nset = n_perm,
                              control = how(blocks = df$CENTER))
  fit <- adonis2(formula, data = df, permutations = permutations,
                 by = "margin", parallel = 1)
  anova <- as.data.frame(fit) %>% tibble::rownames_to_column("term")
  hit <- anova[anova$term %in% c(phenotype, paste0("`", phenotype, "`")), ]
  stopifnot(nrow(hit) == 1, is.finite(hit$`Pr(>F)`))
  result <- tibble(
    phenotype = phenotype, phenotype_type = type, n = nrow(df),
    n_levels = n_distinct(df[[phenotype]]), df = hit$Df,
    r2 = hit$R2, f_stat = hit$F, pval = hit$`Pr(>F)`,
    permutations = n_perm, seed = seed,
    permutation_scheme = "restricted within CENTER",
    formula = paste(deparse(formula), collapse = " "),
    removed_rare_levels = paste(cases$removed, collapse = ";"), status = "ok"
  )
  list(result = result, anova = anova, sample_ids = df$Sample_ID)
}

# --- 2.3 Run or collect ----------------------------------------------------
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
fits <- list()
for (community in communities) {
  for (i in seq_len(nrow(panel))) {
    phenotype <- panel$phenotype[i]
    message(community, " / ", phenotype)
    if (refit) {
      fit <- fit_permanova(distances[[community]], metadata, phenotype,
                           panel$type[i], covar, n_perm, seed)
      fit$result$composition <- community
      fit$result$category <- panel$category[i]
      fit$result$distance <- if (community == "phage") "Bray-Curtis" else "Aitchison (CLR-Euclidean)"
      saveRDS(fit, file.path(res_dir, paste0(community, "_", phenotype, ".rds")))
    } else {
      task <- (match(community, communities) - 1L) * nrow(panel) + i - 1L
      fit <- readRDS(file.path(input_dir, "cache", sprintf("fit_%02d.rds", task)))
      stopifnot(fit$result$phenotype == phenotype,
                fit$result$composition == community,
                fit$result$permutations == n_perm, fit$result$seed == seed)
    }
    fits[[paste(community, phenotype, sep = "/")]] <- fit
  }
}
results <- bind_rows(lapply(fits, `[[`, "result")) %>%
  group_by(composition) %>%
  mutate(qval = p.adjust(pval, "BH"), significant_q10 = qval < 0.1) %>%
  ungroup() %>%
  left_join(display %>% dplyr::select(phenotype, display_name), by = "phenotype")
samples <- bind_rows(lapply(fits, function(fit) {
  tibble(composition = fit$result$composition, phenotype = fit$result$phenotype,
         Sample_ID = fit$sample_ids)
}))
write_tsv(results, file.path(res_dir, "permanova_results.tsv"))
write_tsv(filter(results, significant_q10), file.path(res_dir, "permanova_significant_q10.tsv"))
write_tsv(samples, file.path(res_dir, "samples_per_model.tsv"))
saveRDS(fits, file.path(res_dir, "permanova_models.rds"))

# --- 2.4 Plot marginal R2 -------------------------------------------------
# Marginal R2 is adjusted for the other terms, not degrees-of-freedom adjusted R2.
plot_df <- results %>%
  mutate(label = paste(composition, display_name, sep = " / "),
         label = reorder(label, r2), star = ifelse(significant_q10, "*", ""))
p_permanova <- ggplot(plot_df, aes(x = label, y = 100 * r2, fill = composition)) +
  geom_col(width = 0.7, colour = "white") +
  geom_text(aes(label = star), hjust = -0.2, size = 4) +
  coord_flip() +
  facet_wrap(~ composition, scales = "free_y") +
  scale_fill_manual(values = c(bacteria_species = "#EAAC91", phage = "#B0A0C1")) +
  scale_x_discrete(labels = function(x) sub("^[^/]+ / ", "", x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Marginal variance explained (%)", caption = "* BH q < 0.1") +
  theme_bw() + theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(res_dir, "permanova.pdf"), p_permanova, width = 12, height = 11)
capture.output(sessionInfo(), file = file.path(res_dir, "sessionInfo.txt"))
p_permanova
