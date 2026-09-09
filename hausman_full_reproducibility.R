# ============================================================================
# Hausman Test Monte Carlo -- FULL REPRODUCIBILITY SCRIPT
#
# Reproduces every number reported in the manuscript: the full 192,000-run
# simulation (192 scenarios x 1,000 replications), and the full statistical
# analysis built on it (Tables 1-5 in the published manuscript). This is the
# single script needed to regenerate the paper's results from scratch. Note:
# internal object/file names below (table2, table4, table7, etc.) reflect
# this script's own numbering and do not correspond 1:1 to the manuscript's
# final table numbers after editorial consolidation; see the printed labels
# for what each output actually contains.
#
# Feasibility: the aggregate rank condition for the pooled within-estimator,
# N x (T-1) >= k, holds by at least two orders of magnitude margin for every
# one of the 192 architecture x complexity combinations, so no scenario is
# excluded a priori. (An earlier draft incorrectly excluded scenarios using
# a per-unit T > k rule; that rule was wrong and has been removed. See
# Section 4.1 of the manuscript for the corrected derivation.)
#
# Methodology matches the manuscript exactly (Section 4.5):
#   - ANOVA uses Type II sums of squares (car::Anova), not base R's aov()
#     default (Type I), since Type I is order-dependent and inappropriate
#     for this design.
#   - Pairwise comparisons use two-sample Welch's t-tests (unequal variance),
#     not Bonferroni-corrected tests -- no multiple-comparison correction is
#     applied or claimed.
#   - The architecture x complexity interaction IS included and reported
#     directly, since every combination is now simulated and no cell is
#     structurally excluded.
#   - A logistic regression on the replication-level retain/reject indicator
#     is fit as a robustness check on the ANOVA-on-proportions approach.
#
# Output: written to OUTPUT_DIR. Two categories of files:
#   - Raw simulation results (r_raw_results.csv) -- supplementary
#   - Individual table CSVs -- source data for the manuscript's Tables 1-5
# ============================================================================

suppressMessages({
  library(plm); library(dplyr); library(tidyr); library(purrr); library(car); library(ggplot2)
})

OUTPUT_DIR <- "./results"  # relative to the script location; change if desired
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(1978)
N_REPLICATIONS <- 1000
ALPHA_LEVEL <- 0.05

panel_architectures <- list(
  "Wide Panel" = list(N = 400, T = 4),
  "Square"     = list(N = 200, T = 8),
  "Long Panel" = list(N = 100, T = 16)
)
variable_complexities <- list(
  "Simple" = 1, "Standard" = 3, "Complex" = 6, "High-Dimensional" = 10
)
missingness_mechanisms <- c("Random", "Early Exit", "Late Missing", "Cyclical")
dropout_rates <- c(0.10, 0.20, 0.30, 0.40)

sig_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# ============================================================================
# PART 1: DATA GENERATION AND SIMULATION
# ============================================================================

generate_panel_data <- function(N, T, k) {
  panel_ids <- expand.grid(id = 1:N, time = 1:T)
  alpha <- rnorm(N, mean = 0, sd = 2)
  panel_ids$alpha <- alpha[panel_ids$id]
  x_vars <- matrix(rnorm(N * T * k), nrow = N * T, ncol = k)
  colnames(x_vars) <- paste0("x", 1:k)
  epsilon <- rnorm(N * T, mean = 0, sd = 1.5)
  df <- data.frame(panel_ids, x_vars)
  df$y <- rowSums(df[, paste0("x", 1:k), drop = FALSE]) + df$alpha + epsilon
  df
}

# Missing-data mechanism functions operate on PLAIN data.frame throughout.
# This is the fix: the original implementation carried a plm::pdata.frame
# through these functions, and direct comparison/arithmetic on its pseries
# index columns produced silent errors. Wrapping as pdata.frame happens once,
# immediately before plm::plm() is called.
introduce_random_missing <- function(df, delta) {
  n_obs <- nrow(df)
  missing_indices <- sample(1:n_obs, size = round(delta * n_obs))
  x_cols <- grep("^x", names(df), value = TRUE)
  df[missing_indices, c("y", x_cols)] <- NA
  df
}
introduce_early_exit <- function(df, delta) {
  N <- length(unique(df$id)); T_total <- length(unique(df$time))
  x_cols <- grep("^x", names(df), value = TRUE)
  units_to_drop_count <- max(1, round(2 * N * delta))
  if (delta == 0 || units_to_drop_count == 0) return(df)
  units_to_drop_ids <- sample(unique(df$id), size = units_to_drop_count)
  for (unit_id in units_to_drop_ids) {
    T_exit <- sample(2:T_total, 1)
    df[df$id == unit_id & df$time >= T_exit, c("y", x_cols)] <- NA
  }
  df
}
introduce_late_missing <- function(df, delta) {
  T_total <- length(unique(df$time)); x_cols <- grep("^x", names(df), value = TRUE)
  c_val <- 2 * delta / (T_total + 1)
  df <- df %>% group_by(id) %>%
    mutate(prob_missing = pmin(c_val * time, 0.95), is_missing = runif(n()) < prob_missing) %>%
    ungroup()
  df[df$is_missing, c("y", x_cols)] <- NA
  df$prob_missing <- NULL; df$is_missing <- NULL
  df
}
introduce_cyclical_missing <- function(df, delta) {
  T_total <- length(unique(df$time)); x_cols <- grep("^x", names(df), value = TRUE)
  high_periods <- floor(T_total/2):ceiling(T_total/2 + 1)
  T_high <- length(high_periods); T_low <- T_total - T_high
  if ((3*T_high + T_low) == 0) return(df)
  p_low <- (T_total * delta) / (3*T_high + T_low)
  p_high <- min(3 * p_low, 0.95)
  df$prob <- ifelse(df$time %in% high_periods, p_high, p_low)
  df$is_missing <- runif(nrow(df)) < df$prob
  df[df$is_missing, c("y", x_cols)] <- NA
  df$prob <- NULL; df$is_missing <- NULL
  df
}

run_single_simulation <- function(params) {
  tryCatch({
    complete_data <- generate_panel_data(params$N, params$T, params$k)
    missing_data <- if (params$delta > 0) {
      switch(params$mechanism,
             "Random" = introduce_random_missing(complete_data, params$delta),
             "Early Exit" = introduce_early_exit(complete_data, params$delta),
             "Late Missing" = introduce_late_missing(complete_data, params$delta),
             "Cyclical" = introduce_cyclical_missing(complete_data, params$delta))
    } else complete_data
    obs_per_id <- missing_data %>% filter(!is.na(y)) %>% group_by(id) %>%
      summarise(n_obs = n(), .groups = 'drop')
    valid_ids <- obs_per_id %>% filter(n_obs >= 2) %>% pull(id)  # CORRECTED: a unit needs >=2 observations to contribute any within-unit variation at all; the old n_obs > k threshold was the same per-unit T>k myth applied a second time, post-missingness
    if (length(valid_ids) < 2) {
      return(list(p_value = NA, failure_reason = "Data Failure: Insufficient individuals"))
    }
    final_data <- missing_data %>% filter(id %in% valid_ids) %>%
      plm::pdata.frame(index = c("id", "time"))
    formula <- as.formula(paste("y ~", paste(paste0("x", 1:params$k), collapse = " + ")))
    fe_model <- tryCatch(plm::plm(formula, data = final_data, model = "within"), error = function(e) NULL)
    if (is.null(fe_model)) return(list(p_value = NA, failure_reason = "Model Failure: FE model failed"))
    re_model <- tryCatch(plm::plm(formula, data = final_data, model = "random"), error = function(e) NULL)
    if (is.null(re_model)) return(list(p_value = NA, failure_reason = "Model Failure: RE model failed"))
    ht <- tryCatch(plm::phtest(fe_model, re_model), error = function(e) NULL)
    if (is.null(ht)) return(list(p_value = NA, failure_reason = "Hausman Test Failure"))
    list(p_value = ht$p.value, failure_reason = "Success")
  }, error = function(e) list(p_value = NA, failure_reason = "System Error"))
}

all_scenarios <- expand.grid(
  panel_name = names(panel_architectures), complexity_name = names(variable_complexities),
  mechanism = missingness_mechanisms, delta = dropout_rates, stringsAsFactors = FALSE
) %>% mutate(
  N = map_int(panel_name, ~panel_architectures[[.]]$N),
  T = map_int(panel_name, ~panel_architectures[[.]]$T),
  k = map_int(complexity_name, ~variable_complexities[[.]])
)
feasible_scenarios <- all_scenarios %>% filter(N * (T - 1) >= k)  # CORRECTED: aggregate DOF condition for the pooled within-estimator, not the per-unit T>k heuristic

all_results <- vector("list", nrow(feasible_scenarios))
for (i in 1:nrow(feasible_scenarios)) {
  params <- as.list(feasible_scenarios[i, ])
  scenario_results <- map(1:N_REPLICATIONS, ~ run_single_simulation(params))
  all_results[[i]] <- bind_rows(scenario_results) %>%
    mutate(panel_name = params$panel_name, complexity_name = params$complexity_name,
           mechanism = params$mechanism, delta = params$delta,
           N = params$N, T = params$T, k = params$k,
           specificity = as.numeric(!is.na(p_value) & p_value > ALPHA_LEVEL))
  if (i %% 24 == 0) {
    saveRDS(bind_rows(all_results[1:i]), file.path(OUTPUT_DIR, "r_checkpoint.rds"))
    cat(sprintf("[checkpoint] %d/%d scenarios done\n", i, nrow(feasible_scenarios)))
  }
}
raw <- bind_rows(all_results)
write.csv(raw, file.path(OUTPUT_DIR, "r_raw_results.csv"), row.names = FALSE)
if (file.exists(file.path(OUTPUT_DIR, "r_checkpoint.rds"))) invisible(file.remove(file.path(OUTPUT_DIR, "r_checkpoint.rds")))
cat("Simulation complete:", nrow(raw), "runs.\n\n")

success <- raw %>% filter(failure_reason == "Success")
success$panel_name <- factor(success$panel_name, levels = c("Wide Panel", "Square", "Long Panel"))
success$complexity_name <- factor(success$complexity_name, levels = c("Simple", "Standard", "Complex", "High-Dimensional"))
success$mechanism <- factor(success$mechanism, levels = c("Random", "Early Exit", "Late Missing", "Cyclical"))

# ============================================================================
# PART 2: TABLE 1 -- specificity by mechanism x dropout rate, with per-rate ANOVA
# ============================================================================
table1_rows <- list()
for (d in dropout_rates) {
  row <- list(delta = paste0(d*100, "%"))
  groups <- list()
  for (m in levels(success$mechanism)) {
    sub <- success %>% filter(mechanism == m, delta == d)
    row[[m]] <- sprintf("%.1f%% (n=%d)", mean(sub$specificity)*100, nrow(sub))
    groups[[m]] <- sub$specificity
  }
  fit <- aov(specificity ~ mechanism, data = success %>% filter(delta == d))
  fstat <- summary(fit)[[1]]$`F value`[1]
  pval <- summary(fit)[[1]]$`Pr(>F)`[1]
  row$F_stat <- round(fstat, 2); row$p_value <- round(pval, 3); row$sig <- sig_stars(pval)
  table1_rows[[length(table1_rows)+1]] <- as.data.frame(row, stringsAsFactors = FALSE)
}
table1 <- bind_rows(table1_rows)
write.csv(table1, file.path(OUTPUT_DIR, "table1_mechanism_by_dropout.csv"), row.names = FALSE)

# ============================================================================
# PART 3: TABLE 2 -- overall failure rates by mechanism
# ============================================================================
table2 <- raw %>% group_by(mechanism) %>%
  summarise(total_runs = n(), successful_runs = sum(failure_reason == "Success"),
            failed_runs = total_runs - successful_runs,
            failure_rate = round(failed_runs / total_runs * 100, 2),
            success_rate = round(successful_runs / total_runs * 100, 2), .groups = "drop")
write.csv(table2, file.path(OUTPUT_DIR, "table2_failure_rates.csv"), row.names = FALSE)

# ============================================================================
# PART 4: TABLE 3 -- Long Panel + High-Dimensional outcome breakdown
# ============================================================================
table3 <- raw %>% filter(panel_name == "Long Panel", complexity_name == "High-Dimensional") %>%
  group_by(mechanism, delta) %>%
  summarise(successes = sum(failure_reason == "Success"), n = n(), .groups = "drop")
write.csv(table3, file.path(OUTPUT_DIR, "table3_long_highdim.csv"), row.names = FALSE)

# ============================================================================
# PART 5: TABLE 4 -- specificity by architecture x complexity
# ============================================================================
table4 <- raw %>% group_by(panel_name, complexity_name) %>%
  summarise(N = n(), successful = sum(failure_reason == "Success"),
            failed = N - successful,
            specificity = round(mean(specificity[failure_reason == "Success"]) * 100, 1), .groups = "drop")
write.csv(table4, file.path(OUTPUT_DIR, "table4_architecture_complexity.csv"), row.names = FALSE)

# ============================================================================
# PART 6: TABLE 5 -- 2-way ANOVA (architecture + complexity), Type II SS
# ============================================================================
model5 <- aov(specificity ~ panel_name * complexity_name, data = success)  # UPDATED: includes architecture x complexity interaction, now formally tested per reviewer request
anova5 <- car::Anova(model5, type = 2)
table5 <- data.frame(Effect = rownames(anova5), df = anova5$Df, F_statistic = round(anova5$`F value`, 2),
                      p_value = anova5$`Pr(>F)`, sig = sig_stars(anova5$`Pr(>F)`))
write.csv(table5, file.path(OUTPUT_DIR, "table5_anova_architecture_complexity.csv"), row.names = FALSE)

# ============================================================================
# PART 7: TABLE 6 -- main effects of architecture, with 95% CI
# ============================================================================
table6 <- do.call(rbind, lapply(levels(success$panel_name), function(p) {
  d <- success$specificity[success$panel_name == p]
  m <- mean(d); s <- sd(d); n <- length(d); se <- s / sqrt(n)
  data.frame(panel_name = p, mean_specificity = round(m, 3), sd = round(s, 3), n = n,
             se = round(se, 4), ci_lower = round(m - 1.96*se, 3), ci_upper = round(m + 1.96*se, 3))
}))
write.csv(table6, file.path(OUTPUT_DIR, "table6_main_effects_architecture.csv"), row.names = FALSE)

# ============================================================================
# PART 8: TABLE 7 -- 3-way ANOVA (mechanism + architecture + complexity), Type II SS
# ============================================================================
model7 <- aov(specificity ~ mechanism + factor(delta) + panel_name + complexity_name +
                mechanism:factor(delta) + mechanism:panel_name + mechanism:complexity_name,
              data = success)  # UPDATED: includes mechanism's interactions with rate, architecture, and complexity, per reviewer request
anova7 <- car::Anova(model7, type = 2)
table7 <- data.frame(Effect = rownames(anova7), df = anova7$Df, F_statistic = round(anova7$`F value`, 2),
                      p_value = anova7$`Pr(>F)`, sig = sig_stars(anova7$`Pr(>F)`))
write.csv(table7, file.path(OUTPUT_DIR, "table7_anova_full.csv"), row.names = FALSE)

# ============================================================================
# PART 8b: Logistic regression robustness check on the binary retain/reject outcome
# ============================================================================
model_logit <- glm(specificity ~ mechanism + factor(delta) + panel_name + complexity_name,
                    data = success, family = binomial(link = "logit"))
table7b <- as.data.frame(summary(model_logit)$coefficients)
table7b$term <- rownames(table7b)
write.csv(table7b, file.path(OUTPUT_DIR, "table7b_logistic_robustness.csv"), row.names = FALSE)

# ============================================================================
# PART 9: TABLE 8 -- pairwise comparisons, Welch's t-tests (no correction)
# ============================================================================
pairwise_welch <- function(data, group_col, groups) {
  pairs <- combn(groups, 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pr) {
    da <- data$specificity[data[[group_col]] == pr[1]]
    db <- data$specificity[data[[group_col]] == pr[2]]
    tt <- t.test(da, db)  # Welch's by default (var.equal = FALSE)
    data.frame(comparison = paste(pr[1], "vs", pr[2]),
               p_value = tt$p.value, sig = sig_stars(tt$p.value))
  }))
}
table8_mechanism <- pairwise_welch(success, "mechanism", levels(success$mechanism))
table8_mechanism$factor <- "Missing Data Mechanism"
table8_architecture <- pairwise_welch(success, "panel_name", levels(success$panel_name))
table8_architecture$factor <- "Panel Architecture"
table8_complexity <- pairwise_welch(success, "complexity_name", levels(success$complexity_name))
table8_complexity$factor <- "Model Complexity"
table8 <- bind_rows(table8_mechanism, table8_architecture, table8_complexity) %>%
  select(factor, comparison, p_value, sig)
write.csv(table8, file.path(OUTPUT_DIR, "table8_pairwise_comparisons.csv"), row.names = FALSE)


# ============================================================================
# FIGURES
# ============================================================================
fig_data <- table4
fig_data$panel_name_f <- factor(fig_data$panel_name,
                             levels = c("Long Panel", "Square", "Wide Panel"),
                             labels = c("Long Panel\n(N=100, T=16)", "Square\n(N=200, T=8)", "Wide Panel\n(N=400, T=4)"))
fig_data$complexity_name_f <- factor(fig_data$complexity_name,
                                  levels = c("Simple", "Standard", "Complex", "High-Dimensional"))

# Figure 1: Specificity heatmap (Panel Architecture x Model Complexity)
fig1 <- ggplot(fig_data, aes(x = complexity_name_f, y = panel_name_f, fill = specificity)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", specificity),
                color = specificity < 78),
            fontface = "bold", size = 5.2, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradientn(colors = c("#c0392b", "#f5f2a8", "#1a6e2e"),
                        limits = c(70, 100), name = "Mean\nSpecificity (%)") +
  labs(x = "Model Complexity", y = "Panel Architecture") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 13))
ggsave(file.path(OUTPUT_DIR, "figure1_specificity_heatmap.png"), fig1,
       width = 8.2, height = 4.6, dpi = 300, bg = "white")

# Figure 2: Specificity bar chart (Model Complexity x Panel Architecture)
fig_data$panel_name_short <- factor(fig_data$panel_name,
  levels = c("Wide Panel", "Square", "Long Panel"))
fig_data$complexity_short <- factor(fig_data$complexity_name,
  levels = c("Simple", "Standard", "Complex", "High-Dimensional"),
  labels = c("Simple", "Standard", "Complex", "High-Dim."))

fig2 <- ggplot(fig_data, aes(x = complexity_short, y = specificity, fill = panel_name_short)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%", specificity)),
            position = position_dodge(width = 0.8), vjust = -0.4, angle = 90,
            fontface = "bold", size = 3.4, hjust = 0) +
  scale_fill_manual(values = c("Wide Panel" = "#4472C4", "Square" = "#548235", "Long Panel" = "#C00000")) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20)) +
  labs(x = NULL, y = "Mean Specificity (%)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"))
ggsave(file.path(OUTPUT_DIR, "figure2_specificity_bars.png"), fig2,
       width = 8.2, height = 4.8, dpi = 300, bg = "white")

cat("Both figures saved to", OUTPUT_DIR, "\n\n")

# ============================================================================
# DONE
# ============================================================================
cat("All tables written to", OUTPUT_DIR, "\n\n")
cat("--- Table 2 (mechanism failure rates) ---\n"); print(table2)
cat("\n--- Table 5 (2-way ANOVA) ---\n"); print(table5)
cat("\n--- Table 7 (3-way ANOVA) ---\n"); print(table7)
