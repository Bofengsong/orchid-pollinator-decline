suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(metafor)
  library(purrr)
})

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_file)) {
  dirname(normalizePath(sub("^--file=", "", script_file[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
package_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

meta_path <- file.path(package_root, "data/current/MetaData.csv")
outdir <- file.path(package_root, "outputs/meta_analysis")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

to_pct <- function(x) 100 * (exp(x) - 1)

calc_I2_mv <- function(fit, vi_vec) {
  vi_bar <- mean(vi_vec, na.rm = TRUE)
  100 * sum(fit$sigma2) / (sum(fit$sigma2) + vi_bar)
}

safe_ci <- function(fit) {
  b <- as.numeric(fit$b)
  se <- as.numeric(fit$se)
  list(lb = b - 1.96 * se, ub = b + 1.96 * se)
}

clean_continent <- function(x) {
  x <- str_to_lower(str_squish(as.character(x)))
  case_when(
    str_detect(x, "^asia$") ~ "Asia",
    str_detect(x, "^europe$") ~ "Europe",
    str_detect(x, "^north\\s*america$") ~ "North America",
    str_detect(x, "^south\\s*america$") ~ "South America",
    str_detect(x, "oceania|austral") ~ "Oceania",
    TRUE ~ str_to_title(x)
  )
}

# Broad country-level grouping used for the regional moderator. Central
# American, Caribbean and South American records are combined operationally;
# Mexico remains in North America because the country spans biogeographic
# realms and no Mexican record entered the meta-analysis.
assign_biogeographic_region <- function(country, continent) {
  country <- str_squish(as.character(country))
  out <- as.character(continent)
  out[country %in% c("Costa Rica", "Dominican Republic", "Panama",
                     "Brazil", "Chile", "Ecuador")] <- "Central-South America"
  out[country %in% c("Canada", "Canade", "USA", "USA, Canada",
                     "Canada/USA", "Mexico")] <- "North America"
  out
}

clean_study_method <- function(x) {
  x <- str_to_lower(str_squish(as.character(x)))
  case_when(
    str_detect(x, "field") & str_detect(x, "historical") ~ "Historical + Field",
    str_detect(x, "historical") ~ "Historical",
    str_detect(x, "field") ~ "Field",
    TRUE ~ NA_character_
  )
}

meta0 <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(
    reference_no = as.character(reference_no),
    country = as.character(country),
    pollinator_family = as.character(pollinator_family),
    pollinator_order = str_to_title(str_squish(as.character(pollinator_order))),
    continent = assign_biogeographic_region(country, clean_continent(continent)),
    study_method_clean = clean_study_method(study_method),
    yi = suppressWarnings(as.numeric(yi)),
    vi = suppressWarnings(as.numeric(vi))
  ) %>%
  filter(
    is.finite(yi), is.finite(vi), vi > 0,
    !is.na(reference_no), !is.na(country), !is.na(pollinator_family)
  )

rand_list <- list(~1 | reference_no, ~1 | country, ~1 | pollinator_family)

meta <- meta0 %>%
  mutate(
    reference_no = factor(reference_no),
    country = factor(country),
    pollinator_family = factor(pollinator_family),
    pollinator_order = factor(pollinator_order),
    continent = factor(continent),
    study_method_clean = factor(study_method_clean)
  )

m_main <- rma.mv(yi, vi, random = rand_list, data = meta, method = "REML")
ci_main <- safe_ci(m_main)
I2_main <- calc_I2_mv(m_main, meta$vi)

part1_out <- tibble(
  k = nrow(meta),
  np = n_distinct(meta$reference_no),
  pooled_lnRR = as.numeric(m_main$b),
  ci_lb = ci_main$lb,
  ci_ub = ci_main$ub,
  z = as.numeric(m_main$zval),
  p = as.numeric(m_main$pval),
  pct = to_pct(as.numeric(m_main$b)),
  pct_lb = to_pct(ci_main$lb),
  pct_ub = to_pct(ci_main$ub),
  Q = as.numeric(m_main$QE),
  Q_p = as.numeric(m_main$QEp),
  I2 = I2_main,
  tau2_publication = as.numeric(m_main$sigma2[1]),
  tau2_country = as.numeric(m_main$sigma2[2]),
  tau2_family = as.numeric(m_main$sigma2[3]),
  AIC = AIC(m_main)
)
write_csv(part1_out, file.path(outdir, "Part1_overall_model_summary_corrected.csv"))

qm_categorical <- function(df, varname, rand = rand_list) {
  df1 <- df %>% filter(!is.na(.data[[varname]]))
  if (n_distinct(df1[[varname]]) < 2) return(NULL)
  m_qm <- rma.mv(yi, vi, mods = reformulate(varname), random = rand, data = df1, method = "REML")
  tibble(
    moderator = varname,
    k = nrow(df1),
    np = n_distinct(df1$reference_no),
    Qm = as.numeric(m_qm$QM),
    df = ifelse(length(m_qm$QMdf) == 0, NA, as.numeric(m_qm$QMdf[1])),
    p = as.numeric(m_qm$QMp),
    AIC = AIC(m_qm)
  )
}

subgroup_mv <- function(df, group_var, rand = rand_list, min_k = 2) {
  lvls <- df %>%
    filter(!is.na(.data[[group_var]])) %>%
    pull(.data[[group_var]]) %>%
    unique() %>%
    as.character()

  map_dfr(lvls, function(g) {
    dsub <- df %>% filter(.data[[group_var]] == g)
    if (nrow(dsub) < min_k) {
      return(tibble(
        subgroup_var = group_var, level = g, k = nrow(dsub),
        np = n_distinct(dsub$reference_no), yi = NA_real_, se = NA_real_,
        ci_lb = NA_real_, ci_ub = NA_real_, z = NA_real_, p = NA_real_,
        pct = NA_real_, pct_lb = NA_real_, pct_ub = NA_real_, I2 = NA_real_
      ))
    }
    fit <- rma.mv(yi, vi, random = rand, data = dsub, method = "REML")
    ci <- safe_ci(fit)
    tibble(
      subgroup_var = group_var,
      level = g,
      k = nrow(dsub),
      np = n_distinct(dsub$reference_no),
      yi = as.numeric(fit$b),
      se = as.numeric(fit$se),
      ci_lb = ci$lb,
      ci_ub = ci$ub,
      z = as.numeric(fit$zval),
      p = as.numeric(fit$pval),
      pct = to_pct(as.numeric(fit$b)),
      pct_lb = to_pct(ci$lb),
      pct_ub = to_pct(ci$ub),
      I2 = calc_I2_mv(fit, dsub$vi)
    )
  })
}

qm_tab_meta <- bind_rows(
  qm_categorical(meta, "pollinator_order"),
  qm_categorical(meta, "continent"),
  qm_categorical(meta, "study_method_clean")
)
write_csv(qm_tab_meta, file.path(outdir, "Part2A_Qm_tests_meta_moderators_corrected.csv"))

sub_order <- subgroup_mv(meta, "pollinator_order")
sub_cont <- subgroup_mv(meta, "continent")
sub_method <- subgroup_mv(meta, "study_method_clean")
write_csv(sub_order, file.path(outdir, "Part2B_subgroup_order_corrected.csv"))
write_csv(sub_cont, file.path(outdir, "Part2B_subgroup_continent_corrected.csv"))
write_csv(sub_method, file.path(outdir, "Part2B_subgroup_study_method_corrected.csv"))

target_drivers <- c("Multiple", "Habitat", "Climate", "Urbanization", "Pesticides", "Pathogens", "Invasives")

long0 <- meta0 %>%
  select(
    reference_no, start_year, end_year, publish_year, country, region,
    latitude, longitude, pollinator_order, pollinator_family, pollinator_species,
    study_method, trend, decline_value, sample_size, yi, vi, continent,
    drivers_reported
  ) %>%
  mutate(drivers_reported = as.character(drivers_reported)) %>%
  separate_rows(drivers_reported, sep = ",") %>%
  mutate(driver = str_to_title(str_squish(drivers_reported))) %>%
  filter(driver %in% target_drivers) %>%
  select(-drivers_reported)

write_csv(long0, file.path(outdir, "MainMetaDataLong_corrected.csv"))

long <- long0 %>%
  mutate(
    reference_no = factor(reference_no),
    country = factor(country),
    pollinator_family = factor(pollinator_family),
    driver = factor(driver, levels = target_drivers)
  )

sub_driver <- subgroup_mv(long, "driver")
write_csv(sub_driver, file.path(outdir, "Part2B_subgroup_driver_corrected.csv"))

m_qm_drv <- rma.mv(yi, vi, mods = ~ driver, random = rand_list, data = long, method = "REML")
qm_driver_row <- tibble(
  moderator = "driver",
  k = nrow(long),
  np = n_distinct(long$reference_no),
  Qm = as.numeric(m_qm_drv$QM),
  df = ifelse(length(m_qm_drv$QMdf) == 0, NA, as.numeric(m_qm_drv$QMdf[1])),
  p = as.numeric(m_qm_drv$QMp),
  AIC = AIC(m_qm_drv)
)
write_csv(qm_driver_row, file.path(outdir, "Part2B_Qm_test_driver_corrected.csv"))

png(file.path(outdir, "FigureS12a_Funnel_multilevel_corrected.png"), width = 1600, height = 1200, res = 200)
funnel(m_main, xlab = "Effect size (lnRR)", ylab = "Standard Error",
       main = "Funnel plot (multilevel model, corrected data)", refline = 0)
dev.off()

meta_bias <- meta %>%
  mutate(resid = residuals(m_main), sei = sqrt(vi))
egger_mv <- lm(resid ~ sei, data = meta_bias)
egger_tab <- as.data.frame(summary(egger_mv)$coefficients) %>%
  tibble::rownames_to_column("term") %>%
  as_tibble() %>%
  rename(estimate = Estimate, std_error = `Std. Error`, statistic = `t value`, p_value = `Pr(>|t|)`)
write_csv(egger_tab, file.path(outdir, "TableS4_Egger_like_multilevel_residuals_corrected.csv"))

agg <- meta0 %>%
  group_by(reference_no) %>%
  group_modify(~{
    d <- .x
    if (nrow(d) == 1) {
      tibble(yi = d$yi[1], vi = d$vi[1])
    } else {
      fit_i <- rma(yi, vi, data = d, method = "REML")
      tibble(yi = as.numeric(fit_i$b), vi = as.numeric(fit_i$se)^2)
    }
  }) %>%
  ungroup()

m_study <- rma(yi, vi, data = agg, method = "REML")
egger_study <- regtest(m_study, model = "rma", predictor = "sei")
tf_study <- trimfill(m_study)

tab_tf <- tibble(
  model = c("study_level_before", "study_level_after_trimfill"),
  lnRR = c(as.numeric(m_study$b), as.numeric(tf_study$b)),
  ci_lb = c(as.numeric(m_study$ci.lb), as.numeric(tf_study$ci.lb)),
  ci_ub = c(as.numeric(m_study$ci.ub), as.numeric(tf_study$ci.ub)),
  pct = c(to_pct(as.numeric(m_study$b)), to_pct(as.numeric(tf_study$b))),
  k0 = c(NA_integer_, tf_study$k0)
)
write_csv(tab_tf, file.path(outdir, "TableS5_TrimFill_summary_corrected.csv"))

study_bias_tab <- tibble(
  diagnostic = c("study_level_egger_z", "study_level_egger_p", "trimfill_missing_studies"),
  value = c(as.numeric(egger_study$zval), as.numeric(egger_study$pval), as.numeric(tf_study$k0))
)
write_csv(study_bias_tab, file.path(outdir, "TableS5_study_level_bias_tests_corrected.csv"))

png(file.path(outdir, "FigureS12b_Funnel_trimfill_corrected.png"), width = 2000, height = 900, res = 200)
par(mfrow = c(1, 2))
funnel(m_study, main = "Funnel (study-level, corrected)", xlab = "lnRR", ylab = "SE")
funnel(tf_study, main = "Funnel (trim & fill, corrected)", xlab = "lnRR", ylab = "SE")
par(mfrow = c(1, 1))
dev.off()

inf <- influence(m_study)
png(file.path(outdir, "FigureS13_Influence_diagnostics_corrected.png"), width = 1600, height = 1200, res = 200)
plot(inf)
dev.off()

m_ref <- rma.mv(yi, vi, random = list(~1 | reference_no), data = meta, method = "REML")
m_ref_cty <- rma.mv(yi, vi, random = list(~1 | reference_no, ~1 | country), data = meta, method = "REML")
m_ref_cty_fam <- m_main

sens_re <- tibble(
  model = c("ref_only", "ref+country", "ref+country+family"),
  lnRR = c(as.numeric(m_ref$b), as.numeric(m_ref_cty$b), as.numeric(m_ref_cty_fam$b)),
  ci_lb = c(safe_ci(m_ref)$lb, safe_ci(m_ref_cty)$lb, safe_ci(m_ref_cty_fam)$lb),
  ci_ub = c(safe_ci(m_ref)$ub, safe_ci(m_ref_cty)$ub, safe_ci(m_ref_cty_fam)$ub),
  pct = to_pct(c(as.numeric(m_ref$b), as.numeric(m_ref_cty$b), as.numeric(m_ref_cty_fam$b))),
  AIC = c(AIC(m_ref), AIC(m_ref_cty), AIC(m_ref_cty_fam)),
  tau2_total = c(sum(m_ref$sigma2), sum(m_ref_cty$sigma2), sum(m_ref_cty_fam$sigma2))
)
write_csv(sens_re, file.path(outdir, "TableS2_random_structure_sensitivity_corrected.csv"))

tab_s2_exact <- tibble(
  Model = c("Model 1", "Model 2", "Model 3 (Main)"),
  `Random Effect` = c("Publications", "Publications + Country", "Publications + Country + Family"),
  `Estimate (lnRR)` = c(as.numeric(m_ref$b), as.numeric(m_ref_cty$b), as.numeric(m_main$b)),
  SE = c(as.numeric(m_ref$se), as.numeric(m_ref_cty$se), as.numeric(m_main$se)),
  `z-value` = c(as.numeric(m_ref$zval), as.numeric(m_ref_cty$zval), as.numeric(m_main$zval)),
  AIC = c(AIC(m_ref), AIC(m_ref_cty), AIC(m_main)),
  `p-value` = ifelse(c(as.numeric(m_ref$pval), as.numeric(m_ref_cty$pval), as.numeric(m_main$pval)) < 0.001,
                     "< 0.001",
                     sprintf("%.4f", c(as.numeric(m_ref$pval), as.numeric(m_ref_cty$pval), as.numeric(m_main$pval))))
)
write_csv(tab_s2_exact, file.path(outdir, "TableS2_random_structure_sensitivity_exact_corrected.csv"))

w_orig <- 1 / meta0$vi
cap <- quantile(w_orig, 0.95, na.rm = TRUE)

meta_cap <- meta0 %>%
  mutate(vi_cap = 1 / pmin(1 / vi, cap)) %>%
  filter(is.finite(vi_cap), vi_cap > 0) %>%
  mutate(reference_no = factor(reference_no), country = factor(country), pollinator_family = factor(pollinator_family))

meta_eq <- meta0 %>%
  mutate(vi_eq = 1) %>%
  mutate(reference_no = factor(reference_no), country = factor(country), pollinator_family = factor(pollinator_family))

m_cap <- rma.mv(yi, vi_cap, random = rand_list, data = meta_cap, method = "REML")
m_eq <- rma.mv(yi, vi_eq, random = rand_list, data = meta_eq, method = "REML")

sens_w <- tibble(
  model = c("original", "cap_95pct_weight", "equal_weight"),
  lnRR = c(as.numeric(m_main$b), as.numeric(m_cap$b), as.numeric(m_eq$b)),
  ci_lb = c(safe_ci(m_main)$lb, safe_ci(m_cap)$lb, safe_ci(m_eq)$lb),
  ci_ub = c(safe_ci(m_main)$ub, safe_ci(m_cap)$ub, safe_ci(m_eq)$ub),
  pct = to_pct(c(as.numeric(m_main$b), as.numeric(m_cap$b), as.numeric(m_eq$b))),
  pct_lb = to_pct(c(safe_ci(m_main)$lb, safe_ci(m_cap)$lb, safe_ci(m_eq)$lb)),
  pct_ub = to_pct(c(safe_ci(m_main)$ub, safe_ci(m_cap)$ub, safe_ci(m_eq)$ub)),
  AIC = c(AIC(m_main), AIC(m_cap), AIC(m_eq))
)
write_csv(sens_w, file.path(outdir, "TableS3_weight_sensitivity_corrected.csv"))

combined_subgroups <- bind_rows(
  sub_driver %>% mutate(panel = "Driver"),
  sub_order %>% mutate(panel = "Order"),
  sub_cont %>% mutate(panel = "Continent"),
  sub_method %>% mutate(panel = "Study method")
)
write_csv(combined_subgroups, file.path(outdir, "Figure2_subgroup_summary_all_clean_corrected.csv"))

cat("Corrected meta-analysis outputs saved to: ", outdir, "\n", sep = "")
