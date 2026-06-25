# 04. Multiplicity correction: Benjamini-Hochberg (BH) + Benjamini-Yekutieli (BY)
# Applied within each exposure definition (season x duration x percentile) across
# the 59 autoimmune diseases (Reviewer #3 Major #1, #2).
#   BH: valid under positive regression dependence (PRDS).
#   BY: valid under arbitrary dependence (more conservative).

library(dplyr)

# 1] Load regression estimates ------------------------------------------------
#    One row per disease x exposure definition x subset.
#    Columns used: Outcome, Type (exposure code), RR, LCL, UCL, subset.
df <- read.csv("data/regression_estimates.csv")

# 2] Classify each exposure definition from the exposure code (Type) -----------
df <- df %>%
  mutate(
    RR  = as.numeric(RR),
    LCL = as.numeric(LCL),
    UCL = as.numeric(UCL),
    season = if_else(grepl("cold|cw", Type), "Coldspell", "Heatwave"),
    percentile = case_when(
      grepl("cold1|cw1",   Type) ~ "1%",
      grepl("heat99|hw99", Type) ~ "99%",
      grepl("cold2p5|cw2", Type) ~ "2.5%",
      grepl("heat975|hw97",Type) ~ "97.5%",
      grepl("cold5|cw5",   Type) ~ "5%",
      grepl("heat95|hw95", Type) ~ "95%",
      TRUE ~ NA_character_
    ),
    days = case_when(
      grepl("_3", Type) ~ "3 consecutive days",
      grepl("_2", Type) ~ "2 consecutive days",
      TRUE              ~ "1 day"
    )
  )

# 3] Two-sided p-value from the odds ratio and its 95% CI ----------------------
df <- df %>%
  mutate(
    logOR  = log(RR),
    SE     = (log(UCL) - log(LCL)) / (2 * 1.96),
    pvalue = 2 * (1 - pnorm(abs(logOR / SE)))
  )

# 4] FDR within each exposure definition (BH and BY) ---------------------------
result <- df %>%
  group_by(subset, season, days, percentile) %>%
  mutate(
    qvalue_BH = p.adjust(pvalue, method = "BH"),
    qvalue_BY = p.adjust(pvalue, method = "BY"),
    sig_BH    = qvalue_BH < 0.05,
    sig_BY    = qvalue_BY < 0.05
  ) %>%
  ungroup()

# 5] Save full table and a BH-vs-BY robustness summary -------------------------
write.csv(result, "output/FDR_59_BH_BY.csv", row.names = FALSE)

robustness <- result %>%
  group_by(subset, season) %>%
  summarise(
    n_sig_BH     = sum(sig_BH, na.rm = TRUE),
    n_sig_BY     = sum(sig_BY, na.rm = TRUE),
    n_BH_also_BY = sum(sig_BH & sig_BY, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(robustness, "output/FDR_BH_BY_robustness_summary.csv", row.names = FALSE)

# References:
#   Benjamini Y, Hochberg Y (1995). J R Stat Soc B 57(1):289-300.
#   Benjamini Y, Yekutieli D (2001). Ann Stat 29(4):1165-1188.
