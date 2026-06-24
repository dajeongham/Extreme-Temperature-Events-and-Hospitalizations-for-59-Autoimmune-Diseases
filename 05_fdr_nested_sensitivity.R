# 05. Nested-outcome FDR sensitivity (Reviewer #3 Major #2)
# Crohn's disease (K50) and ulcerative colitis (K51) are nested within IBD
# (K50+K51). Re-applies Benjamini-Hochberg within each exposure definition after
# removing the overlapping unit, under two scenarios:
#   (A) keep IBD,  drop subtypes (CD, UC)
#   (B) keep CD/UC, drop aggregate (IBD)
# Data not included (NHIS restriction); set local paths under data/ and output/.

library(dplyr)

# 1] Load main-analysis estimates ---------------------------------------------
#    One row per disease x exposure definition.
#    Columns used: Outcome, Effect (exposure code), RR, LCL, UCL.
df <- read.csv("data/regression_estimates.csv")

# 2] Classify exposure definition and compute two-sided p-value ----------------
df <- df %>%
  mutate(
    tok = sub(" .*", "", Effect),
    RR  = as.numeric(RR), LCL = as.numeric(LCL), UCL = as.numeric(UCL),
    season = if_else(grepl("cold|cw", tok), "Coldspell", "Heatwave"),
    percentile = case_when(
      grepl("heat975|hw97", tok) ~ "97.5%",
      grepl("heat99|hw99",  tok) ~ "99%",
      grepl("heat95|hw95",  tok) ~ "95%",
      grepl("cold2p5|cw2",  tok) ~ "2.5%",
      grepl("cold5|cw5",    tok) ~ "5%",
      grepl("cold1|cw1",    tok) ~ "1%",
      TRUE ~ NA_character_
    ),
    days = case_when(
      grepl("_3", tok) ~ "3 consecutive days",
      grepl("_2", tok) ~ "2 consecutive days",
      TRUE             ~ "1 day"
    ),
    logOR  = log(RR),
    SE     = (log(UCL) - log(LCL)) / (2 * 1.96),
    pvalue = 2 * (1 - pnorm(abs(logOR / SE)))
  ) %>%
  filter(!is.na(percentile), is.finite(pvalue))

# 3] BH within each exposure definition, optionally excluding nested units -----
apply_BH <- function(data, drop = character(0), qname = "qvalue_BH") {
  data %>%
    filter(!(Outcome %in% drop)) %>%
    group_by(season, days, percentile) %>%
    mutate(!!qname := p.adjust(pvalue, method = "BH"), n_pool = n()) %>%
    ungroup()
}

res_full <- apply_BH(df)                                                # full 59-outcome pool
resA     <- apply_BH(df, drop = c("CD", "UC"), qname = "qvalue_BH_red") # keep IBD
resB     <- apply_BH(df, drop = "IBD",          qname = "qvalue_BH_red") # keep CD/UC

# 4] Save full and reduced-pool results ----------------------------------------
write.csv(res_full, "output/FDR_ALL_fullpool_BH.csv",         row.names = FALSE)
write.csv(resA,     "output/FDR_ALL_dropCDUC_keepIBD_BH.csv", row.names = FALSE)
write.csv(resB,     "output/FDR_ALL_dropIBD_keepCDUC_BH.csv", row.names = FALSE)
