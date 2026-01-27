library(readr)
library(dplyr)
library(stringr)

df_avg <- read_csv("Liver_summary.all_loci.with_tpms.blastp_human.csv", show_col_types = FALSE) %>%
  mutate(
    mean_tpm = sapply(tpms, function(x) {
      if (is.na(x) || str_trim(x) == "") return(NA_real_)
      vals <- as.numeric(str_split(x, ",", simplify = TRUE))
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) NA_real_ else mean(vals)
    })
  ) %>%
  filter(
    max_prob > 0.9,
    is.na(human_best_pident) | human_best_pident < 50,
    mean_tpm >= 5
  )

write_csv(df_avg, "Liver_summary.filtered.csv")
