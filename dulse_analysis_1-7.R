# ============================================================
# PACIFIC DULSE ANALYSIS — CORE PIPELINE (SECTIONS 1–7)
# NOTE: In final version, minerals will precede stoichiometry
# ============================================================

library(dplyr)
library(ggplot2)
library(pheatmap)
library(tibble)

# SECTION 1 — USER INPUT + DATA LOAD
mode <- "monthly"

#mode<- "seasonal"

if(mode == "monthly"){
  data <- read.csv("MONTHLY.csv")
} else if(mode == "seasonal"){
  data <- read.csv("SEASONAL.csv")
} else {
  stop("Invalid mode")
}

data <- data %>% rename_all(tolower)

month_levels <- c("DEC","JAN","FEB","MAR","APR","MAY",
                  "JUN","JUL","AUG","SEP","OCT","NOV")

data$month <- factor(data$month, levels = month_levels)

# SECTION 2 — FAA + C:N
faa_cn_clean <- data %>%
  select(month, faa_mg_kg, cn_ratio) %>%
  rename(faa = faa_mg_kg, cn = cn_ratio)

p_faa <- ggplot(faa_cn_clean, aes(x = cn, y = faa)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "pink") +
  labs(title = "FAA vs C:N",
       x = "C:N ratio",
       y = "Total free amino acids (mg/kg)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures/faa_vs_cn.png", p_faa)

# SECTION 3 — ABIOTIC
plot_abiotic <- function(df, var, ylab, title, filename){
  mean_val <- mean(df[[var]], na.rm = TRUE)
  df <- df %>%
    mutate(group = ifelse(.data[[var]] > mean_val,
                          "Above mean", "Below mean"))

  p <- ggplot(df, aes(x = month, y = .data[[var]], group = 1)) +
    geom_line() +
    geom_point(aes(color = group)) +
    geom_hline(yintercept = mean_val, linetype = "dashed") +
    scale_color_manual(values = c("Above mean" = "red",
                                 "Below mean" = "blue")) +
    labs(title = title, x = "Month", y = ylab) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(paste0("figures/", filename), p)
}

plot_abiotic(data, "temp_c", "Temperature (°C)",
             "Seasonal Temperature", "abiotic_temp.png")

plot_abiotic(data, "par_mol_photons_m2",
             "PAR", "Seasonal PAR", "abiotic_par.png")

# SECTION 4 — AA HEATMAP
aa_df <- data %>%
  select(month, ends_with("mgkg")) %>%
  column_to_rownames("month")

aa_scaled <- scale(aa_df)

pheatmap(aa_scaled)

# SECTION 5 — PCA
pca <- prcomp(aa_scaled)

# SECTION 6 — STOICHIOMETRY
minerals <- read.csv("minerals_raw.csv") %>%
  rename_all(tolower)

phosphorus <- minerals %>%
  select(month, p_mgkg)

stoich <- data %>%
  left_join(phosphorus, by = "month") %>%
  rename(c_pct = `c_%`,
         n_pct = `n_%`) %>%
  mutate(
    C_mol = c_pct / 12.011,
    N_mol = n_pct / 14.007,
    P_mol = (p_mgkg / 10000) / 30.974,
    C_to_P = C_mol / P_mol,
    N_to_P = N_mol / P_mol
  )

# SECTION 7 — CORRELATIONS
corr_df <- stoich %>%
  select(month, C_to_P, N_to_P, cn_ratio) %>%
  left_join(faa_cn_clean %>% select(month, faa), by = "month")

cor_cn <- cor.test(corr_df$faa, corr_df$cn_ratio)

cat("CORE PIPELINE COMPLETE")
