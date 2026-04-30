# ============================================================
# PACIFIC DULSE ANALYSIS — CORE PIPELINE (SECTIONS 1–7)
# ============================================================

# -------------------------
# GLOBAL OUTPUT SETTINGS (MUST BE FIRST)
# -------------------------
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

save_fig <- function(filename, plot = NULL, width = 8, height = 5, dpi = 300){
  path <- file.path(fig_dir, filename)
  if(is.null(plot)){
    ggsave(path, width = width, height = height, dpi = dpi)
  } else {
    ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  }
}

# -------------------------
# LIBRARIES
# -------------------------
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tibble)

# ============================================================
# SECTION 1 — USER INPUT + DATA LOAD
# ============================================================

mode <- "monthly"

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

# ============================================================
# SECTION 2 — FAA + C:N
# ============================================================

faa_cn_clean <- data %>%
  select(month, faa_mg_kg, cn_ratio) %>%
  rename(faa = faa_mg_kg, cn = cn_ratio)

p_faa <- ggplot(faa_cn_clean, aes(x = cn, y = faa)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "pink") +
  labs(title = "FAA vs C:N",
       x = "C:N ratio",
       y = "Total free amino acids (mg/kg)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

save_fig("faa_vs_cn.png", p_faa)

# ============================================================
# SECTION 3 — ABIOTIC
# ============================================================

plot_abiotic <- function(df, var, ylab, title, filename){
  mean_val <- mean(df[[var]], na.rm = TRUE)

  df <- df %>%
    mutate(group = ifelse(.data[[var]] > mean_val,
                          "Above mean", "Below mean"))

  p <- ggplot(df, aes(x = month, y = .data[[var]], group = 1)) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(aes(color = group), size = 3) +
    geom_hline(yintercept = mean_val, linetype = "dashed") +
    scale_color_manual(values = c("Above mean" = "red",
                                 "Below mean" = "blue")) +
    labs(title = title, x = "Month", y = ylab) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))

  save_fig(filename, p)
}

plot_abiotic(data, "temp_c", "Temperature (°C)",
             "Temperature", "abiotic_temp.png")

plot_abiotic(data, "par_mol_photons_m2",
             "PAR", "PAR", "abiotic_par.png")

# ============================================================
# SECTION 4 — AMINO ACID HEATMAP
# ============================================================

aa_df <- data %>%
  select(month, ends_with("mgkg")) %>%
  column_to_rownames("month")

aa_scaled <- scale(aa_df)

pheatmap(
  aa_scaled,
  color = colorRampPalette(c("blue","white","red"))(50),
  main = "Seasonal Amino Acid Profile (Z-scored)",
  filename = file.path(fig_dir, "aa_heatmap.png")
)

# ============================================================
# SECTION 5 — PCA (AMINO ACIDS)
# ============================================================

pca <- prcomp(aa_scaled)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

loadings <- as.data.frame(pca$rotation)

png(file.path(fig_dir, "aa_pca.png"), width = 1000, height = 800)

plot(scores$PC1, scores$PC2,
     xlab = "PC1",
     ylab = "PC2",
     pch = 19)

text(scores$PC1, scores$PC2, labels = scores$month, pos = 3)

arrows(0, 0,
       loadings[,1]*3,
       loadings[,2]*3,
       col = "red", length = 0.1)

text(loadings[,1]*3,
     loadings[,2]*3,
     labels = rownames(loadings),
     col = "blue")

dev.off()

# ============================================================
# SECTION 6 — STOICHIOMETRY (C:N:P)
# ============================================================

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
    N_to_P = N_mol / P_mol,
    CNP_ratio = paste0(round(C_to_P), ":", round(N_to_P), ":1")
  )

write.csv(stoich, file.path("outputs", "stoichiometry_results.csv"), row.names = FALSE)

# ============================================================
# SECTION 7 — CORRELATIONS
# ============================================================

corr_df <- stoich %>%
  select(month, C_to_P, N_to_P, cn_ratio) %>%
  left_join(faa_cn_clean %>% select(month, faa), by = "month")

cor_cn <- cor.test(corr_df$faa, corr_df$cn_ratio)
cor_np <- cor.test(corr_df$faa, corr_df$N_to_P)
cor_cp <- cor.test(corr_df$faa, corr_df$C_to_P)

write.csv(data.frame(
  comparison = c("FAA vs C:N","FAA vs N:P","FAA vs C:P"),
  r = c(cor_cn$estimate, cor_np$estimate, cor_cp$estimate),
  p = c(cor_cn$p.value, cor_np$p.value, cor_cp$p.value)
), file.path("outputs", "correlation_results.csv"), row.names = FALSE)

cat("\nCORE PIPELINE COMPLETE\n")
