# ============================================================
# PACIFIC DULSE — FULL REPRODUCIBLE PIPELINE (SECTIONS 1–8)
# ============================================================

# ============================================================
# GLOBAL SETTINGS (ALWAYS FIRST)
# ============================================================

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

# ============================================================
# LIBRARIES
# ============================================================

library(dplyr)
library(ggplot2)
library(pheatmap)
library(tibble)

# ============================================================
# SECTION 1 — USER INPUT + DATA LOAD
# ============================================================

mode <- "monthly"   # change to "seasonal" if needed

if(mode == "monthly"){
  data <- read.csv("MONTHLY.csv")
} else {
  data <- read.csv("SEASONAL.csv")
}

data <- data %>% rename_all(tolower)

month_levels <- c("DEC","JAN","FEB","MAR","APR","MAY",
                  "JUN","JUL","AUG","SEP","OCT","NOV")

data$month <- factor(data$month, levels = month_levels)

# ============================================================
# SECTION 2 — FAA vs C:N
# ============================================================

faa_cn <- data %>%
  select(month, faa_mg_kg, cn_ratio) %>%
  rename(faa = faa_mg_kg, cn = cn_ratio)

p <- ggplot(faa_cn, aes(x = cn, y = faa)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "FAA vs C:N",
       x = "C:N ratio",
       y = "FAA (mg/kg)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

save_fig("faa_vs_cn.png", p)

# ============================================================
# SECTION 3 — ABIOTIC
# ============================================================

plot_abiotic <- function(df, var, ylab, title, filename){
  mean_val <- mean(df[[var]], na.rm = TRUE)

  df <- df %>%
    mutate(group = ifelse(.data[[var]] > mean_val,
                          "Above mean", "Below mean"))

  p <- ggplot(df, aes(x = month, y = .data[[var]], group = 1)) +
    geom_line(color = "black") +
    geom_point(aes(color = group), size = 3) +
    geom_hline(yintercept = mean_val, linetype = "dashed") +
    scale_color_manual(values = c("Above mean"="red","Below mean"="blue")) +
    labs(title = title, x = "Month", y = ylab) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  save_fig(filename, p)
}

plot_abiotic(data, "temp_c", "Temperature (°C)", "Temperature", "abiotic_temp.png")
plot_abiotic(data, "par_mol_photons_m2", "PAR", "PAR", "abiotic_par.png")

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
  main = "Amino Acid Profile",
  filename = file.path(fig_dir, "aa_heatmap.png")
)

# ============================================================
# SECTION 5 — AMINO ACID PCA
# ============================================================

pca <- prcomp(aa_scaled)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

loadings <- as.data.frame(pca$rotation)

png(file.path(fig_dir, "aa_pca.png"), width = 1000, height = 800)

plot(scores$PC1, scores$PC2, pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "AA PCA")

text(scores$PC1, scores$PC2, labels = scores$month, pos = 3)

arrows(0,0, loadings[,1]*3, loadings[,2]*3, col="red")

dev.off()

# ============================================================
# SECTION 6 — MINERALS (LOAD + PREP)
# ============================================================

minerals <- read.csv("minerals_raw.csv") %>%
  rename_all(tolower)

minerals$month <- factor(minerals$month, levels = month_levels)

mat <- minerals %>%
  arrange(month) %>%
  column_to_rownames("month") %>%
  as.matrix()

mat <- apply(mat, 2, as.numeric)

scale_safe <- function(x){
  s <- sd(x, na.rm = TRUE)
  if(s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

mat_scaled <- apply(mat, 2, scale_safe)
rownames(mat_scaled) <- rownames(mat)

# ============================================================
# SECTION 7 — MINERAL HEATMAP
# ============================================================

pheatmap(
  mat_scaled,
  color = colorRampPalette(c("blue","white","red"))(50),
  main = "Mineral Composition",
  filename = file.path(fig_dir, "mineral_heatmap.png")
)

# ============================================================
# SECTION 8 — MINERAL PCA
# ============================================================

pca <- prcomp(mat_scaled, center = FALSE)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

loadings <- as.data.frame(pca$rotation[,1:2])
loadings$element <- rownames(loadings)

threshold <- 0.25

loadings <- loadings %>%
  filter(abs(PC1) > threshold | abs(PC2) > threshold)

symbol_map <- c(
  b="B", ca="Ca", cu="Cu", fe="Fe", k="K", mg="Mg",
  mn="Mn", mo="Mo", na="Na", p="P", s="S", zn="Zn",
  as="As", cr="Cr", cd="Cd", co="Co", ni="Ni", pb="Pb"
)

loadings$element <- symbol_map[loadings$element]

png(file.path(fig_dir, "mineral_pca.png"), width = 1000, height = 800)

plot(scores$PC1, scores$PC2,
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "Mineral PCA")

text(scores$PC1, scores$PC2,
     labels = scores$month,
     pos = 3)

arrows(0,0,
       loadings$PC1*4,
       loadings$PC2*4,
       col="red")

text(loadings$PC1*4,
     loadings$PC2*4,
     labels=loadings$element,
     col="blue")

dev.off()

cat("\nPIPELINE COMPLETE\n")
