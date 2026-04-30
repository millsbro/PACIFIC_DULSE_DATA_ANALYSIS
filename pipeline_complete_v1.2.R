# ============================================================
# PACIFIC DULSE — FULL PIPELINE (v1.2)
# ============================================================

fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

save_fig <- function(filename, plot, w=8, h=5){
  ggsave(file.path(fig_dir, filename), plot, width=w, height=h, dpi=300)
}

library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(tibble)
library(stringr)
library(grid)

# ============================================================
# LOAD DATA
# ============================================================

mode <- "monthly"

data <- read.csv(ifelse(mode=="monthly","MONTHLY.csv","SEASONAL.csv"))
faa_cn <- read.csv("faa_cn_raw.csv")
minerals <- read.csv("minerals_raw.csv")

clean_names <- function(df){
  names(df) <- names(df) %>%
    tolower() %>%
    gsub("%","pct",.) %>%
    gsub("/","_",.) %>%
    gsub("-","_",.) %>%
    gsub(" ","_",.)
  df
}

data <- clean_names(data)
faa_cn <- clean_names(faa_cn)
minerals <- clean_names(minerals)

month_levels <- c("DEC","JAN","FEB","MAR","APR","MAY",
                  "JUN","JUL","AUG","SEP","OCT","NOV")

data$month <- factor(data$month, levels=month_levels)
faa_cn$month <- factor(faa_cn$month, levels=month_levels)
minerals$month <- factor(minerals$month, levels=month_levels)

# ============================================================
# PCA SAFE SCALING
# ============================================================

scale_safe <- function(x){
  s <- sd(x, na.rm=TRUE)
  if(is.na(s) | s == 0) return(rep(0,length(x)))
  (x-mean(x, na.rm=TRUE))/s
}

# ============================================================
# SECTION 5 — AMINO ACID PCA (FIXED)
# ============================================================

aa <- data %>%
  select(month, contains("mgkg")) %>%
  column_to_rownames("month") %>%
  as.matrix()

aa <- apply(aa,2,as.numeric)

aa_scaled <- apply(aa,2,scale_safe)
aa_scaled[is.na(aa_scaled)] <- 0

pheatmap(aa_scaled, filename=file.path(fig_dir,"aa_heatmap.png"))

pca <- prcomp(aa_scaled)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

loadings <- as.data.frame(pca$rotation)
loadings$var <- rownames(loadings)

var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

p <- ggplot() +
  geom_point(data = scores, aes(PC1, PC2), size = 3) +
  geom_text(data = scores, aes(PC1, PC2, label = month), vjust = -1) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "red") +
  geom_text(data = loadings,
            aes(PC1*3, PC2*3, label = var),
            color = "blue", size = 3) +
  labs(
    title = "Amino Acid PCA",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_classic()

save_fig("aa_pca.png", p)

# ============================================================
# SECTION 7 — MINERAL PCA (FIXED)
# ============================================================

mat <- minerals %>%
  column_to_rownames("month") %>%
  as.matrix()

mat <- apply(mat,2,as.numeric)

mat_scaled <- apply(mat,2,scale_safe)
mat_scaled[is.na(mat_scaled)] <- 0

pheatmap(mat_scaled, filename=file.path(fig_dir,"mineral_heatmap.png"))

pca <- prcomp(mat_scaled)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

loadings <- as.data.frame(pca$rotation)
loadings$var <- rownames(loadings)

var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

p <- ggplot() +
  geom_point(data = scores, aes(PC1, PC2), size = 3) +
  geom_text(data = scores, aes(PC1, PC2, label = month), vjust = -1) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "red") +
  geom_text(data = loadings,
            aes(PC1*3, PC2*3, label = var),
            color = "blue", size = 3) +
  labs(
    title = "Mineral PCA",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_classic()

save_fig("mineral_pca.png", p)

cat("\nPIPELINE COMPLETE v1.2\n")
