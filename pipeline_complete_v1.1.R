# ============================================================
# PACIFIC DULSE — FULL PIPELINE (ANNOTATED VERSION)
# ============================================================
# This script is the canonical, reproducible workflow for:
# - Abiotic analysis
# - FAA vs C:N relationships
# - Macronutrient composition
# - Amino acid profiling (heatmap + PCA)
# - Mineral profiling (heatmap + PCA)
#
# DESIGN PRINCIPLES:
# 1. Fully reproducible (no hidden dependencies)
# 2. Robust to messy column names
# 3. Fail-safe scaling (no NA/Inf crashes)
# 4. All outputs saved automatically to /figures
# ============================================================


# ============================================================
# GLOBAL SETTINGS
# ============================================================

# Folder where all plots will be saved
fig_dir <- "figures"

# Create folder if it doesn't exist
dir.create(fig_dir, showWarnings = FALSE)

# Helper function to save ggplot objects
save_fig <- function(filename, plot, w=8, h=5){
  ggsave(file.path(fig_dir, filename), plot, width=w, height=h, dpi=300)
}


# ============================================================
# LIBRARIES
# ============================================================

library(dplyr)      # data manipulation
library(ggplot2)    # plotting
library(pheatmap)   # heatmaps
library(tidyr)      # pivoting
library(tibble)     # rownames handling
library(stringr)    # pattern detection
library(ggrepel)
library(grid)

# ============================================================
# GLOBAL PLOT THEME (FORCE WHITE BACKGROUND)
# ============================================================

theme_set(
  theme_classic(base_size = 14) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    )
)


# ============================================================
# SECTION 1 — LOAD DATA
# ============================================================

# Toggle between datasets
mode <- "monthly"

# Load datasets
data <- read.csv(ifelse(mode=="monthly","MONTHLY.csv","SEASONAL.csv"))
faa_cn <- read.csv("faa_cn_raw.csv")
minerals <- read.csv("minerals_raw.csv")


# ============================================================
# CLEAN COLUMN NAMES
# ============================================================
# Ensures consistent naming across all datasets
# Handles %, spaces, slashes, etc.

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


# ============================================================
# STANDARDIZE MONTH ORDER
# ============================================================

month_levels <- c("DEC","JAN","FEB","MAR","APR","MAY",
                  "JUN","JUL","AUG","SEP","OCT","NOV")

data$month <- factor(data$month, levels=month_levels)
faa_cn$month <- factor(faa_cn$month, levels=month_levels)
minerals$month <- factor(minerals$month, levels=month_levels)


# ============================================================
# SECTION 2 — ABIOTIC PLOTS
# ============================================================
# Creates 5 plots (temp, PAR, salinity, DO, pH)
# Each includes:
# - mean line
# - above/below mean coloring

plot_abiotic <- function(df, var, ylab, title, fname){

  # Calculate annual mean
  m <- mean(df[[var]], na.rm=TRUE)

  # Classify points relative to mean
  df$group <- ifelse(df[[var]]>m,"Above mean","Below mean")

  # Build plot
  p <- ggplot(df, aes(month, .data[[var]], group=1))+
    geom_line(color="black")+
    geom_point(aes(color=group), size=3)+
    geom_hline(yintercept=m, linetype="dashed")+
    scale_color_manual(values=c("red","blue"))+
    labs(title=title, y=ylab, x="Month")+
    theme_minimal(base_size=14)+
    theme(plot.title=element_text(hjust=0.5))

  # Save
  save_fig(fname,p)
}

plot_abiotic(data,"temp_c","Temperature (°C)","Temperature","abiotic_temp.png")
plot_abiotic(data,"par_mol_photons_m2","PAR","PAR","abiotic_par.png")
plot_abiotic(data,"salinity_ppt","Salinity","Salinity","abiotic_sal.png")
plot_abiotic(data,"do_mgl","DO","Dissolved Oxygen","abiotic_do.png")
plot_abiotic(data,"ph","pH","pH","abiotic_ph.png")


# ============================================================
# SECTION 3 — FAA vs C:N
# ============================================================
# Uses separate FAA dataset (correct scientific separation)

faa <- faa_cn %>%
  rename(
    faa = faa_mg_kg,
    cn  = c_to_n
  )

p <- ggplot(faa, aes(cn,faa))+
  geom_point(size=3)+
  geom_smooth(method="lm", se=TRUE, color="blue")+
  labs(title="FAA vs C:N",
       x="C:N ratio",
       y="FAA (mg/kg)")+
  theme_minimal()

save_fig("faa_vs_cn.png",p)


# ============================================================
# SECTION 4 — MACRONUTRIENTS
# ============================================================
# Automatically detects protein/lipid/carb/ash columns

protein_col <- names(data)[str_detect(names(data),"protein")]
lipid_col   <- names(data)[str_detect(names(data),"lipid")]
carb_col    <- names(data)[str_detect(names(data),"carb")]
ash_col     <- names(data)[str_detect(names(data),"ash")]

macro_long <- data %>%
  select(month,
         protein = all_of(protein_col[1]),
         lipid   = all_of(lipid_col[1]),
         carb    = all_of(carb_col[1]),
         ash     = all_of(ash_col[1])) %>%
  pivot_longer(-month)

p <- ggplot(macro_long, aes(month,value,fill=name))+
  geom_bar(stat="identity")+
  labs(title="Macronutrient Composition")

save_fig("macros_stacked.png",p)


# ============================================================
# SECTION 5 — AMINO ACIDS (UPDATED PCA)
# ============================================================

aa <- data %>%
  select(month, ends_with("mgkg")) %>%
  column_to_rownames("month") %>%
  as.matrix()

aa <- apply(aa,2,as.numeric)

scale_safe <- function(x){
  s <- sd(x, na.rm=TRUE)
  if(is.na(s) | s == 0) return(rep(0,length(x)))
  (x-mean(x, na.rm=TRUE))/s
}

aa_scaled <- apply(aa,2,scale_safe)
aa_scaled[is.na(aa_scaled)] <- 0
aa_scaled[is.infinite(aa_scaled)] <- 0

pheatmap(aa_scaled,
         filename=file.path(fig_dir,"aa_heatmap.png"))


# ============================================================
# PCA — CLEAN (TOP LOADINGS ONLY)
# ============================================================

pca <- prcomp(aa_scaled, center = FALSE, scale. = FALSE)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

loadings <- as.data.frame(pca$rotation)
loadings$var <- rownames(loadings)

# -------------------------
# SELECT TOP CONTRIBUTORS
# -------------------------

top_n <- 8

loadings$importance <- abs(loadings$PC1) + abs(loadings$PC2)

loadings_top <- loadings %>%
    arrange(desc(importance)) %>%
    slice(1:top_n)

# -------------------------
# SCALE ARROWS
# -------------------------

arrow_scale <- 3

# -------------------------
# PLOT
# -------------------------

p <- ggplot(scores, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, color = "black") +
    geom_text_repel(aes(label = month), size = 4) +

    geom_segment(
        data = loadings_top,
        aes(x = 0, y = 0,
            xend = PC1 * arrow_scale,
            yend = PC2 * arrow_scale),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "red"
    ) +

    geom_text_repel(
        data = loadings_top,
        aes(x = PC1 * arrow_scale,
            y = PC2 * arrow_scale,
            label = var),
        color = "blue",
        size = 4
    ) +

    theme_minimal() +
    labs(
        title = "Amino Acid PCA",
        x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)")
    )


save_fig("aa_pca.png", p)



# ============================================================
# SECTION 6 — EAA vs NEAA
# ============================================================

eaa <- c("lysine_mgkg","leucine_mgkg","isoleucine_mgkg","valine_mgkg",
         "methionine_mgkg","phenylalanine_mgkg","tryptophan_mgkg","threonine_mgkg")

aa_eaa <- rowSums(data[,eaa], na.rm=TRUE)
aa_total <- rowSums(aa, na.rm=TRUE)
aa_neaa <- aa_total - aa_eaa

df <- data.frame(month=data$month,EAA=aa_eaa,NEAA=aa_neaa) %>%
  pivot_longer(-month)

p <- ggplot(df,aes(month,value,fill=name))+
  geom_bar(stat="identity")+
  labs(title="EAA vs NEAA")

save_fig("aa_eaa_neaa.png",p)


# ============================================================
# SECTION 7 — MINERALS (UPDATED PCA)
# ============================================================

mat <- minerals %>%
  column_to_rownames("month") %>%
  as.matrix()

mat <- apply(mat,2,as.numeric)

mat_scaled <- apply(mat,2,scale_safe)
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[is.infinite(mat_scaled)] <- 0

pheatmap(mat_scaled,
         filename=file.path(fig_dir,"mineral_heatmap.png"))

# ============================================================
# PCA — MINERALS (CLEAN)
# ============================================================

pca_m <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)

scores_m <- as.data.frame(pca_m$x)
scores_m$month <- rownames(scores_m)

loadings_m <- as.data.frame(pca_m$rotation)
loadings_m$var <- rownames(loadings_m)

# -------------------------
# SELECT TOP CONTRIBUTORS
# -------------------------

top_n <- 8

loadings_m$importance <- abs(loadings_m$PC1) + abs(loadings_m$PC2)

loadings_top_m <- loadings_m %>%
    arrange(desc(importance)) %>%
    slice(1:top_n)

# -------------------------
# SCALE ARROWS
# -------------------------

arrow_scale <- 3

# -------------------------
# PLOT
# -------------------------

p_m <- ggplot(scores_m, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, color = "black") +
    geom_text_repel(aes(label = month), size = 4) +

    geom_segment(
        data = loadings_top_m,
        aes(x = 0, y = 0,
            xend = PC1 * arrow_scale,
            yend = PC2 * arrow_scale),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "red"
    ) +

    geom_text_repel(
        data = loadings_top_m,
        aes(x = PC1 * arrow_scale,
            y = PC2 * arrow_scale,
            label = var),
        color = "blue",
        size = 4
    ) +

    theme_minimal() +
    labs(
        title = "Mineral PCA",
        x = paste0("PC1 (", round(summary(pca_m)$importance[2,1]*100,1), "%)"),
        y = paste0("PC2 (", round(summary(pca_m)$importance[2,2]*100,1), "%)")
    )

save_fig("mineral_pca.png", p_m)

# ============================================================
# SECTION 8 — MINERAL TIME SERIES (FIXED)
# ============================================================

mineral_numeric <- minerals %>%
  select(-month) %>%
  mutate(across(everything(), as.numeric))

mineral_totals <- data.frame(
  month = minerals$month,
  total_mgkg = rowSums(mineral_numeric, na.rm = TRUE)
)

mineral_totals$month <- factor(mineral_totals$month, levels = month_levels)

p <- ggplot(mineral_totals, aes(x = month, y = total_mgkg, group = 1)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1.2) +
  labs(
    title = "Total Mineral Content",
    x = "Month",
    y = "Total Mineral Content (mg/kg)"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(
    limits = c(
      floor(min(mineral_totals$total_mgkg) / 1000) * 1000,
      ceiling(max(mineral_totals$total_mgkg) / 1000) * 1000
    ),
    breaks = scales::pretty_breaks(n = 6)
  )

save_fig("mineral_timeseries.png", p)

cat("\nSECTION 8 COMPLETE\n")



# ============================================================
# SECTION 9 — BIOACTIVES (FINAL)
# ============================================================

bio <- read.csv("bioactives_raw.csv")
bio <- clean_names(bio)

bio$month <- factor(bio$month, levels = month_levels)

bio <- bio %>%
  mutate(
    vitamin_e_mg100g = alpha_t_mg100g + gamma_t_mg100g
  )

vitamins <- c("vitamin_e_mg100g", "k1_ug100g")
pigments <- c("carotenoids_mg100g", "phycoerytherin_g100g")
antiox   <- c("phenolics_mg100g", "tac_mmolkg")

bio_numeric <- bio %>%
  select(month, all_of(c(vitamins, pigments, antiox)))

bio_mat <- bio_numeric %>%
  column_to_rownames("month") %>%
  as.matrix()

bio_mat <- apply(bio_mat, 2, as.numeric)

bio_scaled <- apply(bio_mat, 2, scale_safe)
bio_scaled[is.na(bio_scaled)] <- 0
bio_scaled[is.infinite(bio_scaled)] <- 0

bio_scaled_df <- as.data.frame(bio_scaled)
bio_scaled_df$month <- rownames(bio_scaled_df)

plot_group <- function(vars, title){
  df <- bio_scaled_df %>%
    select(month, all_of(vars)) %>%
    pivot_longer(-month, names_to="compound", values_to="z")

  ggplot(df, aes(month, z, group = compound, color = compound)) +
    geom_line(linewidth = 1) +
    labs(title = title, x = "", y = "Z-score") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    )
}

png(file.path(fig_dir, "bioactives_temporal.png"), width = 1000, height = 1200)

grid::grid.newpage()
pushViewport(grid::viewport(layout = grid::grid.layout(3,1)))

print(plot_group(vitamins, "Vitamins"),
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))

print(plot_group(pigments, "Pigments"),
      vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))

print(plot_group(antiox, "Antioxidants"),
      vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))

dev.off()

pheatmap(
  bio_scaled,
  filename = file.path(fig_dir, "bioactives_heatmap.png"),
  main = "Bioactive Compound Profiles (Z-scored)"
)


# ============================================================
# SECTION 10 — PEARSON CORRELATION
# ============================================================

corr_data <- data %>%
  select(month,
         temp_c,
         salinity_ppt,
         do_mgl,
         par_mol_photons_m2,
         ph)

faa_cn_df <- faa_cn %>%
  select(month, faa_mg_kg, c_to_n) %>%
  rename(faa = faa_mg_kg, cn = c_to_n)

corr_data <- corr_data %>%
  left_join(faa_cn_df, by = "month")

corr_data <- corr_data %>%
  left_join(mineral_totals, by = "month")

bio_corr <- bio %>%
  select(month,
         vitamin_e_mg100g,
         k1_ug100g,
         carotenoids_mg100g,
         phycoerytherin_g100g,
         phenolics_mg100g,
         tac_mmolkg)

corr_data <- corr_data %>%
  left_join(bio_corr, by = "month")

corr_numeric <- corr_data %>%
  select(-month) %>%
  mutate(across(everything(), as.numeric))

cor_matrix <- cor(corr_numeric,
                  use = "pairwise.complete.obs",
                  method = "pearson")

dir.create("outputs", showWarnings = FALSE)

write.csv(cor_matrix,
          file.path("outputs", "pearson_matrix.csv"),
          row.names = TRUE)

pheatmap(
  cor_matrix,
  filename = file.path(fig_dir, "pearson_heatmap.png"),
  main = "Pearson Correlation Matrix",
  color = colorRampPalette(c("blue","white","red"))(100)
)

cat("\n--- FAA vs C:N ---\n")
test <- cor.test(corr_numeric$faa, corr_numeric$cn)
print(test)


# ============================================================
# SECTION 11 — FAA vs C:N (FINAL)
# ============================================================

plot_df <- data.frame(
  faa = corr_numeric$faa,
  cn  = corr_numeric$cn
)

plot_df <- plot_df[complete.cases(plot_df), ]

test <- cor.test(plot_df$faa, plot_df$cn)

r_val <- round(test$estimate, 2)
p_val <- signif(test$p.value, 2)

label_text <- paste0("r = ", r_val, "\n", "p = ", p_val)

p <- ggplot(plot_df, aes(x = cn, y = faa)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "FAA vs C:N",
       x = "C:N ratio",
       y = "Free Amino Acids (mg/kg)") +
  annotate("text",
           x = max(plot_df$cn, na.rm = TRUE),
           y = max(plot_df$faa, na.rm = TRUE),
           label = label_text,
           hjust = 1,
           vjust = 1,
           size = 5) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

save_fig("faa_vs_cn_final.png", p)

cat("\nPIPELINE COMPLETE v1.0\n")


# ============================================================
# SECTION 12 — MACROALGAE INDUSTRY TRENDS (FINAL)
# ============================================================

# -------------------------
# LOAD + CLEAN
# -------------------------

trends <- read.csv("trends.csv")
trends <- clean_names(trends)

# -------------------------
# RENAME COLUMNS (EXACT MATCH)
# -------------------------

trends <- trends %>%
  rename(
    food        = `seaweed.food`,
    agriculture = `x.seaweed.fertilizer`,
    biofuel     = `x.algae.biofuel`,
    cosmetics   = `x.seaweed.cosmetics`,
    pharma      = `x.algae.pharmaceutical`
  )

# -------------------------
# FORMAT TIME
# -------------------------

trends$time <- as.Date(trends$time)
trends$year <- as.numeric(format(trends$time, "%Y"))

# -------------------------
# AGGREGATE TO YEARLY MEAN
# -------------------------

trends_yearly <- trends %>%
  group_by(year) %>%
  summarise(
    across(food:pharma, mean, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------
# LONG FORMAT
# -------------------------

trends_long <- trends_yearly %>%
  pivot_longer(
    cols = -year,
    names_to = "sector",
    values_to = "interest"
  )

# -------------------------
# PLOT (SMOOTHED — FINAL)
# -------------------------

p <- ggplot(trends_long, aes(x = year, y = interest, color = sector)) +
  
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.5, span = 0.4) +
  
  scale_color_manual(
    values = c(
      food = "black",
      agriculture = "#2E8B57",
      biofuel = "#E74C3C",
      cosmetics = "#8E44AD",
      pharma = "#F39C12"
    ),
    labels = c(
      food = "Food",
      agriculture = "Agriculture",
      biofuel = "Biofuel",
      cosmetics = "Cosmetics",
      pharma = "Pharma"
    )
  ) +
  
  labs(
    title = "Trends in Interest in Macroalgal Applications by Industry",
    x = "Year",
    y = "Relative Search Interest",
    color = ""
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

# -------------------------
# SAVE
# -------------------------

save_fig("macroalgae_trends.png", p)

cat("\nPIPELINE COMPLETE v1.1\n")
