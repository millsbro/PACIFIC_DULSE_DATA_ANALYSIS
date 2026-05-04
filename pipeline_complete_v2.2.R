# ============================================================
# PACIFIC DULSE — FULL PIPELINE 
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
# !!! USER INPUT REQUIRED !!!
# SELECT DATA MODE BEFORE RUNNING SCRIPT
# ============================================================
cat("\nSELECT DATA MODE:\n")
cat("1: Monthly\n")
cat("2: Seasonal\n")

choice <- readline(prompt = "Enter 1 or 2: ")

mode <- ifelse(choice == "2", "seasonal", "monthly")

method_choice <- if (mode == "seasonal") "lm" else "loess"

cat(paste0("\nRunning in ", toupper(mode), " mode\n\n"))


# ============================================================
# GLOBAL SETTINGS
# ============================================================

# Folder where all plots will be saved
fig_dir <- "figures"

# Create folder if it doesn't exist
dir.create(fig_dir, showWarnings = FALSE)

# Helper function to save ggplot objects
save_fig <- function(filename, plot, w=8, h=5){
  ggsave(
    file.path(fig_dir, filename),
    plot = plot,
    width = w,
    height = h,
    dpi = 300,
    bg = "white"
  )
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
library(grid)       # plotting controls
library(ggrepel)    # plotting controls






# ============================================================
# SECTION 1 — LOAD DATA
# ============================================================

# -------------------------
# LOAD DATASETS 
# -------------------------

data <- read.csv(ifelse(mode=="monthly","MONTHLY.csv","SEASONAL.csv"))

bio <- read.csv(ifelse(mode=="monthly",
                       "bioactives_raw.csv",
                       "bioactives_seasonal.csv"))

minerals <- read.csv(ifelse(mode=="monthly",
                            "minerals_raw.csv",
                            "minerals_seasonal.csv"))

faa_cn <- read.csv(ifelse(mode=="monthly",
                          "faa_cn_raw.csv",
                          "faa_seasonal.csv"))

trends <- read.csv("trends.csv")


# -------------------------
# CLEAN COLUMN NAMES
# -------------------------

clean_names <- function(df){
  names(df) <- names(df) %>%
    tolower() %>%
    gsub("%","pct",.) %>%
    gsub("/","_",.) %>%
    gsub("-","_",.) %>%
    gsub(" ","_",.)
  df
}

data     <- clean_names(data)
bio      <- clean_names(bio)
minerals <- clean_names(minerals)
faa_cn   <- clean_names(faa_cn)
trends   <- clean_names(trends)

if(mode == "seasonal"){

  bio$month <- dplyr::case_when(
    bio$month %in% c("DEC","JAN","FEB") ~ "Winter",
    bio$month %in% c("MAR","APR","MAY") ~ "Spring",
    bio$month %in% c("JUN","JUL","AUG") ~ "Summer",
    bio$month %in% c("SEP","OCT","NOV") ~ "Fall"
  )

}


# -------------------------
# STANDARDIZE TIME FACTOR
# -------------------------

if(mode == "monthly"){

  time_levels <- c("DEC","JAN","FEB","MAR","APR","MAY",
                   "JUN","JUL","AUG","SEP","OCT","NOV")

} else {

  time_levels <- c("Winter","Spring","Summer","Fall")

}

data$month     <- factor(data$month, levels = time_levels)
bio$month      <- factor(bio$month, levels = time_levels)
minerals$month <- factor(minerals$month, levels = time_levels)
faa_cn$month   <- factor(faa_cn$month, levels = time_levels)


# ================================
# FAA PCA (creates aa_pc1)
# ================================

faa_cols <- c(
  "glycine_mgkg","glutamine_mg.kg","serine_mgkg","asparagine_mgkg",
  "glutamicacid_mgkg","asparticacid_mgkg","histidine_mgkg","arginine_mgkg",
  "lysine_mgkg","cystine_mgkg","tryptophan_mgkg","phenylalanine_mgkg",
  "leucine_mgkg","isoleucine_mgkg","methionine_mgkg","valine_mgkg",
  "proline_mgkg","tyrosine_mgkg","cysteine_mgkg","alanine_mgkg",
  "threonine_mgkg","taurine_mgkg"
)

faa_data <- data[, faa_cols]

# Remove zero-variance columns
faa_data <- faa_data[, apply(faa_data, 2, var, na.rm = TRUE) != 0]

# Keep track of complete rows
complete_idx <- complete.cases(faa_data)

# Scale + PCA
faa_scaled <- scale(faa_data[complete_idx, ])
aa_pca <- prcomp(faa_scaled, center = TRUE, scale. = TRUE)

# Initialize + assign PC1
data$aa_pc1 <- NA
data$aa_pc1[complete_idx] <- aa_pca$x[,1]




# ============================================================
# SECTION 2 — MACROALGAE INDUSTRY TRENDS 
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
  
geom_smooth(method = method_choice, se = FALSE, linewidth = 1.5) +
  
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






# ============================================================
# SECTION 3 — ABIOTIC PLOTS
# ============================================================

plot_abiotic <- function(df, var, ylab, title, fname){

  # Calculate annual mean
  m <- mean(df[[var]], na.rm=TRUE)

  # Classify points relative to mean
  df$group <- ifelse(df[[var]]>m,"Above mean","Below mean")

  # Build plot
  p <- ggplot(df, aes(month, .data[[var]], group=1))+
    geom_line(color="black", na.rm = TRUE)+
    geom_point(aes(color=group), size=3, na.rm = TRUE)+
    geom_hline(yintercept=m, linetype="dashed")+
    scale_color_manual(values=c("red","blue"))+
    labs(title=title, y=ylab, x="Month")+
    theme_classic(base_size=14)+
    theme(
      plot.title = element_text(hjust=0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  # Save
  save_fig(fname,p)
}

plot_abiotic(data,"temp_c","Temperature (°C)","Temperature","abiotic_temp.png")
plot_abiotic(data,"par_mol_photons_m2","PAR","PAR","abiotic_par.png")
plot_abiotic(data,"salinity_ppt","Salinity","Salinity","abiotic_sal.png")
plot_abiotic(data,"do_mgl","DO","Dissolved Oxygen","abiotic_do.png")
plot_abiotic(data,"ph","pH","pH","abiotic_ph.png")







# ============================================================
# SECTION 4 — MACRONUTRIENTS
# ============================================================

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
# SECTION 5 — AMINO ACIDS 
# ============================================================

# HEAT MAP
# ============================================================
aa <- data %>%
  select(month, ends_with("mgkg")) %>%
  column_to_rownames("month") %>%
  as.matrix()

aa <- as.matrix(aa)
storage.mode(aa) <- "numeric"

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


# PCA — CLEAN (TOP LOADINGS ONLY)
# ============================================================

pca <- prcomp(aa_scaled, center = FALSE, scale. = FALSE)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

loadings <- as.data.frame(pca$rotation)
loadings$var <- rownames(loadings)

loadings_aa <- loadings %>%
  mutate(importance = abs(PC1) + abs(PC2))


# SELECT TOP CONTRIBUTORS
# -------------------------

top_n=8

if(mode == "monthly"){

  loadings_plot_aa <- loadings_aa %>%
    arrange(desc(importance)) %>%
    slice(1:top_n)

} else {

  loadings_plot_aa <- loadings_aa

}


# SCALE ARROWS
# -------------------------

arrow_scale <- 3


# PLOT
# -------------------------

p <- ggplot(scores, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, color = "black") +
    geom_text_repel(aes(label = month), size = 4) +

    geom_segment(
        data = loadings_plot_aa,
        aes(x = 0, y = 0,
            xend = PC1 * arrow_scale,
            yend = PC2 * arrow_scale),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "red"
    ) +

    geom_text_repel(
        data = loadings_plot_aa,
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
# SECTION 7 — FAA vs C:N
# ============================================================

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
# SECTION 8 — MINERALS 
# ============================================================

mat <- minerals %>%
  column_to_rownames("month") %>%
  as.matrix()

mat <- as.matrix(mat)
storage.mode(mat) <- "numeric"

mat_scaled <- apply(mat,2,scale_safe)
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[is.infinite(mat_scaled)] <- 0

pheatmap(mat_scaled,
         filename=file.path(fig_dir,"mineral_heatmap.png"))


# PCA 
# ============================================================

pca_m <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)

scores_m <- as.data.frame(pca_m$x)
scores_m$month <- rownames(scores_m)

loadings_m <- as.data.frame(pca_m$rotation)
loadings_m$var <- rownames(loadings_m)


# SELECT TOP CONTRIBUTORS
# -------------------------

top_n=8

loadings_top_m <- loadings_m %>%
  arrange(desc(importance)) %>%
  slice(1:top_n)

if(mode == "monthly"){

  loadings_plot_m <- loadings_m %>%
    arrange(desc(loadings_m$importance)) %>%
    slice(1:top_n)

} else {

  loadings_plot_m <- loadings_m

}


# SCALE ARROWS
# -------------------------

arrow_scale <- 3


# PLOT
# -------------------------

p_m <- ggplot(scores_m, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, color = "black") +
    geom_text_repel(aes(label = month), size = 4) +

    geom_segment(
        data = loadings_plot_m,
        aes(x = 0, y = 0,
            xend = PC1 * arrow_scale,
            yend = PC2 * arrow_scale),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "red"
    ) +

    geom_text_repel(
        data = loadings_plot_m,
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


# TIME SERIES 
# ============================================================

mineral_numeric <- minerals %>%
  select(-month) %>%
  mutate(across(everything(), as.numeric))

mineral_totals <- data.frame(
  month = minerals$month,
  total_mgkg = rowSums(mineral_numeric, na.rm = TRUE)
)

mineral_totals$month <- factor(mineral_totals$month, levels = time_levels)

p <- ggplot(mineral_totals, aes(x = month, y = total_mgkg, group = 1)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = method_choice, se = FALSE, color = "red", linewidth = 1.2) +
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







# ============================================================
# SECTION 9 — BIOACTIVES 
# ============================================================

bio <- clean_names(bio)

bio$month <- factor(bio$month, levels = time_levels)

bio <- bio %>%
  arrange(month)

bio <- bio %>%
  mutate(
    vitamin_e_mg100g = alpha_t_mg100g + gamma_t_mg100g
  )

vitamins <- c("vitamin_e_mg100g", "k1_ug100g")
pigments <- c("carotenoids_mg100g", "phycoerythrin_g100g")
antiox   <- c("phenolics_mg100g", "tac_mmolkg")
sulfpoly <- c("sulfated_polysaccharides_mggdw")

bio_scaled_df <- bio %>%
  dplyr::select(month, all_of(c(vitamins, pigments, antiox, sulfpoly))) %>%
  mutate(across(-month, scale)) %>%
  as.data.frame()

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

# SULF POLY PANEL 
# -------------------------

plot_sulf <- function(){
  df <- bio_scaled_df %>%
    select(month, all_of(sulfpoly)) %>%
    pivot_longer(-month, names_to="compound", values_to="z")

  ggplot(df, aes(month, z, group = compound)) +
    geom_line(linewidth = 1, color = "purple") +
    geom_point(size = 2, color = "purple") +
    labs(title = "Sulfated Polysaccharides", x = "", y = "Z-score") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
}

# EXPORT PANEL FIGURE
# -------------------------

png(file.path(fig_dir, "bioactives_temporal.png"), width = 1000, height = 1400, bg = "white")

grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(4,1)))

print(plot_group(vitamins, "Vitamins"),
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))

print(plot_group(pigments, "Pigments"),
      vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))

print(plot_group(antiox, "Antioxidants"),
      vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))

print(plot_sulf(),
      vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1))

dev.off()


# HEATMAP 
# -------------------------

mat <- bio_scaled_df %>%
  select(-month) %>%
  as.matrix()

rownames(mat) <- paste0(bio_scaled_df$month, "_", seq_len(nrow(mat)))

pheatmap(
  mat,
  filename = file.path(fig_dir, "bioactives_heatmap.png"),
  main = "Bioactive Compound Profiles (Z-scored)"
)







# ============================================================
# SECTION 10 — PEARSON CORRELATION 
# ============================================================

# -------------------------
# BASE DATA (ABIOTIC + MACRO + N)
# -------------------------

corr_base <- data %>%
  select(
    month,
    temp_c,
    salinity_ppt,
    do_mgl,
    par_mol_photons_m2,
    ph,
    cn_ratio,
    protein_.,
    lipid_.,
    carb_.,
    ash_.,
    moisture_.
  )

# -------------------------
# BIOACTIVES 
# -------------------------

corr_bio <- bio %>%
  select(
    month,
    phenolics_mg100g,
    tac_mmolkg,
    carotenoids_mg100g,
    phycoerythrin_g100g,
    vitamin_e_mg100g,
    k1_ug100g,
    sulfated_polysaccharides_mggdw
  )

# -------------------------
# CALL AMINO ACID PCA (PC1)
# -------------------------

aa_scores <- data %>%
  select(month, aa_pc1)

# -------------------------
# MERGE EVERYTHING
# -------------------------

corr_data <- corr_base %>%
  left_join(corr_bio, by = "month") %>%
  left_join(aa_scores, by = "month")

# -------------------------
# NUMERIC ONLY
# -------------------------

corr_numeric <- corr_data %>%
  select(where(is.numeric))

# -------------------------
# CORRELATION MATRIX
# -------------------------

corr_matrix <- cor(corr_numeric, use = "pairwise.complete.obs")

#print(round(corr_matrix, 3))

write.csv(
  corr_matrix,
  file.path(fig_dir, "pearson_correlation_matrix.csv"),
  row.names = TRUE
)

# -------------------------
# HEATMAP 
# -------------------------

pheatmap(
  corr_matrix,
  filename = file.path(fig_dir, "pearson_heatmap.png"),
  main = "Pearson Correlation Matrix (Full System)",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  na_col = "white"
)

cat("\n--- PEARSON TEST: SULFATED POLYSACCHARIDES vs C:N ---\n\n")

x <- corr_numeric$sulfated_polysaccharides_mggdw
y <- corr_numeric$cn_ratio

valid <- is.finite(x) & is.finite(y)

if(sum(valid) >= 3){
  print(cor.test(x[valid], y[valid]))
} else {
  message("Pearson test skipped: not enough valid points")
}






# ============================================================
# SECTION 11 — FAA vs C:N (FINAL)
# ============================================================

# -------------------------
# IDENTIFY AA COLUMNS (AA ONLY)
# -------------------------

aa_cols <- names(data)[grepl("_mgkg$", names(data)) & 
                       !grepl("^(ca|fe|zn|mg|na|k|p|s|mn|cu|mo|b|as|cr|cd|co|ni|pb)_", names(data))]

# -------------------------
# CALCULATE TOTAL AA
# -------------------------

total_aa <- rowSums(data[, aa_cols], na.rm = TRUE)

# -------------------------
# BUILD DATASET
# -------------------------

plot_df <- data.frame(
  faa = total_aa,
  cn  = data$cn_ratio
)

plot_df <- na.omit(plot_df)

# -------------------------
# CORRELATION
# -------------------------

cor_res <- cor.test(plot_df$faa, plot_df$cn)

r_val <- round(cor_res$estimate, 2)
p_val <- round(cor_res$p.value, 3)

# -------------------------
# PLOT
# -------------------------

p <- ggplot(plot_df, aes(x = cn, y = faa)) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  annotate(
    "text",
    x = max(plot_df$cn),
    y = max(plot_df$faa),
    label = paste0("r = ", r_val, "\np = ", p_val),
    hjust = 1,
    vjust = 1,
    size = 5
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "FAA vs C:N",
    x = "C:N ratio",
    y = "Free Amino Acids (mg/kg)"
  )

# -------------------------
# SAVE
# -------------------------

save_fig("faa_vs_cn_final.png", p)

# -------------------------
# PRINT RESULT
# -------------------------

cat("\n--- PEARSON TEST: TOTAL AA (FAA PROXY) vs C:N ---\n\n")

print(cor_res)




cat(paste0("\nPIPELINE COMPLETE v2.2 (", toupper(mode), " MODE)\n"))
cat("\nALL OUTPUTS SAVED TO /figures\n") 
