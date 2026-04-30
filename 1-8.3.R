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
# SECTION 5 — AMINO ACIDS
# ============================================================
# Includes:
# - Z-scored heatmap
# - PCA

aa <- data %>%
  select(month, contains("mgkg")) %>%
  column_to_rownames("month") %>%
  as.matrix()

aa <- apply(aa,2,as.numeric)

# Safe scaling function (prevents NA/Inf crashes)
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

pca <- prcomp(aa_scaled)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

png(file.path(fig_dir,"aa_pca.png"),1000,800)
plot(scores$PC1,scores$PC2,pch=19)
text(scores$PC1,scores$PC2,labels=scores$month,pos=3)
dev.off()


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
# SECTION 7 — MINERALS
# ============================================================
# Same approach as amino acids

mat <- minerals %>%
  column_to_rownames("month") %>%
  as.matrix()

mat <- apply(mat,2,as.numeric)

mat_scaled <- apply(mat,2,scale_safe)
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[is.infinite(mat_scaled)] <- 0

pheatmap(mat_scaled,
         filename=file.path(fig_dir,"mineral_heatmap.png"))

pca <- prcomp(mat_scaled)

scores <- as.data.frame(pca$x)
scores$month <- rownames(scores)

png(file.path(fig_dir,"mineral_pca.png"),1000,800)
plot(scores$PC1,scores$PC2,pch=19)
text(scores$PC1,scores$PC2,labels=scores$month,pos=3)
dev.off()


# ============================================================
# COMPLETE
# ============================================================

cat("\nPIPELINE COMPLETE (ANNOTATED VERSION)\n")



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
