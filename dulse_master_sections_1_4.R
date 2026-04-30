# ============================================================
# SECTION 1 — SETUP
# ============================================================

library(tidyverse)
library(ggplot2)
library(dplyr)

dir.create("figures", showWarnings = FALSE)

theme_clean <- function() {
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey80", size = 0.6),
      panel.grid.minor = element_line(color = "grey90", size = 0.3),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    )
}

# ============================================================
# SECTION 2 — DATA INPUT
# ============================================================

# data <- seasonal
data <- monthly

data <- data %>%
  rename(
    temp = temp_c,
    sal  = salinity_ppt,
    do   = do_mgl,
    par  = par_mol_photons_m2
  )

# ============================================================
# SECTION 3 — ABIOTIC FIGURES
# ============================================================

plot_above_below <- function(df, var, ylab, title, filename) {
  mean_val <- mean(df[[var]], na.rm = TRUE)
  df <- df %>%
    mutate(group = ifelse(.data[[var]] > mean_val, "Above mean", "Below mean"))
  p <- ggplot(df, aes(x = month, y = .data[[var]], group = 1)) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(aes(color = group), size = 3) +
    geom_hline(yintercept = mean_val, linetype = "dashed") +
    scale_color_manual(values = c("Above mean" = "red", "Below mean" = "blue")) +
    labs(title = title, x = "Month", y = ylab, color = "") +
    theme_clean()
  ggsave(paste0("figures/", filename), p, width = 8, height = 5, dpi = 300)
  return(p)
}

p_temp <- plot_above_below(data, "temp","Temperature (°C)","Temperature","abiotic_temp.png")
p_sal  <- plot_above_below(data, "sal","Salinity (ppt)","Salinity","abiotic_sal.png")
p_do   <- plot_above_below(data, "do","Dissolved Oxygen (mg/L)","Dissolved Oxygen","abiotic_do.png")
p_par  <- plot_above_below(data, "par",expression(PAR~(mol~m^{-2}~d^{-1})),"PAR","abiotic_par.png")

mean_ph <- mean(ph$ph, na.rm = TRUE)
ph_plot <- ph %>%
  mutate(group = ifelse(ph > mean_ph, "Above mean", "Below mean"))

p_ph <- ggplot(ph_plot, aes(x = month, y = ph, group = 1)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(aes(color = group), size = 3) +
  geom_errorbar(aes(ymin = ph - sd, ymax = ph + sd),
                width = 0.2, color = "orange") +
  geom_hline(yintercept = mean_ph, linetype = "dashed") +
  scale_color_manual(values = c("Above mean" = "red", "Below mean" = "blue")) +
  labs(title = "pH", x = "Month", y = "pH", color = "") +
  theme_clean()

ggsave("figures/abiotic_ph.png", p_ph, width = 8, height = 5, dpi = 300)

# ============================================================
# SECTION 4 — MACROMOLECULES
# ============================================================

data <- data %>%
  rename(
    protein  = `protein_%`,
    lipid    = `lipid_%`,
    carb     = `carb_%`,
    ash      = `ash_%`,
    moisture = `moisture_%`
  )

macro_long <- data %>%
  select(month, protein, lipid, carb, ash, moisture) %>%
  pivot_longer(cols = -month, names_to = "component", values_to = "value")

macro_long$component <- factor(
  macro_long$component,
  levels = c("carb","protein","lipid","ash","moisture")
)

p_macro <- ggplot(macro_long,
                  aes(x = month, y = value, fill = component)) +
  geom_bar(stat = "identity") +
  labs(title = "Macromolecules", x = "Month", y = "Percent (%)", fill = "") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/macros_stacked.png",
       p_macro, width = 8, height = 5, dpi = 300)

p_temp; p_sal; p_do; p_par; p_ph; p_macro

cat("\nSECTIONS 1–4 COMPLETE\n")
