# ============================================================
# DULSE MASTER SCRIPT — SECTIONS 1–3 (CANONICAL)
# ============================================================

library(tidyverse)
library(ggplot2)
library(dplyr)

# -------------------------
# USER INPUT — SELECT MODE
# -------------------------
# data <- seasonal
data <- monthly

# -------------------------
# STANDARDIZE COLUMN NAMES
# -------------------------
data <- data %>%
  rename(
    temp = temp_c,
    sal  = salinity_ppt,
    do   = do_mgl,
    par  = par_mol_photons_m2
  )

# -------------------------
# COMMON THEME
# -------------------------
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

# -------------------------
# OUTPUT FOLDER
# -------------------------
dir.create("figures", showWarnings = FALSE)

# -------------------------
# HELPER FUNCTION
# -------------------------
plot_above_below <- function(df, var, ylab, title, filename) {

  mean_val <- mean(df[[var]], na.rm = TRUE)

  df <- df %>%
    mutate(group = ifelse(.data[[var]] > mean_val, "Above mean", "Below mean"))

  p <- ggplot(df, aes(x = month, y = .data[[var]], group = 1)) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(aes(color = group), size = 3) +
    geom_hline(yintercept = mean_val, linetype = "dashed") +
    scale_color_manual(values = c("Above mean" = "red", "Below mean" = "blue")) +
    labs(
      title = title,
      x = "Month",
      y = ylab,
      color = ""
    ) +
    theme_clean()

  ggsave(paste0("figures/", filename), p, width = 8, height = 5, dpi = 300)

  return(p)
}

# -------------------------
# ABIOTIC PLOTS
# -------------------------

p_temp <- plot_above_below(
  data, "temp",
  "Temperature (°C)",
  "Seasonal Temperature Relative to Annual Mean",
  "abiotic_temp.png"
)

p_sal <- plot_above_below(
  data, "sal",
  "Salinity (ppt)",
  "Seasonal Salinity Relative to Annual Mean",
  "abiotic_sal.png"
)

p_do <- plot_above_below(
  data, "do",
  "Dissolved Oxygen (mg/L)",
  "Seasonal DO Relative to Annual Mean",
  "abiotic_do.png"
)

p_par <- plot_above_below(
  data, "par",
  expression(PAR~(mol~m^{-2}~d^{-1})),
  "Seasonal PAR Relative to Annual Mean",
  "abiotic_par.png"
)

# -------------------------
# pH (CONSISTENT STYLE)
# -------------------------

mean_ph <- mean(ph$ph, na.rm = TRUE)

ph_plot <- ph %>%
  mutate(group = ifelse(ph > mean_ph, "Above mean", "Below mean"))

p_ph <- ggplot(ph_plot, aes(x = month, y = ph, group = 1)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(aes(color = group), size = 3) +
  geom_errorbar(
    aes(ymin = ph - sd, ymax = ph + sd),
    width = 0.2,
    color = "orange"
  ) +
  geom_hline(yintercept = mean_ph, linetype = "dashed") +
  scale_color_manual(values = c("Above mean" = "red", "Below mean" = "blue")) +
  labs(
    title = "Seasonal pH (High Tide Filtered)",
    x = "Month",
    y = "pH",
    color = ""
  ) +
  theme_clean()

ggsave("figures/abiotic_ph.png", p_ph, width = 8, height = 5, dpi = 300)

# -------------------------
# DISPLAY
# -------------------------

p_temp
p_sal
p_do
p_par
p_ph

cat("\nSECTIONS 1–3 COMPLETE\n")
