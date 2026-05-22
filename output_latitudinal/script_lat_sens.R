##################################################################
# Latitudinal sensitivity simulations                            #
# - Large warm-adapted diatom thermal curve & latitudinal runs  #
# - Prochlorococcus alpha (ctag_max) sensitivity                #
# - Comparison: Pro vs Syn / Pro vs warm-adapted diatom         #
# - Relative proteome investments by latitudinal zone           #
##################################################################

rm(list = ls())

library(ggplot2)
library(patchwork)
library(dplyr)
library(zoo)
library(gridExtra)
library(viridis)
library(tidyr)

# Plotting settings ####
mytheme <- theme(panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.background = element_rect(colour = "black", fill = "white", size = 1.0),
                 axis.text = element_text(size = 13, colour = "black"),
                 axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                 axis.text.y = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                 axis.title = element_text(colour = "black"),
                 axis.title.y = element_text(size = 14, margin = unit(c(-0.2, -0.2, -0.2, -0.2), "cm")),
                 axis.title.x = element_text(size = 13, margin = unit(c(-0.2, -0.2, -0.2, -0.2), "cm")),
                 axis.ticks.length = unit(-0.25, "cm"),
                 plot.margin = margin(20, 5, 10, 15),
                 legend.title = element_blank(),
                 legend.direction = "vertical",
                 legend.text = element_text(size = 12, colour = "black"),
                 legend.key = element_blank(),
                 legend.background = element_rect(fill = NA, colour = NA),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 13, color = "black"))


###################################################
# Figure 1: Thermal growth curve + validation (A+B) #
###################################################

dl <- read.csv("output_warm_diatom_thermal_curve_large.csv")  # Large warm-diatom thermal curve
dl$phenotype <- "Warm-adapted diatom LARGE"
dl$DIN2 <- factor(dl$DIN, levels = c("20", "0.1"), labels = c("N-replete", "N-deplete"))

dl2 <- dl %>%
  filter(case_when(
    DIN2 == "N-deplete" ~ NC < 0.3 & NC > 0.04 & size > 0.4 & size < 100  & mu > 1e-3 & mu < 2.0,
    DIN2 == "N-replete" ~ NC < 0.3 & NC > 0.01 & size > 0.4 & size < 2000 & mu > 7e-2 & mu < 2.0
  ))
dl2 <- dl2[-c(74, 80), ]

d1_large <- dl2 %>%
  group_by(DIN2) %>%
  filter(mu == max(mu, na.rm = TRUE))

p5 <- ggplot(dl2, aes(x = temp, y = mu, linetype = DIN2)) +
  geom_line(size = 1) +
  facet_wrap(~ phenotype) +
  scale_y_continuous(breaks = c(0.0, 0.4, 0.8), labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(7, 33)) +
  ylab(expression(paste("Growth rate (", day^-1, ")", sep = ""))) +
  xlab(expression(paste("Temperature  (", degree, "C)"))) +
  mytheme +
  theme(legend.position = "top") +
  geom_point(data = d1_large, aes(x = temp, y = mu), size = 5, color = "#762A83") +
  geom_vline(data = d1_large, aes(xintercept = temp, linetype = DIN2),
             color = "#762A83", size = 1.0, show.legend = FALSE)

p5

tbl_data <- data.frame(
  "Parameter"   = c("μmax (d⁻¹)", "Slope (d⁻¹°C⁻¹)", "Diameter (μm)"),
  "C. sp."      = c("0.31",  "0.05", "-"),
  "C. wailesii" = c("0.76",  "0.03", "-"),
  "Model"       = c("0.479", "0.05", "120"),
  check.names = FALSE
)

italic_cols <- c(FALSE, TRUE, TRUE, FALSE)

col_theme <- ttheme_default(
  base_size = 11,
  core = list(fg_params = list(hjust = 0.5, x = 0.5),
              bg_params = list(fill = "white", col = "black")),
  colhead = list(fg_params = list(hjust = 0.5, x = 0.5,
                                  fontface = ifelse(italic_cols, "italic", "bold")),
                 bg_params = list(fill = "white", col = "black"))
)

tbl_grob <- tableGrob(tbl_data, rows = NULL, theme = col_theme)

p5_AB <- p5 + wrap_elements(tbl_grob) +
  plot_layout(widths = c(1.1, 1.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

p5_AB

ggsave("figure_large_diatom_thermal_growth_AB.png", p5_AB, device = "png", width = 9, height = 4, dpi = 600)

#####################################
# Load and process latitudinal data #
#####################################

cesm <- na.omit(read.csv("model_CESM.csv"))

dlf <- read.csv("output_warm_diatom_latitudinal_future_LARGE.csv")   # Large warm-diatom future
dlp <- read.csv("output_warm_diatom_latitudinal_present_LARGE.csv")  # Large warm-diatom present

dlat <- rbind(dlf, dlp)
dlat$model     <- c(rep("future", nrow(dlf)), rep("present", nrow(dlp)))
dlat$phenotype <- "Large warm-\nadapted diatom"
dlat$latitude  <- rep(cesm$latitude, 2)
dlat$istress   <- dlat$ilpb * dlat$itag

dlat <- dlat %>%
  filter(NC < 0.3, NC > 0.04, size > 5, size < 500, mu > 0.001, mu < 1.0) %>%
  filter(!(model == "present" & latitude == 41.5))

dlat$model <- factor(dlat$model, levels = c("present", "future"))

latitude_seq <- seq(-75.5, 75.5, by = 0.5)

complete_grid <- expand.grid(
  latitude  = latitude_seq,
  model     = unique(dlat$model),
  phenotype = unique(dlat$phenotype)
)

dlat_approx <- complete_grid %>%
  left_join(dlat, by = c("latitude", "model", "phenotype")) %>%
  group_by(model, phenotype) %>%
  mutate(mu      = zoo::na.approx(mu,      na.rm = FALSE),
         size    = zoo::na.approx(size,    na.rm = FALSE),
         itr     = zoo::na.approx(itr,     na.rm = FALSE),
         icbd    = zoo::na.approx(icbd,    na.rm = FALSE),
         ichl    = zoo::na.approx(ichl,    na.rm = FALSE),
         istress = zoo::na.approx(istress, na.rm = FALSE)) %>%
  ungroup() %>%
  mutate(across(.cols = -c(phenotype, model, latitude),
                .fns  = ~ replace(.x, !is.na(mu) & mu == 0, NA_real_)))

#############################################
# Figure 3: Cell size by latitude (future) #
#############################################

selected_latitudes <- round(seq(-70, 70, length.out = 12))

dlat_size5 <- dlat_approx %>%
  filter(latitude %in% selected_latitudes) %>%
  select(phenotype, latitude, size, model) %>%
  group_by(latitude, phenotype) %>%
  mutate(
    size_future  = size[model == "future"][1],
    size_present = size[model == "present"][1],
    size_diff = if (!is.na(size_future[1]) && !is.na(size_present[1])) {
      size_future - size_present
    } else {
      NA_real_
    },
    size_change = case_when(
      size_diff > 0  ~ "increase",
      size_diff < 0  ~ "decrease",
      size_diff == 0 ~ "no change",
      TRUE           ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(-size_future, -size_present)

psize_large <- ggplot(dlat_size5 %>% filter(model == "future"),
                      aes(x = latitude, y = phenotype, size = size, shape = size_change)) +
  geom_point(fill = "#762A83", color = "black") +
  coord_flip() +
  xlab("Latitude") +
  ylab("") +
  scale_x_continuous(breaks = c(-70, -35, 0, 35, 70), labels = scales::number_format(accuracy = 1),
                     limits = c(-70, 70), expand = c(0.08, 0.08)) +
  scale_size_area(breaks = c(45, 80, 120), labels = c("45", "80", "120"), max_size = 10,
                  guide = guide_legend(override.aes = list(shape = 21, fill = NA, color = "black", stroke = 0.8))) +
  scale_shape_manual(values = c("decrease" = 1, "increase" = 19, "no change" = 21), guide = "none") +
  mytheme +
  theme(legend.position = "right",
        plot.margin = margin(40, 5, 10, 15),
        axis.text.x = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15))

psize_large

ggsave("figure_large_diatom_size_lat.png", psize_large, device = "png", width = 3.3, height = 6, dpi = 600)

##################################################
# Figure 4: Pro alpha sensitivity - growth rate  #
##################################################

dpaf <- read.csv("output_pro_alpha_sens_future_ws.csv")   # Pro alpha sensitivity future (warm-start)
dpap <- read.csv("output_pro_alpha_sens_present_ws.csv")  # Pro alpha sensitivity present (warm-start)

dpaf$latitude <- rep(cesm$latitude, length(unique(dpaf$ctag_max)))
dpap$latitude <- rep(cesm$latitude, length(unique(dpap$ctag_max)))

# Filter bad solutions
dpaf <- dpaf %>% filter(NC < 0.3, NC > 0.05, size > 0.6, size <= 6, mu > 0.001, mu < 0.60, temp >= 10)
dpap <- dpap %>% filter(NC < 0.3, NC > 0.05, size > 0.6, size <= 4, mu > 0.04,  mu < 0.60, temp >= 10)

dpaf$model <- "future"
dpap$model <- "present"
dpas <- rbind(dpaf, dpap)
dpas$model <- factor(dpas$model, levels = c("present", "future"))

# Build complete grid and interpolate across latitudes
dpas_grid <- expand.grid(
  latitude = latitude_seq,
  model    = levels(dpas$model),
  ctag_max = unique(dpaf$ctag_max)
)

dpas_approx <- dpas_grid %>%
  left_join(dpas %>% select(latitude, model, ctag_max, mu), by = c("latitude", "model", "ctag_max")) %>%
  group_by(model, ctag_max) %>%
  arrange(latitude) %>%
  mutate(mu = zoo::na.approx(mu, na.rm = FALSE, maxgap = 18)) %>%
  ungroup()

alpha_vals   <- sort(unique(dpas$ctag_max))
alpha_colors <- setNames(
  colorRampPalette(c("#FFF3B0", "#E09200", "#C1440E"))(length(alpha_vals)),
  as.character(alpha_vals)
)
alpha_labels <- setNames(
  c("2e6", "4e6", "6e6", "8e6", "1e7", "1.2e7", "1.4e7", "1.6e7", "1.8e7", "2e7"),
  as.character(alpha_vals)
)

p_pro_alpha <- ggplot(dpas_approx, aes(y = mu, x = latitude, col = factor(ctag_max), group = factor(ctag_max))) +
  geom_line(size = 0.8, na.rm = TRUE) +
  facet_wrap(~ model) +
  scale_color_manual(values = alpha_colors, labels = alpha_labels,
                     name = expression(atop(italic("Prochlorococcus"), alpha))) +
  xlab("Latitude") +
  ylab(expression(paste("Growth rate (", day^-1, ")", sep = ""))) +
  scale_y_continuous(breaks = c(0.1, 0.3, 0.5), labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.6), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-70, -35, 0, 35, 70), labels = scales::number_format(accuracy = 1), limits = c(-70, 70), expand = c(0, 0)) +
  mytheme +
  theme(legend.position = "right",
        legend.title = element_text(size = 13, colour = "black"),
        legend.text  = element_text(size = 11, colour = "black")) +
  coord_flip(clip = "off")

p_pro_alpha

######################################################
# Figure 4 overlay: add Synechococcus + warm diatom  #
######################################################

dsf <- read.csv("output_syn_latitudinal_future.csv")           # Synechococcus future
dsp <- read.csv("output_syn_latitudinal_present.csv")          # Synechococcus present
dwf <- read.csv("output_warm_diatom_latitudinal_future.csv")   # Warm-adapted diatom future
dwp <- read.csv("output_warm_diatom_latitudinal_present.csv")  # Warm-adapted diatom present

dsf$latitude <- cesm$latitude; dsf$model <- "future"
dsp$latitude <- cesm$latitude; dsp$model <- "present"
dwf$latitude <- cesm$latitude; dwf$model <- "future"
dwp$latitude <- cesm$latitude; dwp$model <- "present"

dsyn <- rbind(dsf, dsp) %>%
  filter(NC < 0.3, NC > 0.05, size > 0.6, size < 28, mu > 0.04, mu < 0.75) %>%
  select(latitude, model, mu)
dsyn$model <- factor(dsyn$model, levels = c("present", "future"))

dwd <- rbind(dwf, dwp) %>%
  filter(NC < 0.3, NC > 0.05, size > 0.6, size < 28, mu > 0.04, mu < 0.75) %>%
  filter(!(model == "present" & latitude == -39.5),
         !(model == "present" & latitude ==   5.5)) %>%
  select(latitude, model, mu)
dwd$model <- factor(dwd$model, levels = c("present", "future"))

overlay_grid <- expand.grid(latitude = latitude_seq, model = c("present", "future"))
overlay_grid$model <- factor(overlay_grid$model, levels = c("present", "future"))

dsyn_approx <- overlay_grid %>%
  left_join(dsyn, by = c("latitude", "model")) %>%
  group_by(model) %>% arrange(latitude) %>%
  mutate(mu = zoo::na.approx(mu, na.rm = FALSE, maxgap = 15)) %>%
  ungroup()

dwd_approx <- overlay_grid %>%
  left_join(dwd, by = c("latitude", "model")) %>%
  group_by(model) %>% arrange(latitude) %>%
  mutate(mu = zoo::na.approx(mu, na.rm = FALSE, maxgap = 15)) %>%
  ungroup()

p_pro_alpha_overlay <- p_pro_alpha +
  scale_y_continuous(breaks = c(0.1, 0.3, 0.5), labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.7), expand = c(0, 0)) +
  geom_line(data = dsyn_approx %>% mutate(species = "Synechococcus"),
            aes(x = latitude, y = mu, group = interaction(model, species), linetype = species),
            color = "#00441BFF", size = 1.2, na.rm = TRUE, inherit.aes = FALSE) +
  geom_line(data = dwd_approx %>% mutate(species = "Warm-adapted diatom"),
            aes(x = latitude, y = mu, group = interaction(model, species), linetype = species),
            color = "#762A83", size = 1.2, na.rm = TRUE, inherit.aes = FALSE) +
  scale_linetype_manual(
    values = c("Synechococcus" = "solid", "Warm-adapted diatom" = "solid"),
    labels = c("Synechococcus"       = expression(italic("Synechococcus") ~ alpha == "1e7"),
               "Warm-adapted diatom" = expression("Warm-adapted diatom" ~ alpha == "2e7")),
    name  = NULL,
    guide = guide_legend(order = 2, override.aes = list(color = c("#00441BFF", "#762A83"), linewidth = 1.2))
  ) +
  guides(colour = guide_legend(order = 1))

p_pro_alpha_overlay

ggsave("figure_pro_alpha_sens_lat_overlay.png", p_pro_alpha_overlay, device = "png", width = 7, height = 5, dpi = 600)

#############################################################
# Figure 5: Comparison Pro vs Syn / Pro vs warm diatom (AB) #
#############################################################

# Interpolation helper
interp_series <- function(df, type_label, maxgap = 15) {
  data.frame(latitude = latitude_seq) %>%
    left_join(df, by = "latitude") %>%
    arrange(latitude) %>%
    mutate(mu   = zoo::na.approx(mu,   na.rm = FALSE, maxgap = maxgap),
           size = zoo::na.approx(size, na.rm = FALSE, maxgap = maxgap),
           type = type_label)
}

# Filter to single alpha values from already-loaded dpaf (future only)
pro_1e7 <- dpaf %>%
  filter(dplyr::near(ctag_max, 1e7)) %>%
  select(latitude, mu, size)

pro_2e7 <- dpaf %>%
  filter(dplyr::near(ctag_max, 2e7)) %>%
  select(latitude, mu, size)

# Filter Syn and warm diatom (future only, reuse already-loaded dsf / dwf)
syn_f <- dsf %>%
  filter(NC < 0.3, NC > 0.05, size > 0.6, size < 28, mu > 0.04, mu < 0.75) %>%
  select(latitude, mu, size)

wd_f <- dwf %>%
  filter(NC < 0.3, NC > 0.05, size > 0.6, size < 28, mu > 0.04, mu < 0.75) %>%
  select(latitude, mu, size)

# maxgap=10 for Pro 1e7: keeps equatorial no-growth gap (13 NAs) but bridges noise gaps
# maxgap=18 for Pro 2e7: bridges 35N bad-solution cluster; no equatorial gap for this alpha
pro1_int <- interp_series(pro_1e7, "Pro (alpha=1e7)", maxgap = 10)
syn_int  <- interp_series(syn_f,   "Synechococcus",   maxgap = 15)
pro2_int <- interp_series(pro_2e7, "Pro (alpha=2e7)", maxgap = 18)
wd_int   <- interp_series(wd_f,    "Warm-adapted diatom", maxgap = 15)

# Panel A: Pro alpha=1e7 vs Synechococcus ####
dA <- bind_rows(pro1_int %>% mutate(series = "Pro (alpha=1e7)"),
                syn_int  %>% mutate(series = "Synechococcus"))
dA$series <- factor(dA$series, levels = c("Pro (alpha=1e7)", "Synechococcus"))

phen_labels_A <- c("Pro (alpha=1e7)" = expression(italic("Prochlorococcus") ~ (alpha == "1e7")),
                   "Synechococcus"   = expression(italic("Synechococcus")))

pA <- ggplot(dA, aes(x = latitude, y = mu, col = series, linetype = series)) +
  geom_line(size = 1.2, na.rm = TRUE) +
  scale_color_manual(values = c("Pro (alpha=1e7)" = "#5AAE61", "Synechococcus" = "#00441BFF"),
                     labels = phen_labels_A) +
  scale_linetype_manual(values = c("Pro (alpha=1e7)" = "solid", "Synechococcus" = "solid"),
                        labels = phen_labels_A) +
  xlab("Latitude") +
  ylab(expression(paste("Growth rate (", day^-1, ")", sep = ""))) +
  scale_y_continuous(breaks = c(0.1, 0.3, 0.5), labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.65), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-70, -35, 0, 35, 70), labels = scales::number_format(accuracy = 1), limits = c(-70, 70), expand = c(0, 0)) +
  ggtitle("A") +
  mytheme +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0, face = "bold", size = 16)) +
  coord_flip(clip = "off")

# Panel B: Pro alpha=2e7 vs Warm-adapted diatom ####
dB <- bind_rows(pro2_int %>% mutate(series = "Pro (alpha=2e7)"),
                wd_int   %>% mutate(series = "Warm-adapted diatom"))
dB$series <- factor(dB$series, levels = c("Pro (alpha=2e7)", "Warm-adapted diatom"))

phen_labels_B <- c("Pro (alpha=2e7)"     = expression(italic("Prochlorococcus") ~ (alpha == "2e7")),
                   "Warm-adapted diatom" = "Warm-adapted diatom")

pB <- ggplot(dB, aes(x = latitude, y = mu, col = series, linetype = series)) +
  geom_line(size = 1.2, na.rm = TRUE) +
  scale_color_manual(values = c("Pro (alpha=2e7)" = "#5AAE61", "Warm-adapted diatom" = "#762A83"),
                     labels = phen_labels_B) +
  scale_linetype_manual(values = c("Pro (alpha=2e7)" = "solid", "Warm-adapted diatom" = "solid"),
                        labels = phen_labels_B) +
  xlab("Latitude") +
  ylab(expression(paste("Growth rate (", day^-1, ")", sep = ""))) +
  scale_y_continuous(breaks = c(0.1, 0.3, 0.5), labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.65), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-70, -35, 0, 35, 70), labels = scales::number_format(accuracy = 1), limits = c(-70, 70), expand = c(0, 0)) +
  ggtitle("B") +
  mytheme +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0, face = "bold", size = 16)) +
  coord_flip(clip = "off")

p_combined <- pA + pB + plot_layout(ncol = 2)

ggsave("figure_pro_alpha_comparison_AB.png", p_combined, device = "png", width = 6.5, height = 5, dpi = 600)

