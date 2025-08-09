#####################################################################
# Figures: thermal growth curves for the different functional types #
#####################################################################

# cleaning the workspace
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(grid)

# Plotting settings ####
mytheme <- theme(legend.position = "top",
                 plot.margin = margin(12, 5, 10, 15),
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.background = element_rect(colour = "black", fill = "white", 
                                                 size = 1.0),
                 axis.text = element_text(size = 15, colour ="black"),
                 axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
                 axis.text.y = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),                   
                 axis.title = element_text(colour ="black"),
                 axis.title.y = element_text(size = 17, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                 axis.title.x = element_text(size = 17, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                 axis.ticks.length = unit(-0.25, "cm"),
                 legend.title = element_blank(),
                 legend.direction = "vertical",
                 legend.text = element_text(size=15, colour ="black"),
                 legend.key = element_blank(),
                 plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                 strip.background = element_rect(color = "black", fill = "white"),
                 strip.text = element_text(size = 15, color = "black"))

#######################################################################
# Thermal curves for N-replete and N-limited for different phenotypes #
#######################################################################

dcy <- read.csv("output_pro_thermal_curve.csv")
dcy$phenotype <- rep("Prochlorococcus", nrow(dcy))
dcy2 <- read.csv("output_syn_thermal_curve.csv")
dcy2$phenotype <- rep("Synechococcus", nrow(dcy2))
dw <- read.csv("output_warm_diatom_thermal_curve.csv")
dw$phenotype <- rep("Warm-adapted diatom", nrow(dw))
dc <- read.csv("output_cold_diatom_thermal_curve.csv")
dc$phenotype <- rep("Cold-adapted diatom", nrow(dc))

df <- rbind(dc,dw, dcy, dcy2)

df$DIN2 <- factor(df$DIN, levels = c("20", "0.1"), labels = c("N-replete","N-deplete"))

# Filtering out non-optimal solutions
df2 <- df %>% 
  filter(case_when(
    DIN2 == "N-deplete" ~ NC < 0.3 & NC > 0.04 & size > 0.4 & size < 50 & mu > 1e-3 & mu < 2.0,
    DIN2 == "N-replete" ~ NC < 0.3 & NC > 0.04 & size > 0.4 & size < 40 & mu > 7e-2 & mu < 2.0
  ))

df2 <- df2 %>%
  filter(!(phenotype == "Cold-adapted diatom" & temp > 23),
         !(temp > 33.3),
         !(phenotype == "Warm-adapted diatom" & DIN2 == "N-replete" & temp == 27.25),
         !(phenotype == "Prochlorococcus" & DIN2 == "N-replete" & temp %in% c(31)),
         !(phenotype == "Prochlorococcus" & DIN2 == "N-deplete" & temp %in% c(20.5, 25.5, 29, 30.5)),
         !(phenotype == "Synechococcus" & DIN2 == "N-replete" & temp %in% c(10.5)),
         !(phenotype == "Synechococcus" & DIN2 == "N-deplete" & temp %in% c(12.25, 13.5, 21)))

# estimate nitrogen uptake affinity
df2$radius <- df2$size/2
df2$cellvol <- 4*pi*(df2$radius^3)/3
df2$n <- df2$ctr * df2$cellvol
s_tr <- 1.26e-5      # specific surface area of membrane transporter molecules (um2)
s <- sqrt(s_tr/4*pi) # the radius of a membrane transporter (um)
D <- 1.7e-9          # diffusion rate of a nutrient molecule in seawater
df2$p <- df2$n*pi*s^2/(4*pi*df2$radius^2)
df2$aff <- 4*pi*D*df2$radius*(df2$n*s)/(df2$n*s + pi*df2$radius*(1-df2$p))
max(df2$aff)
min(df2$aff)
df2$saff <- df2$aff/df2$cellvol
df2$istress <- df2$itag * df2$ilpb

d1 <- df2 %>% group_by(DIN2, phenotype) %>% filter(mu == max(mu, na.rm = TRUE))

# Let's make sure genera appear in italic in our plots
phen_labels <- c(
  "Prochlorococcus" = expression(italic("Prochlorococcus")),
  "Synechococcus"   = expression(italic("Synechococcus")),
  "Warm-adapted diatom" = "Warm-adapted diatom",
  "Cold-adapted diatom" = "Cold-adapted diatom")

p1 <- ggplot(data = df2 %>% filter(phenotype == "Warm-adapted diatom"), aes(x = temp, y = mu, linetype = DIN2)) +
  geom_line(size = 1) +
  facet_wrap(~ phenotype) +
  scale_y_continuous(breaks = c(0.0, 0.7, 1.4), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.4)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(7,33)) +
  scale_linewidth_manual(values = c(1, 2)) +
  ylab("") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = "", axis.text.y = element_blank()) +
  geom_point(data = d1 %>% filter(phenotype == "Warm-adapted diatom"), aes(x = temp, y = mu), size = 5, color = "#762A83") +
  geom_vline(data = d1 %>% filter(phenotype == "Warm-adapted diatom"), aes(xintercept = temp, linetype = DIN2), color = "#762A83", size = 1.0, show.legend = FALSE)
p1

p2 <- ggplot(data = df2 %>% filter(phenotype == "Cold-adapted diatom"), aes(x = temp, y = mu, linetype = DIN2)) +
  geom_line(size = 1) +
  facet_wrap(~ phenotype) +
  scale_y_continuous(breaks = c(0.0, 0.7, 1.4), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.4)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(7,33)) +
  scale_linewidth_manual(values = c(1, 2)) +
  ylab("") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = "", axis.text.y = element_blank()) +
  geom_point(data = d1 %>% filter(phenotype == "Cold-adapted diatom"), aes(x = temp, y = mu), size = 5, color = "#C2A5CF") +
  geom_vline(data = d1 %>% filter(phenotype == "Cold-adapted diatom"), aes(xintercept = temp, linetype = DIN2), color = "#C2A5CF", size = 1.0, show.legend = FALSE)
p2

p3 <- ggplot(data = df2 %>% filter(phenotype == "Prochlorococcus"), aes(x = temp, y = mu, linetype = DIN2)) +
  geom_line(size = 1) +
  facet_wrap(~ phenotype) +
  scale_y_continuous(breaks = c(0.0, 0.7, 1.4), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.4)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(7,33)) +
  scale_linewidth_manual(values = c(1, 2)) +
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = "",
        strip.text     = element_text(face = "italic")) +
  geom_point(data = d1 %>% filter(phenotype == "Prochlorococcus"), aes(x = temp, y = mu), size = 5, color = "#5AAE61") +
  geom_vline(data = d1 %>% filter(phenotype == "Prochlorococcus"), aes(xintercept = temp, linetype = DIN2), color = "#5AAE61", size = 1.0, show.legend = FALSE)
p3

p4 <- ggplot(data = df2 %>% filter(phenotype == "Synechococcus"), aes(x = temp, y = mu, linetype = DIN2)) +
  geom_line(size = 1) +
  facet_wrap(~ phenotype) +
  scale_y_continuous(breaks = c(0.0, 0.7, 1.4), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.4)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(7,33)) +
  scale_linewidth_manual(values = c(1, 2)) +
  ylab("") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  mytheme +
  theme(legend.position = "", axis.text.y = element_blank(), strip.text     = element_text(face = "italic")) +
  geom_point(data = d1 %>% filter(phenotype == "Synechococcus"), aes(x = temp, y = mu), size = 5, color = "#00441BFF") +
  geom_vline(data = d1 %>% filter(phenotype == "Synechococcus"), aes(xintercept = temp, linetype = DIN2), color = "#00441BFF", size = 1.0, show.legend = FALSE)
p4

p1 + p2 + p3 + p4

ggsave("figureS1.png", (p3 + p4 + p1 + p2) + plot_layout(nrow = 1), device = "png", width = 12, height = 3.5, dpi = 600) 

######################################################
# Model validation against Anderson et al. 2021 data #
######################################################

# Import data from Anderson et al.
da <- read.csv("data_anderson.csv")

da2 <- da %>% group_by(Group) %>% summarize(min_topt = min(t_opt))  

da <- da %>% mutate(adapt = ifelse(t_opt < 20.6, "cold", "warm")) %>%
  filter(Group %in% c("cyanobacteria", "diatoms")) %>%
  mutate(type = ifelse(Group == "diatoms" & adapt == "warm", "DW",
                       ifelse(Group == "diatoms" & adapt == "cold", "DC",
                              ifelse(Group == "cyanobacteria" & adapt == "warm", "C", "other"))))
da2 <- da %>%
  mutate(
    Group2 = case_when(str_starts(Strain, "Prochlorococcus") ~ "Pro",
                       str_starts(Strain, "Synechococcus")  ~ "Syn", 
                       TRUE ~ Group),
    type2 = case_when(str_starts(Strain, "Prochlorococcus") ~ "Pro",
                       str_starts(Strain, "Synechococcus")  ~ "Syn", 
                       TRUE ~ type),
  )

da %>% count(Group)
da %>% count(type)
da2 %>% count(Group2)
da2 %>% count(type2)

da2$type2 <- factor(da2$type2, levels = c("Pro", "Syn", "DW", "DC"), labels = c("Pro", "Syn", "W-diatom", "C-diatom"))

da2median <- da2 %>% group_by(type2) %>% summarize(median_umax = median(max_GR),
                                                   median_slope = median(left_slope))
da2median

# now retrieve summary data from model simulations for validation

df2_mu <- df2 %>%
  group_by(phenotype, DIN2) %>%
  filter(mu == max(mu, na.rm = TRUE)) %>% 
  select(mu, phenotype, temp, size)

df2_mu

df2_replete <- df2 %>% filter(DIN2 == "N-replete")

data_slope <- data.frame(type = df2_replete$phenotype, 
                         Temp = df2_replete$temp,
                         GR = df2_replete$mu) %>%
  mutate(type = factor(type, levels = c("Cold-adapted diatom", "Warm-adapted diatom", "Prochlorococcus", "Synechococcus")), Curve = as.integer(factor(type))) %>%
  relocate(Curve, Temp, GR, .before = type)

# sub-setting model thermal curve data to fit the data to the Norberg equation (just like we did for the Anderson et al. dataset)
data_slope2 <- data_slope %>%
  group_by(type) %>%
  arrange(Temp) %>%
  mutate(
    temp_maxGR = Temp[which.max(GR)],
    GR_max     = max(GR),
    GR_min_low  = if (any(Temp < temp_maxGR)) min(GR[Temp < temp_maxGR]) else NA_real_,
    GR_min_high = if (any(Temp > temp_maxGR)) min(GR[Temp > temp_maxGR]) else NA_real_,
    mid_low  = (GR_max + GR_min_low) / 2,
    mid_high = (GR_max + GR_min_high) / 2
  ) %>%
  filter(
    GR == GR_max |
      GR == GR_min_low |
      GR == GR_min_high |
      (Temp < temp_maxGR & !is.na(mid_low) &
         abs(GR - mid_low) == min(abs(GR[Temp < temp_maxGR] - mid_low))) |
      (Temp > temp_maxGR & !is.na(mid_high) &
         abs(GR - mid_high) == min(abs(GR[Temp > temp_maxGR] - mid_high)))
  ) %>%
  ungroup()

write.csv(data_slope2, "sample_data.csv", row.names = FALSE)

source("estimating_growth_slope_from_Anderson_dataset.R")

dsim <- data.frame(type = c("Pro","Syn", "W-diatom", "C-diatom"),
                   umax = c(0.532, 0.738, 1.26, 0.550),
                   left_slope = c(0.056, 0.079, 0.13, 0.061),     
                   delta_topt = c(0.2, 1.3, 2.3, -3.5),
                   size = c(0.647, 0.685, 5.13, 5.02),
                   sizeNrep = c(0.786, 1.37, 9.89, 8.92))

dsim$type <- factor(dsim$type, levels = c("Pro","Syn", "W-diatom", "C-diatom"))


p1 <- ggplot(da2, aes(x = type2, y = left_slope, fill = type2)) +
  geom_boxplot(width = 0.5) +
  ylab("Growth slope") +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF", "#762A83", "#C2A5CF")) +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2), labels = scales::number_format(accuracy = 0.1), limits = c(0,0.2)) +
  xlab("") +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  geom_point(data = dsim, aes(x = type, y = left_slope, fill = type),
             shape = 24, size = 4, color = "black", stroke = 1.2, 
             position = position_nudge(x = 0.4)) +
  mytheme +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none",
        plot.margin = margin(12, 5, 1, 15))
p1


p2 <- ggplot(da2, aes(x = type2, y = max_GR, fill = type2)) +
  geom_boxplot(width = 0.5, show.legend = FALSE) +
  ylab(expression(paste(,italic(mu[max]) ," (", day^-1,")", sep=""))) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF"), guide = "none") +
  scale_y_continuous(breaks = c(0.0, 1.5, 3.0), labels = scales::number_format(accuracy = 0.1), limits = c(0,3)) +
  xlab("") +
  geom_jitter(aes(shape = "Anderson et al. 2021"), width = 0.1, size = 1, alpha = 0.3, show.legend = TRUE) +
  geom_point(data = dsim, aes(x = type, y = umax, shape = "Model", fill = type),
             size = 4, color = "black", stroke = 1.2, 
             position = position_nudge(x = 0.4), show.legend = TRUE) +
  scale_shape_manual(values = c("Anderson et al. 2021" = 21, "Model" = 24)) +
  mytheme +
  ggtitle("Model validation")+
  theme(
    legend.position = c(0.3,0.85),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    plot.title = element_text(face = "plain", size = 18)) +
  coord_cartesian(clip = "off") +
  annotation_custom(grob = textGrob("A", x = unit(-0.25, "npc"), y = unit(1.15, "npc"),
                    just = "left", gp = gpar(col = "black", fontsize = 26))) +
  guides(shape = guide_legend(nrow = 3, override.aes = list(fill = "black")))
p2

pcol1 <- p2/p1
pcol1

p3 <- ggplot(dsim, aes(x = type, y = delta_topt, fill = type)) +
  geom_point(shape = 24, size = 5, color = "black", stroke = 1.2) +
  ylab(expression(paste(,Delta * italic(T[opt]) ," (", degree,"C)", sep=""))) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  scale_y_continuous(breaks = c(-2.5, 0, 2.5), labels = scales::number_format(accuracy = 0.1), limits = c(-3.6,2.5)) +
  xlab("")+
  geom_hline(yintercept = 0, linetype = "dotted") +
  mytheme +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none") 
p3

p4 <- ggplot(dsim, aes(x = type, y = size, fill = type)) +
  geom_point(shape = 24, size = 5, color = "black", stroke = 1.2) +
  ylab(expression(paste("Cell length (", mu,"m)", sep=""))) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  scale_y_continuous(breaks = c(0.0, 3.0, 6.0), labels = scales::number_format(accuracy = 0.1), limits = c(0,6)) +
  xlab("")+
  ggtitle("Emergent traits")+
  mytheme +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none", axis.text.x = element_blank(),
        plot.margin = margin(10, 5, 10, 15), plot.title = element_text(face = "plain", size = 18)) +
  annotation_custom(grob = textGrob("B", x = unit(-0.25, "npc"), y = unit(1.15, "npc"),
                                    just = "left", gp = gpar(col = "black", fontsize = 26)))
p4

pcol2b <- p4 / p3

pcol1 | pcol2b

p4s <- ggplot(dsim, aes(x = type, y = sizeNrep, fill = type)) +
  geom_point(shape = 24, size = 4, color = "black", stroke = 1.2) +
  ylab(expression(paste("Cell length (", mu,"m)", sep=""))) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  scale_y_continuous(breaks = c(0.0, 7.0, 14.0), labels = scales::number_format(accuracy = 0.1), limits = c(0,14)) +
  xlab("")+
  mytheme +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none",
        plot.margin = margin(15, 5, 10, 15))

p4s

ggsave("figureS2.png", p4s, device = "png", width = 5, height = 4, dpi = 600) 

# Nitrogen uptake
dvtr <- df2 %>%
  group_by(phenotype, DIN2) %>%
  mutate(vtr_scaled = vtr / max(vtr, na.rm = TRUE)) %>%
  ungroup()

dvtr$phenotype <- factor(dvtr$phenotype, levels = c("Prochlorococcus", "Synechococcus", "Warm-adapted diatom", "Cold-adapted diatom"))

pnup <- ggplot(dvtr %>% filter(DIN2 == "N-deplete"), aes(x = temp, y = vtr_scaled, col = phenotype)) +
  geom_line(size = 2) +
  ylab(expression(paste("Scaled ", italic(N[upt]), sep = ""))) +
  xlab("") +
  scale_y_continuous(breaks = c(0.4, 0.7, 1.0), labels = scales::number_format(accuracy = 0.1), limits = c(0.3,1.1)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,33)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF"), labels = phen_labels) +
  mytheme +
  ggtitle("Emergent traits")+
  theme(strip.text = element_blank(), legend.position = "",
        axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(), plot.margin = margin(30, 5, 10, 15),
        plot.title = element_text(face = "plain", size = 18)) +
  coord_cartesian(clip = "off") +
  annotation_custom(grob = textGrob("C", x = unit(-0.25, "npc"), y = unit(1.15, "npc"),
                                    just = "left", gp = gpar(col = "black", fontsize = 26)))
pnup

df2$phenotype <- factor(df2$phenotype, levels = c("Prochlorococcus", "Synechococcus", "Warm-adapted diatom", "Cold-adapted diatom"))

# Carbon use efficiency
pcue <- ggplot(df2 %>% filter(DIN2 == "N-deplete"), aes(x = temp, y = 1 - (resp/phot), col = phenotype)) +
  geom_line(size = 2) +
  ylab("CUE") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,33)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF"), labels = phen_labels) +
  mytheme +
  coord_cartesian(clip = "off") +
  theme(strip.text = element_blank(), legend.position = c(0.35,0.3), legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_text(size = 16, colour ="black"))
pcue

pcol3 <- pnup / pcue
pcol3

pcol1 | pcol2b | pcol3

ggsave("figure3.png", pcol1 | pcol2b | pcol3, device = "png", width = 15, height = 8.5, dpi = 600) 

###############################################################
# Figure S2: proteome investments N-replete versus N-limiting #
###############################################################

p1 <- ggplot(df2, aes(x = temp, y = itr, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("transporters") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(n.breaks = 4) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF"), labels = phen_labels) +
  mytheme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(strip.background = element_blank(), legend.position = "top",
        axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.justification = c(0, 1))
p1

p2 <- ggplot(df2, aes(x = temp, y = icbd, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("glycolysis") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2), labels = scales::number_format(accuracy = 0.1), limits = c(0,.2)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.text = element_blank(), strip.background = element_blank(), 
        legend.position = "", axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
p2

p3 <- ggplot(df2, aes(x = temp, y = istress, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("heat-mitigation") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(n.breaks = 3) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.text = element_blank(), strip.background = element_blank(), 
        legend.position = "", axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
p3

p4 <- ggplot(df2, aes(x = temp, y = size, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab(expression(paste("cell length (", mu,"m)", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_log10()+
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.text = element_blank(), legend.position = "",
        axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
p4

p5 <- ggplot(df2, aes(x = temp, y = ichl, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("photosystems") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(breaks = c(0, 0.15, 0.3), labels = scales::number_format(accuracy = 0.01), limits = c(0,0.3)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.text = element_blank(), strip.background = element_blank(), 
        legend.position = "", axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
p5

p6 <- ggplot(df2, aes(x = temp, y = irub, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("rubisco") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), labels = scales::number_format(accuracy = 0.01), limits = c(0,0.2)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.text = element_blank(), strip.background = element_blank(), 
        legend.position = "", axis.title.x = element_text(size = 17, colour ="black"),
        axis.text.x = element_text(size = 16, colour ="black"))
p6

p7 <- ggplot(df2, aes(x = temp, y = irib, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("ribosomes") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(n.breaks = 4) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.background = element_blank(), 
        legend.position = "", axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
p7

p8 <- ggplot(df2, aes(x = temp, y = ilpb, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("lipid synthesis") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), labels = scales::number_format(accuracy = 0.01), limits = c(0,0.2)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.text = element_blank(), strip.background = element_blank(), 
        legend.position = "", axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
p8

p9 <- ggplot(df2, aes(x = temp, y = icha, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("repair") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(breaks = c(0, 0.04, 0.08), labels = scales::number_format(accuracy = 0.01), limits = c(0,0.08)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.text = element_blank(), strip.background = element_blank(), 
        legend.position = "", axis.title = element_text(size = 18, colour ="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
p9

p10 <- ggplot(df2, aes(x = temp, y = ctag, col = phenotype)) +
  facet_wrap(~DIN2) +
  geom_line(size = 2) +
  ylab("antioxidants") +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  scale_y_continuous(n.breaks = 4) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(8,31)) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF", "#762A83","#C2A5CF")) +
  mytheme +
  theme(strip.text = element_blank(), strip.background = element_blank(), 
        legend.position = "", axis.title.x = element_text(size = 17, colour ="black"),
        axis.text.x = element_text(size = 16, colour ="black"))
p10

pall <- (p1 + p7) / (p2 + p8) / (p3 + p9) / (p5 + p4) / (p6 + p10)

ggsave("figureS3.png", pall, device = "png", width = 11, height = 11, dpi = 600) 

