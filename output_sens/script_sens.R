############################################################
# Figure 2: growth trade-offs and shifts in thermal traits #
############################################################

# cleaning the workspace
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(RColorBrewer)

# Plotting settings
mytheme <-   theme(panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_rect(colour = "black", fill = "white", 
                                                   size = 1.5),
                   axis.text = element_text(size = 15, colour ="black"),
                   axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
                   axis.text.y = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),                   
                   axis.title = element_text(size = 17, colour ="black"),
                   axis.title.y = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                   axis.title.x = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                   axis.ticks.length = unit(0.25, "cm"),
                   plot.margin = unit(c(1, 1.5, 0.5, 0.5), "lines"),
                   legend.position = "top",
                   legend.title = element_text(size=15, colour ="black"),
                   legend.text = element_text(size=12, colour ="black"),
                   legend.key = element_rect(fill = ""),
                   plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                   strip.background = element_rect(color = "black", fill = "white"),
                   strip.text = element_text(size = 15, color = "black")
)

#########################################################################
# Individual effects of thermal dependency and heat-mitigation capacity #
#########################################################################

# Model output for variable Ea
dea <- read.csv("output_sensitivity_analysis_Ea.csv", header = T)

# First filter pass to remove local minima
dea <- dea %>% filter(NC < 0.5 & NC > 0.008 & size > 8 & size < 40 & mu > 0.1)

dea$panel <- rep("thermal dependency (Ea)", nrow(dea))

# Plot with facet label in italic
p_Ea <- ggplot(dea, aes(x = temp, y = mu, col = Ea)) +
  geom_line(aes(group = factor(Ea)), size = 1) +
  facet_wrap(~ panel, labeller = label_value) +  
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  xlim(6, 31) +
  scale_y_continuous(breaks = c(0.0, 0.8, 1.6), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.6)) +
  scale_color_gradient(low = "gray", high = "black",
                       name = NULL,
                       breaks = seq(min(dea$Ea), max(dea$Ea), length.out = 2),
                       labels = c("low", "high"),
                       guide = guide_colorbar(barwidth = 1, barheight = 3.5)) +
  mytheme +
  theme(legend.position = c(0.18, 0.78))
p_Ea

# Model output for variable alpha
da <- read.csv("output_sensitivity_analysis_alpha.csv", header = T)

da <- da %>% filter(ctag_max %in% c(4e6, 6e6, 9600000, 13200000, 16800000, 20400000))

da$panel <- rep("heat-mitigation (α)", nrow(da))

# First filter pass to remove local minima
da <- da %>% filter(NC < 0.5 & NC > 0.008 & size > 8 & size < 200 & mu > 0.01)

# Plot with facet label in italic
p_alpha <- ggplot(da, aes(x = temp, y = mu, col = ctag_max)) +
  geom_line(aes(group = factor(ctag_max)), size = 1) +
  facet_wrap(~ panel, labeller = label_value) +  
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  xlab(expression(paste("Temperature  (", degree,"C)"))) +
  xlim(6, 33) +
  scale_y_continuous(breaks = c(0.0, 0.8, 1.6), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.6)) +
  scale_color_gradient(low = "gray", high = "black",
                       name = NULL,
                       breaks = seq(min(da$ctag_max), max(da$ctag_max), length.out = 2),
                       labels = c("low", "high"),
                       guide = guide_colorbar(barwidth = 1, barheight = 3.5)) +
  mytheme +
  theme(legend.position = c(0.18, 0.78),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
p_alpha

prow1 <- p_Ea + p_alpha
prow1

#################
# Topt analyses #
#################

df <- read.csv("output_heatmaps.csv", header = T)    
df$DIN2 <- factor(df$DIN, levels = c("20", "0.1"), labels = c("high","low"))

# First filter pass to remove local minima
df2 <- df %>% 
  filter(case_when(
    DIN2 == "low" ~ NC < 0.3 & NC > 0.05 & size > 0.0009 & size < 6 & mu > 3e-3,
    DIN2 == "high" ~ NC < 0.5 & NC > 0.008 & size > 8 & size < 40 & mu > 0.1
  ))

# Checking for more local minima
ggplot(df2 %>% filter(DIN2 == "high"), aes(x = temp, y = mu)) +
  facet_grid(ctag_max ~ Ea) +
  geom_line()

ggplot(df2 %>% filter(DIN2 == "low"), aes(x = temp, y = mu)) +
  facet_grid(ctag_max ~ Ea) +
  geom_line()

ggplot(df2 %>% filter(DIN2 == "high"), aes(x = temp, y = vtr)) +
  facet_grid(ctag_max ~ Ea) +
  geom_line()

ggplot(df2 %>% filter(DIN2 == "low"), aes(x = temp, y = vtr)) +
  facet_grid(ctag_max ~ Ea) +
  geom_line()

# removing any leftover local minima
df2 <- df2 %>% filter(ctag_max > 4e6,
                      Ea > 0.5,
                      !(DIN2 == "high" & Ea == 0.5 & ctag_max == 4e6 & temp == 26.5),
                      !(DIN2 == "low" & Ea == 0.55 & ctag_max == 8e6 & temp == 29),
                      !(DIN2 == "low" & Ea == 0.72 & ctag_max == 8e6 & temp == 20.75),
                      !(DIN2 == "low" & Ea == 0.88 & ctag_max == 1e7 & temp == 25.75),
                      !(DIN2 == "low" & Ea == 0.77 & ctag_max == 6e6 & temp == 22),
                      !(DIN2 == "low" & Ea == 0.5 & ctag_max == 2.2e7 & temp == 18.5))

# Now let's interpolate data across temperatures with the optimal solutions we have

# Generate a complete grid of temperature values
temp_seq <- seq(18, 29, by = 0.25)

# Create a complete grid with all combinations of temperature other factor variables
complete_grid <- expand.grid(
  temp = temp_seq,
  DIN2 = unique(df2$DIN2),  
  Ea = unique(df2$Ea),
  ctag_max = unique(df2$ctag_max))

# Join the complete grid with our original data
df3 <- complete_grid %>%
  left_join(df2, by = c("temp", "DIN2", "Ea", "ctag_max"))

# Interpolate missing values based on temperature
df3_approx <- df3 %>%
  group_by(DIN2, Ea, ctag_max) %>%
  mutate(mu = zoo::na.approx(mu, na.rm = FALSE),
         vtr = zoo::na.approx(vtr, na.rm = FALSE))

# Plot maximum growth rates under nutrient-limited conditions
d1 <- df3_approx %>% group_by(DIN2, Ea, ctag_max) %>% filter(mu == max(mu, na.rm = TRUE))

d1b <- d1 %>% filter(DIN2 == "low") 

pumax <- ggplot(d1b, aes(x = factor(Ea), y = factor(ctag_max), fill = mu)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = brewer.pal(n = 11, name = "YlOrBr"),
    limits = c(0.32,0.55),  
    breaks = c(0.32, 0.42, 0.52),
    name = expression(paste(italic(mu[max])," (", day^-1, ")", sep = "")),
    guide = guide_colorbar(barwidth = 16, barheight = 1, title.position = "top", title.hjust = 0.5, frame.colour = "black",frame.linewidth = 0.5)) +
  coord_cartesian(expand = FALSE) +
  ylab(expression(paste("heat-mitigation (" ,italic(α), ")", sep = ""))) +
  xlab(expression(paste("thermal dependency (", italic(E[a]),")", sep = ""))) +
  scale_x_discrete(breaks = c("0.55", "1.055"), labels = c("low", "high")) + 
  scale_y_discrete(breaks = c("6e+06", "2.4e+07"), labels = c("low", "high")) +  
  annotate("text", x = 2, y = 9.7, label = "TT", size = 10, color = "black", vjust = 1.5, fontface = "bold") +
  annotate("text", x = 8.8, y = 3.1, label = "NE", size = 10, color = "black", vjust = 1.5, fontface = "bold") +
  theme(
    legend.title = element_text(size = 15, angle = 0),
    legend.text = element_text(size = 13),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15, colour = "black"),
    axis.title.x = element_text(size = 17),
    legend.position = "top",       
    legend.direction = "horizontal",
    panel.border = element_rect(color = "black", fill = NA, size = 1))
pumax

# Plot delta Topt (Topt N-limited - Topt N-replete)
df_ratio <- d1 %>%
  group_by(Ea, ctag_max) %>%
  summarise(delta_temp = temp[DIN2 == "low"] - temp[DIN2 == "high"],
            delta_vtr = vtr[DIN2 == "low"] - vtr[DIN2 == "high"])

# ensure inaccurate delta_temp are removed from the dataset
# this can happen if there were multiple consecutive missing optimal solutions
df_ratio2 <- df_ratio %>%
  arrange(Ea, ctag_max) %>%
  group_by(Ea) %>%
  mutate(prev = lag(delta_temp),
         nextt = lead(delta_temp),
         delta_temp = case_when(delta_temp >= 0 & delta_temp > coalesce(prev, delta_temp) ~ coalesce(prev, delta_temp),
                                delta_temp < 0 & (delta_temp < coalesce(prev, delta_temp) & delta_temp < coalesce(nextt, delta_temp)) ~ coalesce(nextt, delta_temp),
                                TRUE ~ delta_temp)) %>%
  select(-prev, -nextt) %>%
  ungroup()

df_ratio3 <- df_ratio2 %>%
  arrange(ctag_max, Ea) %>%         
  group_by(ctag_max) %>%            
  mutate(prev = lag(delta_temp),
         nextt = lead(delta_temp),
         delta_temp = case_when(delta_temp >= 0 & delta_temp < coalesce(prev, delta_temp) ~ coalesce(prev, delta_temp),
                                delta_temp < 0 & (delta_temp > coalesce(prev, delta_temp) & delta_temp > coalesce(nextt, delta_temp)) ~ coalesce(nextt, delta_temp),
                                TRUE ~ delta_temp)) %>%
  select(-prev, -nextt) %>%
  ungroup()

pdelta <- ggplot(df_ratio3, aes(x = factor(Ea), y = factor(ctag_max), fill = delta_temp)) +
  geom_raster()  +
  scale_fill_gradientn(
    colours = rev(brewer.pal(n = 11, name = "RdBu")),
    na.value = "grey50",  
    name = expression(paste(italic(Delta * T[opt])," (", degree, "C)", sep = "")),
    guide = guide_colorbar(barwidth = 16, barheight = 1, title.position = "top", title.hjust = 0.5, frame.colour = "black",frame.linewidth = 0.5),
    limits = c(-5, 5), 
    values = scales::rescale(c(-5, 0, 5))  
  ) +
  coord_cartesian(expand = FALSE) +
  ylab(expression(paste("heat-mitigation (", italic(α), ")", sep = ""))) +
  xlab(expression(paste("thermal dependency (", italic(E[a]), ")", sep = ""))) +
  scale_x_discrete(breaks = c("0.55", "1.055"), labels = c("low", "high")) + 
  scale_y_discrete(breaks = c("6e+06", "2.4e+07"), labels = c("low", "high")) +  
  theme(
    legend.title = element_text(size = 15, angle = 0),
    legend.text = element_text(size = 13),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 17),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",       
    legend.direction = "horizontal",
    panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  annotate("text", x = 2, y = 9.7, label = "TT", size = 10, color = "black", vjust = 1.5, fontface = "bold") +
  annotate("text", x = 8.8, y = 3.1, label = "NE", size = 10, color = "black", vjust = 1.5, fontface = "bold")

pdelta

# Plot temperature at maximum nutrient uptake rate under nutrient-limited conditions
d2 <- df3_approx %>% filter(DIN2 == "low") %>% group_by(DIN2, Ea, ctag_max) %>% filter(vtr == max(vtr, na.rm = TRUE))

d2b <- d2 %>%
  arrange(ctag_max, Ea) %>%           
  group_by(ctag_max) %>%              
  mutate(prev_temp = lag(temp),
         temp = if_else(!is.na(prev_temp) & temp < prev_temp, prev_temp, temp)) %>%
  select(-prev_temp) %>%              
  ungroup()

pnut <- ggplot(d2b, aes(x = factor(Ea), y = factor(ctag_max), fill = temp)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = brewer.pal(n = 11, name = "YlOrBr"),
    limits = c(20,30), 
    breaks = c(20, 25, 30),
    name = expression(paste(italic(T)," (", degree, "C)", " at max ", italic(N[upt]), sep = "")),
    guide = guide_colorbar(barwidth = 16, barheight = 1, title.position = "top", title.hjust = 0.5, frame.colour = "black",frame.linewidth = 0.5)) +
  coord_cartesian(expand = FALSE) +
  ylab(expression(paste("heat-mitigation (" ,italic(α), ")", sep = ""))) +
  xlab(expression(paste("thermal dependency (", italic(E[a]),")", sep = ""))) +
  scale_x_discrete(breaks = c("0.5", "1"), labels = c("low", "high")) + 
  scale_y_discrete(breaks = c("6e+06", "2.4e+07"), labels = c("low", "high")) +  
  theme(
    legend.title = element_text(size = 15, angle = 0),
    legend.text = element_text(size = 13),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",       
    legend.direction = "horizontal",
    panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  annotate("text", x = 2, y = 9.7, label = "TT", size = 10, color = "black", vjust = 1.5, fontface = "bold") +
  annotate("text", x = 8.8, y = 3.1, label = "NE", size = 10, color = "black", vjust = 1.5, fontface = "bold")

pnut

#################
# Tmax analyses #
#################

dmax <- read.csv("output_heatmap_Tmax.csv", header = T)    
dmax$DIN2 <- factor(dmax$DIN, levels = c("20", "0.1"), labels = c("high","low"))

# First filter pass to remove local minima
dmax2 <- dmax %>% 
  filter(case_when(
    DIN2 == "low" ~ NC < 0.2 & NC > 0.05 & size > 0.1 & size < 10 & mu > 0.01 & mu < 1 & vtr < 4000))

# Checking for more local minima
ggplot(dmax2 %>% filter(DIN2 == "low"), aes(x = temp, y = mu)) +
  facet_grid(ctag_max ~ Ea) +
  geom_line() +
  ylim(0,0.6)

ggplot(dmax2 %>% filter(DIN2 == "low"), aes(x = temp, y = vtr)) +
  facet_grid(ctag_max ~ Ea) +
  geom_line()

# removing any leftover local minima
dmax2 <- dmax2 %>% filter(!(Ea == 0.77 & ctag_max == 6e6 & temp == 32.75),
                          !(Ea == 0.55 & ctag_max == 6e6 & temp == 31),
                          !(Ea == 0.83 & ctag_max == 8e6 & temp %in% c(33.25, 33.75)),
                          !(Ea == 0.88 & ctag_max == 8e6 & temp == 33.75),
                          !(Ea == 1.05 & ctag_max == 8e6 & temp == 34),
                          !(Ea == 0.61 & ctag_max == 1.2e7 & temp == 33.5),
                          !(Ea == 0.94 & ctag_max == 1.4e7 & temp == 35.5),
                          !(Ea == 0.94 & ctag_max == 1.4e7 & temp == 35.5),
                          !(Ea == 0.99 & ctag_max == 1.4e7 & temp == 34.5),
                          !(Ea == 1.05 & ctag_max == 1.4e7 & temp == 35.25),
                          !(Ea == 0.66 & ctag_max == 1.6e7 & temp == 34.5),
                          !(Ea == 0.88 & ctag_max == 2e7 & temp == 33.5),
                          !(Ea == 1.05 & ctag_max == 2e7 & temp == 35))

dmax3 <- dmax2 %>% group_by(ctag_max, Ea) %>% filter(temp == max(temp))

# With the plot below, we find that Tmax is approx. constant across Ea values. 
ggplot(dmax3, aes(x = factor(Ea), y = factor(ctag_max), fill = temp)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = brewer.pal(n = 11, name = "YlOrBr"),
    limits = c(32,36), 
    breaks = c(32, 34, 36),
    name = expression(paste(italic(T[max])," (", degree, "C)", sep = "")),
    guide = guide_colorbar(barwidth = 16, barheight = 1, title.position = "top", title.hjust = 0.5, frame.colour = "black",frame.linewidth = 0.5)) +
  coord_cartesian(expand = FALSE) +
  ylab(expression(paste("heat-mitigation (" ,italic(α), ")", sep = ""))) +
  xlab(expression(paste("thermal dependency (", italic(E[a]),")", sep = ""))) +
  scale_x_discrete(breaks = c("0.55", "1.05"), labels = c("low", "high")) + 
  scale_y_discrete(breaks = c("6e+06", "2.4e+07"), labels = c("low", "high")) +  
  theme(
    legend.title = element_text(size = 15, angle = 0),
    legend.text = element_text(size = 13),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 17),
    axis.ticks = element_blank(),
    legend.position = "top",       
    legend.direction = "horizontal",
    panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  annotate("text", x = 2, y = 9.7, label = "TT", size = 10, color = "black", vjust = 1.5, fontface = "bold") +
  annotate("text", x = 8, y = 3, label = "NE", size = 10, color = "black", vjust = 1.5, fontface = "bold")

# However, it is hard to determine accurate Tmax values for every thermal curve (optimizer might not find optimal solutions when growth is close to zero).
# We thus assume that 1) Tmax is equal across Ea values and 2) Tmax is identified as the mode across Ea values.
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

Tmax_mode <- dmax3 %>%
  group_by(ctag_max) %>%
  summarise(temp_mode = get_mode(temp), .groups = "drop")
Tmax_mode

dmax4 <- dmax3 %>%
  left_join(Tmax_mode, by = "ctag_max") %>%
  mutate(temp = temp_mode) %>%
  select(-temp_mode)

# Plot maximal survival temperature under nutrient-limited conditions
ptmax <- ggplot(dmax4, aes(x = factor(Ea), y = factor(ctag_max), fill = temp)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = brewer.pal(n = 11, name = "YlOrBr"),
    limits = c(31.7,34.3), 
    breaks = c(32, 33, 34),
    name = expression(paste(italic(T[max])," (", degree, "C)", sep = "")),
    guide = guide_colorbar(barwidth = 16, barheight = 1, title.position = "top", title.hjust = 0.5, frame.colour = "black",frame.linewidth = 0.5)) +
  coord_cartesian(expand = FALSE) +
  ylab(expression(paste("heat-mitigation (" ,italic(α), ")", sep = ""))) +
  xlab(expression(paste("thermal dependency (", italic(E[a]),")", sep = ""))) +
  scale_x_discrete(breaks = c("0.55", "1.05"), labels = c("low", "high")) + 
  scale_y_discrete(breaks = c("6e+06", "2.4e+07"), labels = c("low", "high")) +  
  theme(
    legend.title = element_text(size = 15, angle = 0),
    legend.text = element_text(size = 13),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 17),
    axis.ticks = element_blank(),
    legend.position = "top",       
    legend.direction = "horizontal",
    panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  annotate("text", x = 2, y = 9.7, label = "TT", size = 10, color = "black", vjust = 1.5, fontface = "bold") +
  annotate("text", x = 8.8, y = 3.1, label = "NE", size = 10, color = "black", vjust = 1.5, fontface = "bold")
ptmax

prow3 <- pdelta + pnut

prow4 <- ptmax + pumax

prow3 / prow4

##############################################################################
# Two mutually exclusive strategies: nutrient-efficient and thermal-tolerant #
##############################################################################

dtt <- read.csv("output_sens_thermal_curve_TT.csv", header = T) 
dne <- read.csv("output_sens_thermal_curve_NE.csv", header = T)    
ds  <- rbind(dtt, dne)
ds$DIN2 <- factor(ds$DIN, levels = c("20", "0.1"), labels = c("N-replete","N-limited"))

# First filter pass to remove local minima
ds2 <- ds %>% 
  filter(case_when(
    DIN2 == "N-limited" ~ NC < 0.3 & NC > 0.05 & size > 0.0009 & size < 40 & mu > 3e-3,
    DIN2 == "N-replete" ~ NC < 0.5 & NC > 0.008 & size > 8 & size < 300 & mu > 0.001
  ))

# Label phenotypes by strategy
ds2 <- ds2 %>%
  mutate(strategy = case_when(
    Ea == 1.05 & ctag_max == 6e6   ~ "nutrient-efficient (NE)",
    Ea == 0.661 & ctag_max == 2.4e7 ~ "thermal-tolerant (TT)",
    TRUE                           ~ NA_character_))

# Checking for more local minima
ggplot(ds2, aes(x = temp, y = mu)) +
  facet_grid(DIN2 ~ strategy) +
  geom_line()

# Find Topt
dtopt <- ds2 %>% group_by(DIN2, strategy) %>% filter(mu == max(mu, na.rm = TRUE)) %>% select(DIN2, temp, mu, strategy)
ds_ratio <- dtopt %>% group_by(strategy) %>% summarise(delta_temp = temp[DIN2 == "N-limited"] - temp[DIN2 == "N-replete"])

arrows_ds2 <- data.frame(x = c(26, 24), xend = c(28, 21), y = c(0.9, 0.6), yend = c(0.9, 0.6), delta_temp = c(2, -3),
                        strategy = c("nutrient-efficient (NE)", "thermal-tolerant (TT)"))

prow2 <- ggplot(ds2, aes(x = temp, y = mu)) +
  geom_line(aes(linetype = factor(DIN2)), size = 1) +
  facet_wrap(~ strategy) +
  geom_point(data = dtopt, aes(x = temp, y = mu, fill = DIN2), size = 4, shape = 21) +
  scale_fill_manual(values = c("black", "white"), guide = "none") +
  scale_y_continuous(breaks = c(0.0, 0.8, 1.6), labels = scales::number_format(accuracy = 0.1), limits = c(0, 1.6)) +
  scale_x_continuous(breaks = c(10, 20, 30), labels = scales::number_format(accuracy = 1), limits = c(9, 33)) +
  ylab(expression(paste("Growth rate (", day^-1, ")", sep = ""))) +
  xlab(expression(paste("Temperature  (", degree, "C)"))) +
  scale_color_gradientn(
    colors = rev(brewer.pal(11, "RdBu")),
    limits = c(-7, 6),
    guide = "none") +
  geom_segment(
    data = arrows_ds2,
    aes(x = x, xend = xend, y = y, yend = yend, color = delta_temp),
    arrow = arrow(length = unit(0.35, "cm")),
    size = 2.5) +
  mytheme +
  theme(
    legend.position = c(0.1, 0.75),
    strip.text = element_text(size = 15),
    legend.title = element_blank(),
    legend.key = element_blank()
  )

prow2

# Plot final figure
pcombined <- prow1 / prow2 / prow3 / prow4 +
  plot_layout(heights = c(1, 0.8, 1.2, 1.2)) +
  plot_annotation(tag_levels = 'A')

pcombined <- pcombined & theme(plot.tag = element_text(size = 25))

ggsave("figure2.png", pcombined, device = "png", width = 9, height = 17, dpi = 600) 


