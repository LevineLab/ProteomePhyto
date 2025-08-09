###################################
# Figures latitudinal simulations #
###################################

rm(list = ls())

library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(tidyr)
library(viridis)
library(ggpubr)
library(grid)

# Plotting settings ####
mytheme <-   theme(panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_rect(colour = "black", fill = "white", 
                                                   size = 1.5),
                   axis.text = element_text(size = 13, colour ="black"),
                   axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
                   axis.text.y = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),                   
                   axis.title = element_text(size = 17, colour ="black"),
                   axis.title.y = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                   axis.title.x = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                   axis.ticks.length = unit(-0.25, "cm"),
                   plot.margin = unit(c(1, 1.5, 0.5, 0.5), "lines"),
                   legend.position = c(0.15, 0.75),
                   legend.title = element_blank(),
                   legend.text = element_text(size=12, colour ="black"),
                   legend.key = element_blank(),
                   plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                   strip.background = element_rect(color = "black", fill = "white"),
                   strip.text = element_text(size = 15, color = "black")
)


###########################
# CESM environmental data #
###########################

cesm <- read.csv("model_CESM.csv")
cesm <- na.omit(cesm)

cesm_env <- na.omit(data.frame(temp = c(cesm$X20m_temperature_1980_2000, cesm$X20m_temperature_2080_2100),
                               nit = c(cesm$X20m_nitrate_1980_2000, cesm$X20m_nitrate_2080_2100),
                               time = c(rep("1980-2000", nrow(cesm)), rep("2080-2100", nrow(cesm))),
                               latitude = rep(cesm$latitude, 2)))

cesm$tdiff <- cesm$X20m_temperature_2080_2100 - cesm$X20m_temperature_1980_2000
cesm$tdiff <- (cesm$X20m_nitrate_1980_2000 - cesm$X20m_nitrate_2080_2100)/cesm$X20m_nitrate_1980_2000 * 100
mean(cesm$tdiff)
sd(cesm$tdiff)
min(cesm$tdiff)
max(cesm$tdiff)

p_lat_temp_cesm <- ggplot(cesm_env, aes(y = temp, x = latitude, linetype = time)) +
  geom_line(col = "black", size = 1) +
  scale_x_continuous(breaks = c(-70, -35, 0, 35 ,70), labels = scales::number_format(accuracy = 1), limits = c(-70,70), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 15, 30), labels = scales::number_format(accuracy = 1), limits = c(0,34), expand = c(0,0)) +
  ylab(expression("Temperature (" * degree * "C)")) +
  xlab("Latitude") +
  theme(legend.position = "none",
        plot.margin = margin(12, 5, 10, 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.0),
        axis.text = element_text(size = 13, colour ="black"),
        axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),                   
        axis.title = element_text(colour ="black"),
        axis.title.y = element_text(size = 13, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
        axis.title.x = element_text(size = 12, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
        axis.ticks.length = unit(-0.25, "cm"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.text = element_text(size=12, colour ="black"),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.text = element_text(size = 15, color = "black")) +
  coord_flip(clip = "off") 

p_lat_temp_cesm

p_lat_nit_cesm <- ggplot(cesm_env, aes(y = nit*1000, x = latitude, linetype = time)) +
  geom_line(col = "black", size = 1) +
  mytheme +
  ylab(expression("Nitrate (" * mu * "M)")) +
  xlab("") +
  scale_x_continuous(breaks = c(-70, -35, 0, 35 ,70), labels = scales::number_format(accuracy = 1), limits = c(-70,70), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 10, 20), labels = scales::number_format(accuracy = 1), limits = c(0,25), expand = c(0,0)) +
  theme(legend.position = "top",
        plot.margin = margin(12, 5, 10, 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", 
                                        size = 1.0),
        axis.text = element_text(size = 13, colour ="black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),                   
        axis.title = element_text(colour ="black"),
        axis.title.y = element_text(size = 13, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
        axis.title.x = element_text(size = 12, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
        axis.ticks.length = unit(-0.25, "cm"),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.text = element_text(size=12, colour ="black"),
        legend.key = element_blank(),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.text = element_text(size = 15, color = "black")) +
  coord_flip(clip = "off") 

p_lat_nit_cesm

pmenv <- p_lat_temp_cesm + p_lat_nit_cesm
pmenv

#################################################
# Mapping the latitudinal transect in the globe #
#################################################

library(sf)
library(ggplot2)
library(rnaturalearth)

# Load world map in sf format
world_sf <- ne_countries(scale = "medium", returnclass = "sf")

cesm <- data.frame(longitude = rep(-30, nrow(cesm)), 
                   latitude = cesm$latitude)

latitudes <- seq(-90, 90, length.out = 10000)
longitudes <- seq(-180, 180, length.out = 20000)

gtop <- data.frame(longitude = longitudes, 
                   latitude = 90)

gbottom <- data.frame(longitude = longitudes, 
                      latitude = -90)

gright <- data.frame(longitude = 180, 
                     latitude = latitudes)

gleft <- data.frame(longitude = -180, 
                    latitude = latitudes)

# Convert points to sf object and transform to Robinson projection
cesm_sf <- st_as_sf(cesm, coords = c("longitude", "latitude"), crs = 4326)
cesm_robin <- st_transform(cesm_sf, crs = "+proj=robin")

gtop_sf <- st_as_sf(gtop, coords = c("longitude", "latitude"), crs = 4326)
gtop_robin <- st_transform(gtop_sf, crs = "+proj=robin")

gbottom_sf <- st_as_sf(gbottom, coords = c("longitude", "latitude"), crs = 4326)
gbottom_robin <- st_transform(gbottom_sf, crs = "+proj=robin")

gright_sf <- st_as_sf(gright, coords = c("longitude", "latitude"), crs = 4326)
gright_robin <- st_transform(gright_sf, crs = "+proj=robin")

gleft_sf <- st_as_sf(gleft, coords = c("longitude", "latitude"), crs = 4326)
gleft_robin <- st_transform(gleft_sf, crs = "+proj=robin")

# Alternatively

bound <- st_sf(geometry = st_sfc(
  st_polygon(x = list(cbind(c(-180, rep(180, 100), rep(-180, 100)),
                            c(-90, seq(-90, 90, length = 100), 
                              seq(90, -90, length = 100))))),
  crs = 4326))


# Plot the map and points
pglobe <- ggplot(data = world_sf) +
  geom_sf(fill = "darkgray", color = "darkgray") + 
  geom_sf(data = cesm_robin, color = "red", size = 3) + 
  geom_sf(data = gtop_robin, color = "black", size = 5) +
  geom_sf(data = gbottom_robin, color = "black", size = 5) +
  geom_sf(data = gright_robin, color = "black", size = 3) +
  geom_sf(data = gleft_robin, color = "black", size = 3) +
  coord_sf(crs = "+proj=robin", expand = FALSE) +      
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),                            
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
pglobe

ggsave("p_globe.png", pglobe, device = "png", width = 11, height = 6, dpi = 600)


#########################
# CESM Latitudinal runs #
#########################

# Import model output for present and future conditions for each phenotype

dwf <- read.csv("output_warm_diatom_latitudinal_future.csv")   # Warm-diatom future
dwp <- read.csv("output_warm_diatom_latitudinal_present.csv")  # Warm-diatom present
dcf <- read.csv("output_cold_diatom_latitudinal_future.csv")   # Cold-diatom future
dcp <- read.csv("output_cold_diatom_latitudinal_present.csv")  # Cold-diatom present
dcyf <- read.csv("output_pro_latitudinal_future.csv")          # Prochlorococcus future
dcyp <- read.csv("output_pro_latitudinal_present.csv")         # Prochlorococcus present
dsf <- read.csv("output_syn_latitudinal_future.csv")           # Synechococcus future
dsp <- read.csv("output_syn_latitudinal_present.csv")          # Synechococcus present

dlat <- rbind(dwf, dwp, dcf, dcp, dcyf, dcyp, dsf, dsp)

dlat$model <- c(rep("future", nrow(dwf)),
                rep("present", nrow(dwp)),
                rep("future", nrow(dcp)),
                rep("present", nrow(dcf)),
                rep("future", nrow(dcyf)),
                rep("present", nrow(dcyp)),
                rep("future", nrow(dsf)),
                rep("present", nrow(dsp)))
dlat$phenotype <- c(rep("Warm-adapted diatom", 2*nrow(dwf)),
                    rep("Cold-adapted diatom", 2*nrow(dcp)),
                    rep("Prochlorococcus", 2*nrow(dcyf)),
                    rep("Synechococcus", 2*nrow(dsf)))
dlat$latitude <- rep(cesm$latitude, 8)

dlat$istress <- dlat$ilpb * dlat$itag

# Filter pass to remove non-optimal solutions
dlat <- dlat %>% filter(NC < 0.3 & NC > 0.05 & size > 0.6 & size < 28 & mu > 0.04 & mu < 0.75,
                        !(phenotype == "Prochlorococcus" & size > 6))

# Prochlrococcus does not survive below 10 degrees C, so exclude these rows from our data
dlat <- dlat %>% filter(!(phenotype == "Prochlorococcus" & temp < 10))

dlat$phenotype <- factor(dlat$phenotype, levels = c("Prochlorococcus", "Synechococcus", "Warm-adapted diatom", "Cold-adapted diatom"))

dlat$model <- factor(dlat$model, levels = c("present", "future"))

ggplot(dlat, aes(y = mu, x = latitude, col = factor(phenotype), fill = factor(phenotype), shape = factor(phenotype))) + 
  geom_line(size = 1) +
  geom_point(col = "black", size = 3) + 
  facet_wrap(~model) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  scale_shape_manual(values = c(21,22,24,23)) +
  coord_flip() +
  xlab("Latitude") +
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  ylim(0,0.7) +
  xlim(-70,70) +
  mytheme +
  theme(legend.position = "right", legend.direction = "vertical", legend.key = element_blank()) +
  guides(size = "none") 

# with the plot above, we can see non-optimal solutions, so we remove these manually:
dlat <- dlat %>%
  filter(!(phenotype == "Warm-adapted diatom" & model == "present" & latitude == -39.5),
         !(phenotype == "Warm-adapted diatom" & model == "present" & latitude == 5.5),
         !(phenotype == "Cold-adapted diatom" & model == "present" & latitude == -45.5),
         !(phenotype == "Cold-adapted diatom" & model == "future" & latitude == 35.5),
         !(phenotype == "Cold-adapted diatom" & model == "future" & latitude == 35.5),
         !(phenotype == "Cold-adapted diatom" & model == "future" & latitude == -70.5))

# Plot again to check if data points were correctly removed
ggplot(dlat, aes(y = mu, x = latitude, col = factor(phenotype), fill = factor(phenotype), shape = factor(phenotype))) + 
  geom_line(size = 1) +
  geom_point(col = "black", size = 3) + 
  facet_wrap(~model) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  scale_shape_manual(values = c(21,22,24,23)) +
  coord_flip() +
  xlab("Latitude") +
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  ylim(0,0.7) +
  xlim(-70,70) +
  mytheme +
  theme(legend.position = "right", legend.direction = "vertical", legend.key = element_blank()) +
  guides(size = "none") 

# Plot also cell diameter
ggplot(dlat, aes(y = size, x = latitude, col = factor(phenotype), fill = factor(phenotype), shape = factor(phenotype))) + 
  geom_line(size = 1) +
  geom_point(col = "black", size = 3) + 
  facet_wrap(~model) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  scale_shape_manual(values = c(21,22,24,23)) +
  coord_flip() +
  xlab("Latitude") +
  ylab("Cell diameter (um)") +
  scale_y_log10(limits = c(0.5, 15)) +
  mytheme +
  theme(legend.position = "right", legend.direction = "vertical", legend.key = element_blank()) +
  guides(size = "none") 

# with the plot above, we can see non-optimal solutions where Prochlorochococcus is unrealistically large, so we remove these:
dlat <- dlat %>%
  filter(!(phenotype == "Prochlorococcus" & model == "future" & latitude >= -5 & latitude <= 8))

# Now let's interpolate data across latitudes with the optimal solutions we have

# Generate a complete grid of latitude values with 0.5-degree intervals
latitude_seq <- seq(-75.5, 75.5, by = 0.5)

# Create a complete grid with all combinations of latitude and other factor variables
complete_grid <- expand.grid(
  latitude = latitude_seq,
  model = unique(dlat$model),  
  phenotype = unique(dlat$phenotype))

# Join the complete grid with your original data
dlat_expanded <- complete_grid %>%
  left_join(dlat, by = c("latitude", "model", "phenotype"))

# Interpolate missing values based on latitude
dlat_approx <- dlat_expanded %>%
  group_by(model, phenotype) %>%
  mutate(mu = zoo::na.approx(mu, na.rm = FALSE),
         itr = zoo::na.approx(itr, na.rm = FALSE),
         icbd = zoo::na.approx(icbd, na.rm = FALSE),
         ichl = zoo::na.approx(ichl, na.rm = FALSE),
         istress = zoo::na.approx(istress, na.rm = FALSE),
         size = zoo::na.approx(size, na.rm = FALSE))

dlat_approx2 <- dlat_approx %>%
  mutate(mu = ifelse(phenotype == "Cold-adapted diatom" & model == "present" & mu < 0.2, 0,
                     ifelse(phenotype == "Cold-adapted diatom" & model == "future" & latitude > -40 & latitude < 35 & mu < 0.4, 0, mu)))

# Let's make sure genera appear in italic in our plots
phen_labels <- c(
  "Prochlorococcus" = expression(italic("Prochlorococcus")),
  "Synechococcus"   = expression(italic("Synechococcus")),
  "Warm-adapted diatom" = "Warm-adapted diatom",
  "Cold-adapted diatom" = "Cold-adapted diatom")

master <- c("Prochlorococcus", "Synechococcus", "Warm-adapted diatom", "Cold-adapted diatom")

labels <- ifelse(master %in% c("Prochlorococcus", "Synechococcus"),
                 paste0("italic('", master, "')"),
                 paste0("'", master, "'"))

dlat_approx2 <- dlat_approx2 %>%
  mutate(
    phenotype = factor(phenotype, levels = master),
    phenotype_label = factor(
      case_when(
        phenotype == "Prochlorococcus" ~ labels[1],
        phenotype == "Synechococcus" ~ labels[2],
        phenotype == "Warm-adapted diatom" ~ labels[3],
        phenotype == "Cold-adapted diatom" ~ labels[4]
      ),
      levels = labels
    )
  )

ggplot(dlat_approx2, aes(y = mu, x = latitude, col = factor(phenotype), fill = factor(phenotype), shape = factor(phenotype))) + 
  geom_line(size = 2) +
  facet_wrap(~model) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF"), labels = phen_labels) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF"), labels = phen_labels) +
  xlab("Latitude") +
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6), labels = scales::number_format(accuracy = 0.1), limits = c(0,0.7), expand = c(0,0)) +
  scale_x_continuous(breaks = c(-70, -35, 0, 35 ,70), labels = scales::number_format(accuracy = 1), limits = c(-70,70), expand = c(0,0)) +
  theme(
    legend.position = "top",
    plot.margin = margin(20, 5, 10, 15),  
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white", size = 1.0),
    axis.text = element_text(size = 13, colour ="black"),
    axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y = element_blank(),                   
    axis.title = element_text(colour ="black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.text = element_text(size=12, colour ="black"),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, color = "black"),
    legend.background = element_rect(fill = NA, colour = NA)) +
  coord_flip(clip = "off")

# No optimal solutions were found for Prochlorococcus near the equator under future conditions because growth is predicted to be very low and size unrealistically high. Thus, let's set these as regions with growth = 0
dlat_approx2 <- dlat_approx2 %>% ungroup() %>% mutate(mu = if_else(phenotype == "Prochlorococcus" & between(latitude, -5, 8) & model == "future", NA_real_, mu),
                                      across(.cols = -c(latitude, phenotype, phenotype_label, model, mu),
                                             .fns = ~ replace(.x, phenotype == "Prochlorococcus" & between(latitude, -5, 8), NA_real_)))

# Now replace where growth is zero with NA values across the dataset so that we do not plot these in our figure
dlat_approx2 <- dlat_approx2 %>% ungroup() %>% mutate(across(.cols = -c(phenotype, model, phenotype_label, latitude),
                                             .fns = ~ replace(.x, mu == 0, NA_real_)))

pl2 <- ggplot(dlat_approx2, aes(y = mu, x = latitude, col = factor(phenotype), fill = factor(phenotype), shape = factor(phenotype), linetype = factor(phenotype))) + 
  geom_line(size = 1.2, na.rm = FALSE) +
  facet_wrap(~model) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF"), labels = phen_labels) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF"), labels = phen_labels) +
  scale_linetype_manual(values = c("solid", "solid", "22", "22"), labels = phen_labels) +
  xlab("Latitude") +
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6), labels = scales::number_format(accuracy = 0.1), limits = c(0,0.7), expand = c(0,0)) +
  scale_x_continuous(breaks = c(-70, -35, 0, 35 ,70), labels = scales::number_format(accuracy = 1), limits = c(-70,70), expand = c(0,0)) +
  theme(
    legend.position = "top",
    plot.margin = margin(20, 5, 10, 15),  
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white", size = 1.0),
    axis.text = element_text(size = 13, colour ="black"),
    axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y = element_blank(),                   
    axis.title = element_text(colour ="black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.text = element_text(size=12, colour ="black"),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, color = "black"),
    legend.background = element_rect(fill = NA, colour = NA)) +
  coord_flip(clip = "off")

pl2

pl2sep <- ggplot(dlat_approx2 %>% filter(model == "present"), aes(y = mu, x = latitude, col = factor(phenotype), fill = factor(phenotype), shape = factor(phenotype))) + 
  geom_line(size = 2) +
  facet_wrap(~ phenotype_label, nrow = 1, labeller = label_parsed) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  scale_fill_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  xlab("Latitude") +
  ylab(expression(paste("Growth rate (", day^-1,")", sep=""))) +
  scale_y_continuous(breaks = c(0.2, 0.5, 0.8), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.0), expand = c(0,0)) +
  scale_x_continuous(breaks = c(-70, -35, 0, 35 ,70), labels = scales::number_format(accuracy = 1), limits = c(-70,70), expand = c(0,0)) +
  theme(
    legend.position = "none",
    plot.margin = margin(20, 5, 10, 15),  # Increase the top margin for extra space
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white", size = 1.0),
    axis.text = element_text(size = 13, colour ="black"),
    axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
    axis.text.y = element_blank(),                   
    axis.title = element_text(colour ="black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.text = element_text(size=12, colour ="black"),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, color = "black")
  ) +
  coord_flip(clip = "off")
pl2sep

ggsave("figureS4.png", pl2sep, device = "png", width = 8, height = 4.5, dpi = 600)

# Now let's plot the phenotype with the highest growth rate at any given latitude
winner2 <- dlat_approx2 %>%
  group_by(latitude, model) %>%
  arrange(desc(mu), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

pl2b <- ggplot(winner2, aes(x = latitude, fill = phenotype)) +
  geom_bar(stat = "count", show.legend = FALSE, width = 0.8) +  
  facet_wrap(~model) +  
  scale_fill_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +  
  xlab("Latitude") +
  ylab("Dominant type") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(-70, -35, 0, 35 ,70), labels = scales::number_format(accuracy = 1), limits = c(-70,70), expand = c(0,0)) +
  theme(
    legend.position = "top",
    plot.margin = margin(40, 5, 10, 15),  
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white", size = 1.0),
    axis.text = element_text(size = 13, colour ="black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),                   
    axis.title = element_text(colour ="black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
    axis.ticks.length = unit(-0.25, "cm"),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.text = element_text(size=12, colour ="black"),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, color = "black", angle = 45)
  ) +
  coord_flip(clip = "off")

pl2b

z_lat2 <- p_lat_temp_cesm + p_lat_nit_cesm + pl2 + pl2b + plot_layout(widths = c(1, 1, 3, 0.5))

z_lat2

##################################################
# Relative proteome investments across latitudes #
##################################################

# Growth rates

ppmu <- ggplot(dlat_approx2, aes(y = mu, x = latitude, group = interaction(phenotype, model),linetype = factor(model), color = phenotype)) +
  geom_path(size = 1, na.rm = FALSE) +
  facet_wrap(~ phenotype, ncol = 1) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  coord_flip() +
  xlab("Latitude") +
  ylab(expression(paste("growth rate (", day^-1,")", sep=""))) +
  scale_x_continuous(breaks = c(-50, 0, 50), labels = scales::number_format(accuracy = 1), limits = c(-80,80)) +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.0)) +
  mytheme +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 10),
        strip.text = element_blank()) 
ppmu

# Transporters

pptr <- ggplot(dlat_approx2, aes(y = itr, x = latitude, group = interaction(phenotype, model),linetype = factor(model), color = phenotype)) +
  geom_path(size = 1, na.rm = FALSE) +
  facet_wrap(~ phenotype, ncol = 1) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  coord_flip() +
  xlab("Latitude") +
  ylab("transporters") +
  scale_x_continuous(breaks = c(-50, 0, 50), labels = scales::number_format(accuracy = 1), limits = c(-80,80)) +
  scale_y_log10() +
  mytheme +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank()) 
pptr

# Glycolysis

ppgl <- ggplot(dlat_approx2, aes(y = icbd, x = latitude, group = interaction(phenotype, model),linetype = factor(model), color = phenotype)) +
  geom_path(size = 1, na.rm = FALSE) +
  facet_wrap(~ phenotype, ncol = 1) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  coord_flip() +
  xlab("Latitude") +
  ylab("glycolysis") +
  scale_x_continuous(breaks = c(-50, 0, 50), labels = scales::number_format(accuracy = 1), limits = c(-80,80)) +
  scale_y_continuous(breaks = c(0.1, 0.15, 0.2), labels = scales::number_format(accuracy = 0.01), limits = c(0.08,0.21)) +
  mytheme +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank()) 
ppgl

# Heat-mitigation

pphm <- ggplot(dlat_approx2, aes(y = istress, x = latitude, group = interaction(phenotype, model),linetype = factor(model), color = phenotype)) +
  geom_path(size = 1, na.rm = FALSE) +
  facet_wrap(~ phenotype, ncol = 1) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  coord_flip() +
  xlab("Latitude") +
  ylab("heat-mitigation") +
  scale_x_continuous(breaks = c(-50, 0, 50), labels = scales::number_format(accuracy = 1), limits = c(-80,80)) +
  scale_y_log10() +
  mytheme +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank()) 
pphm

# Photosystems

ppph <- ggplot(dlat_approx2, aes(y = ichl, x = latitude, group = interaction(phenotype, model),linetype = factor(model), color = phenotype)) +
  geom_path(size = 1, na.rm = FALSE) +
  facet_wrap(~ phenotype_label, ncol = 1, strip.position = "right", labeller = label_parsed) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  coord_flip() +
  xlab("Latitude") +
  ylab("photosystems") +
  scale_x_continuous(breaks = c(-50, 0, 50), labels = scales::number_format(accuracy = 1), limits = c(-80,80)) +
  scale_y_continuous(breaks = c(0.00, 0.15, 0.3), labels = scales::number_format(accuracy = 0.01), limits = c(0,0.3)) +
  mytheme +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(size = 12),
        strip.background = element_blank()) 
ppph

# Cell diameter

ppsi <- ggplot(dlat_approx2, aes(y = size, x = latitude, group = interaction(phenotype, model),linetype = factor(model), color = phenotype)) +
  geom_path(size = 1, na.rm = FALSE) +
  facet_wrap(~ phenotype, ncol = 1) +
  scale_color_manual(values = c("#5AAE61", "#00441BFF","#762A83", "#C2A5CF")) +
  coord_flip() +
  xlab("Latitude") +
  ylab(expression(paste("cell diameter (", mu,"m)", sep=""))) +
  scale_x_continuous(breaks = c(-50, 0, 50), labels = scales::number_format(accuracy = 1), limits = c(-80,80)) +
  scale_y_log10() +
  mytheme +
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_blank()) 
ppsi

ggsave("figureS6.png", ppmu + ppsi + pptr + ppgl + pphm + ppph + plot_layout(nrow = 1), device = "png", width = 12 + 2.4, height = 9, dpi = 600) 

####################################################
# Changes in cell size summarized across latitudes #
####################################################

dlat_approx3 <- dlat_approx2 %>%
  mutate(phenotype2 = case_when(phenotype %in% c("Cold-adapted diatom", "Warm-adapted diatom") ~ "diatom",
                                phenotype == "Prochlorococcus" ~ "Pro",
                                phenotype == "Synechococcus" ~ "Syn"))

dlat_approx3 <- dlat_approx3 %>% filter(
  !(phenotype == "warm-adapted diatom" & (latitude > 50 | latitude < -50)),
  !(phenotype == "cold-adapted diatom" & latitude >= -50 & latitude <= 50))

selected_latitudes <- round(seq(-70, 70, length.out = 12))

# Filter the dataframe to keep only the rows where latitude matches the selected values
dlat_approx4 <- dlat_approx3 %>%
  filter(latitude %in% selected_latitudes) %>%
  select(phenotype2, latitude, size, model)

dlat_approx4$vol <- (4*pi*(dlat_approx4$size/2)^3)/3

# Compute the changes in cell size and the direction of change in cell size
dlat_approx5 <- dlat_approx4 %>%
  group_by(latitude, phenotype2) %>%
  mutate(
    size_future = size[model == "future"][1],
    size_present = size[model == "present"][1],
    size_diff = if (!is.na(size_future) && !is.na(size_present)) {
      size_future - size_present
    } else {
      NA_real_
    },
    size_change = case_when(
      size_diff > 0 ~ "increase",
      size_diff < 0 ~ "decrease",
      size_diff == 0 ~ "no change",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(-size_future, -size_present)

dlat_approx5$phenotype2 <- factor(dlat_approx5$phenotype2, levels = c("Pro", "Syn", "diatom"))

# Plot summary for cell size results
psize <- ggplot(dlat_approx5 %>% filter(latitude > -70 & model == "future"), aes(x = latitude, y = phenotype2, size = size, shape = size_change)) +
  facet_wrap(~ model) +
  geom_point() +
  coord_flip() +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = c(-70, -35, 0, 35 ,70), labels = scales::number_format(accuracy = 1), limits = c(-70,70), expand = c(0,0)) +
  scale_size_area(breaks = c(0.8, 3.0, 8.0),labels = c("0.8", "3.0", "8.0"), max_size = 8,
                  guide = guide_legend(override.aes = list(shape = 21, size = c(2.4, 4.5, 7)))) +
  scale_shape_manual(values = c(1, 19, 21), guide = FALSE) +
  theme(
    legend.position = "top",
    plot.margin = margin(40, 5, 10, 15),  
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white", size = 1.0),
    axis.text = element_text(size = 13, colour ="black"),
    axis.text.y = element_blank(), 
    axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 45, hjust = 1),                     
    axis.title = element_text(colour ="black"),
    axis.title.y = element_text(size = 15, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
    axis.title.x = element_text(size = 15, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.text = element_text(size=12, colour ="black"),
    legend.key = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, color = "black")) 

psize

pcombined <- p_lat_temp_cesm + p_lat_nit_cesm + pl2 + pl2b + psize + plot_layout(widths = c(0.9, 0.9, 2.2, 0.5, 0.9))
pcombined
ggsave("figure4.png", pcombined, device = "png", width = 9, height = 6.5, dpi = 600) 


##########################################################
# Relative proteome investments for poles versus tropics #
##########################################################

dlat_eq <- dlat_approx3 %>% 
  filter(latitude > -22 & latitude < 22, phenotype %in% c("Warm-adapted diatom", "Synechococcus", "Prochlorococcus")) %>%
  mutate(lat = "tropics")
dlat_eq %>% count(model, phenotype)

dlat_te <- dlat_approx3 %>% 
  filter(latitude < 70 & latitude > 50 | (latitude < -50 & latitude > -70), phenotype %in% c("Cold-adapted diatom", "Synechococcus")) %>%
  mutate(lat = "poles")
dlat_te %>% count(model, phenotype)

dlat_ter <- dlat_approx3 %>% 
  filter(latitude < 47 & latitude > 22 | (latitude > -32 & latitude < -12), phenotype %in% c("Warm-adapted diatom", "Prochlorococcus", "Synechococcus")) %>%
  mutate(lat = "temperate")
dlat_ter %>% count(model, phenotype)

dlat_eq_te <- rbind(dlat_eq, dlat_te, dlat_ter)

dlat_eq_te %>% count(model, phenotype)

dlat_eq_te <- dlat_eq_te %>%
  drop_na(itr, ichl, icbd, istress)

dlat_eq_te <- dlat_eq_te %>% group_by(latitude,model,phenotype,lat) %>% mutate(sumi = sum(itr, icbd, istress, ichl, na.rm = T),
                                                                               itr2 = itr/sumi,
                                                                               istress2 = istress/sumi,
                                                                               icbd2 = icbd/sumi,
                                                                               ichl2 = ichl/sumi,
                                                                               sumi2 = sum(itr2, istress2, icbd2, ichl2, na.rm = T))

dpie <- dlat_eq_te %>% group_by(model, phenotype, lat) %>% 
  summarize(mitr = mean(itr2), mistress = mean(istress2), micbd = mean(icbd2), michl = mean(ichl2)) %>%
  group_by(model, phenotype, lat) %>%
  mutate(total = sum(mitr, mistress,micbd,michl)) %>%
  pivot_longer(cols = c(mitr, mistress, micbd, michl), 
               names_to = "inv", 
               values_to = "value")

dpie$inv <- factor(dpie$inv, levels = c("mitr", "mistress", "micbd", "michl"), labels = c("transporters", "heat-mitigation", "respiration", "photosystems"))

dpie2 <- dpie %>%
  mutate(phenotype2 = case_when(
    phenotype %in% c("Warm-adapted diatom", "Cold-adapted diatom") ~ "diatom",
    TRUE ~ phenotype
  ))

# Create a new variable combining 'phenotype2' and 'inv'
dpie2 <- dpie2 %>%
  mutate(phenotype_inv = interaction(phenotype, inv))

# Create the plot

dpie2$lat <- factor(dpie2$lat, levels = c("tropics", "temperate", "poles"))
dpie2$phenotype2 <- factor(dpie2$phenotype2, levels = c("Prochlorococcus", "Synechococcus", "diatom"))

pinv <- ggplot(dpie2, aes(x = factor(model), y = value, fill = inv)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7, color = "black") +
  facet_grid(lat ~ phenotype2) +
  theme_bw() +
  ylab("Proteome investments") + 
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.0)) +
  xlab("") +
  mytheme +
  theme(legend.position = "right")

pinv

pinvf <- ggplot(dpie2, aes(x = factor(model), y = value, fill = inv)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.9, color = "black") +
  facet_grid(lat ~ phenotype2) +
  theme_bw() +
  ylab("Proteome investments") + 
  scale_fill_grey(start = 0.3, end = 0.8, 
                  guide = guide_legend(nrow = 4, ncol = 1)) +  
  scale_y_continuous(breaks = c(0.2, 0.5, 0.8), labels = scales::number_format(accuracy = 0.1), limits = c(0,1.0), expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  xlab("") +
  theme(
    legend.position = "top",
    plot.margin = margin(40, 5, 10, 15),  
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white", size = 1.0),
    axis.text = element_text(size = 13, colour ="black"),
    axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm"), angle = 45, hjust = 1), 
    axis.text.y = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),                     
    axis.title = element_text(colour ="black"),
    axis.title.y = element_text(size = 13, margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
    axis.title.x = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
    axis.ticks.length = unit(-0.25, "cm"),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.text = element_text(size=12, colour ="black"),
    legend.key = element_blank(),
    legend.key.size = unit(0.8, "lines"),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, color = "black")
  ) +
  coord_cartesian(clip = "off")

pinvf

ggsave("figureS7.png", pinvf, device = "png", width = 5, height = 6.5, dpi = 600) 
