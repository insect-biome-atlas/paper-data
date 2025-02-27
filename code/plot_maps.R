
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(MetBrewer)
library(terra)
library(tidyterra)
library(geodata)
library(scales)


# Sweden samples ----------------------------------------------------------

meta_se <- read_tsv("data/Sweden/sites_metadata_SE.tsv") 
meta_mg <- read_tsv("data/Madagascar/sites_metadata_MG.tsv")

# water bodies --------------------------------------------------------------------------------
# map of Sweden
sweden <- ne_countries(country = "sweden"  , scale = "large", returnclass = "sf") 

# GLobal lakes and rivers
lakes <-  ne_download(scale = 10, type = "lakes", category = "physical" , returnclass = "sf") %>%
          st_make_valid() %>% st_transform(st_crs(sweden))

# Swedish lakes 
lakes_sweden   <- st_intersection(lakes , sweden) 

# hill shaded elevation map -----------------------------------------------

# Get elevation, slope, and aspect data
swe_elev        <- elevation_30s("SWE" , path = ".") 
swe_slope       <- terrain(swe_elev, "slope", unit = "radians")
swe_aspect      <- terrain(swe_elev, "aspect", unit = "radians")

# Decide on which slopes to shade
swe_hshade        <- shade(swe_slope, swe_aspect, 30, 270)
names(swe_hshade) <- "shades" # add names to index variable
pal_greys <- hcl.colors(1000, "Grays")


# Get an index of which values to shade and to what degree
index <- swe_hshade %>%
          mutate(index_col = scales::rescale(shades, to = c(1, length(pal_greys)))) %>% 
          mutate(index_col = round(index_col)) %>%
          pull(index_col)

# Get a palette of grey colours to act as the "shade"
shade_cols <- pal_greys[index]

# Make the plot
swe_hshade_plot <- ggplot() +
  geom_spatraster(data = swe_hshade, fill = shade_cols, maxcell = Inf,alpha = 1)   # Avoid resampling with maxcell

# Scale min max values of elevation raster
elev_limits <- minmax(swe_elev) %>% as.vector() # Get min max
elev_limits <- c(floor(elev_limits[1] / 500), ceiling(elev_limits[2] / 500)) * 500 # Rounded to lower and upper 500
elev_limits <- pmax(elev_limits, 0) # Set min to 0




# same but for madagascar -------------------------------------------------

# Get elevation, slope, and aspect data
mad_elev        <- elevation_30s("MDG" , path = ".") 
mad_slope       <- terrain(mad_elev, "slope", unit = "radians")
mad_aspect      <- terrain(mad_elev, "aspect", unit = "radians")

# Decide on which slopes to shade
mad_hshade        <- shade(mad_slope, mad_aspect, 30, 270)
names(mad_hshade) <- "shades" # add names to index variable

# Get an index of which values to shade and to what degree
index <- mad_hshade %>%
          mutate(index_col = scales::rescale(shades, to = c(1, length(pal_greys)))) %>% 
          mutate(index_col = round(index_col)) %>%
          pull(index_col)

# Get a palette of grey colours to act as the "shade"
shade_cols <- pal_greys[index]

# Make the plot
mad_hshade_plot <- ggplot() +
  geom_spatraster(data = mad_hshade, fill = shade_cols, maxcell = Inf,alpha = 1)   # Avoid resampling with maxcell

# Scale min max values of elevation raster
elev_limits <- minmax(mad_elev) %>% as.vector() # Get min max
elev_limits <- c(floor(elev_limits[1] / 500), ceiling(elev_limits[2] / 500)) * 500 # Rounded to lower and upper 500
elev_limits <- pmax(elev_limits, 0) # Set min to 0


# plots ---------------------------------------------------------------------------------------

# Format meta-data
meta_mg <- meta_mg  |>
            mutate(trap_habitat=recode(trap_habitat , 
                                       "Dry_Forest"          = "Dry Forest",
                                        "Montane_Rainforest" = "Wet Forest",
                                        "Rainforest"         = "Wet Forest"),
                              malaise_trap_type = recode(malaise_trap_type , "Single_trap" = "Single trap"))

# Format meta-data
meta_se <- meta_se  |>
            mutate(malaise_trap_type = recode(malaise_trap_type , "Single_trap" = "Single trap"),
                   trap_habitat = str_to_title(trap_habitat)) |> 
            drop_na(NILS_mhabitat)

tcols <- met.brewer("Veronese" , n = length(unique(meta_se$NILS_mhabitat)))

p1 <- swe_hshade_plot + # PLot hillshaded map
  geom_spatraster(data = swe_elev,na.rm = TRUE, maxcell = Inf, show.legend = TRUE) +
  ggspatial::annotation_scale(location = 'tl',width_hint = .6,text_cex = 2)+
  scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.4,direction = 1)+ # Add tinted alpha layer to hill shaded layer
  geom_sf(data=lakes_sweden , fill  = "blue", alpha = .3 ,colour = NA)+
  geom_point(data=meta_se , aes(longitude_WGS84 , latitude_WGS84 , 
                                shape = malaise_trap_type),
                                size=6, colour = "white")+
  geom_point(data=meta_se , aes(longitude_WGS84 , latitude_WGS84 , 
                                colour = trap_habitat, shape = malaise_trap_type),
                              size=5)+
  scale_colour_manual(values = tcols)+
  theme_linedraw(base_size = 25)+
  labs(x = "Longitude" ,
                        y = "Latitude" , 
                        fill = "Elevation (m)",
                        shape = "Trap-type",
                        colour = "Trap habitat")+
  theme(legend.key = element_rect(fill = "white"),
        legend.position = "left",
        legend.direction = "vertical",
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(1.5, "cm"))+
  guides(color=guide_legend(override.aes=list(fill=NA,size=6)),
         shape=guide_legend(override.aes=list(fill=NA,size=6)))

p2 <- mad_hshade_plot + # PLot hillshaded map
  geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = TRUE) +
  ggspatial::annotation_scale(location = 'tl',width_hint = .6,text_cex = 2)+
  geom_point(data=meta_mg , aes(longitude_WGS84 , latitude_WGS84 , 
                                shape = malaise_trap_type),
             size=6, colour = "white")+
  geom_point(data=meta_mg , aes(longitude_WGS84 , latitude_WGS84 , 
                                colour = trap_habitat, shape = malaise_trap_type),size=5)+ 
  scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.4,direction = 1)+
  scale_colour_viridis_d(option="mako", end=.8)+
  theme_linedraw(base_size = 25)+
  labs(x = "Longitude" ,
       y = "Latitude" , 
       fill = "Elevation",
       shape = "Trap-type",
       colour = "Trap habitat")+
  theme(legend.key = element_rect(fill = "white"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(1.5, "cm"))+
  guides(color=guide_legend(override.aes=list(fill=NA,size=6)),
         shape=guide_legend(override.aes=list(fill=NA,size=6)))



# plot ----------------------------------------------------------------------------------------
tiff("figures/maps.tiff", width = 2000, height = 1000, compression = "lzw")
p1 + p2  & labs(x = "Longitude" , y = "Latitude" , fill = "Elevation", shape = "Trap-type") 
dev.off()

browseURL("figures/maps.tiff")



