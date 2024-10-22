# Include libraries 
library(ggplot2)
library(ggpubr)
library(grid)
library(tidyverse)
library(shadowtext)

# Set theme
theme_set(theme_bw() + theme(legend.position="top"))

# Set paths to processed data and metadata directories
path_se <- "~/dev/figshare-repos/iba/processed_data/v2.1_cleaned_noise_filtered/"
path_mg <- "~/dev/figshare-repos/iba/processed_data/v2.1_cleaned_noise_filtered/"
path_metadata <- "~/dev/figshare-repos/iba/raw_data/v4/"

# Include function for getting data
old_dir <- getwd()
setwd("~/dev/ms-repos-iba/utils/")
source("get_iba_co1_data_fxn.R")
setwd(old_dir)

# Define a cleaning function

clean <- function(X) {
    X <- X[X$Phylum=="Arthropoda",]
#    X <- X[!grepl("_X",X$Family),]
#    X <- X[!grepl("unclassified",X$Family),]
#    X <- X[!grepl("unresolved",X$Family),]
    X[X$Species!="Zoarces gillii",] 
}

# Define a function that makes a bar plot
make_barplot <- function(data, ntax_insects, ntax_other, max_tick, div) {

    BLUE <- "#076fa2"
    RED <- "#E3120B"
    BLACK <- "#202020"
    GREY <- "grey50"

    # Make a data frame based on number of OTUs and add Class
    df <- data.frame(table(data$Order))
    colnames(df) <- c("Order","OTUs")
    df$Class <- data$Class[match(df$Order,data$Order)]

    # Make a new data frame with the top 10 among hexapods and other arthropods
    df_ins <- df[df$Class %in% c("Insecta","Collembola","Diplura","Protura"),]
    df_oth <- df[!(df$Class %in% c("Insecta","Collembola","Diplura","Protura")),]

    df_ins <- df_ins[order(df_ins$OTUs, decreasing=TRUE),]
    df_oth <- df_oth[order(df_oth$OTUs, decreasing=TRUE),]

    df_ins <- df_ins[1:ntax_insects,]
    df_ins <- df_ins[order(df_ins$OTUs),]

    df_oth <- df_oth[1:ntax_other,]
    df_oth <- df_oth[order(df_oth$OTUs),]

    df <- rbind(df_oth,df_ins)

    df$Order <- factor(df$Order, levels = df$Order)
    df$y <- seq(length(df$Order)) * 0.9
    max_x <- max(df$OTUs)
    df$label_cutoff <- max_x * as.numeric(sapply(as.character(df$Order),nchar)) / 46

    # Draw basic plot
    plt <- ggplot(df) + geom_col(aes(OTUs,Order), fill=BLUE, width=0.7)

    # Refine plot
    plt <- plt +
        scale_x_continuous(
            limits = c(0, max_x),
            breaks = seq(0, max_tick, by = max_tick/div),
            position = "top"  # Labels are located on the top
        ) +
        # The vertical axis extends upwards and downwards
        scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
        geom_hline(yintercept=5.5, linewidth=0.3, linetype="dashed", color=BLACK) +
        annotate("text", x=0.88*max_x, y=6.5, label="Hexapods", size=4, color=BLACK) +
        annotate("text", x=0.80*max_x, y=1.5, label="Other arthropods", size=4, color=BLACK) +
        theme(
            # Set background color to white
            panel.background = element_rect(fill = "white"),
            # Set the color and the width of the grid lines for the horizontal axis
            panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
            panel.grid.minor.x = element_line(color = "white", linewidth = 0.3),
            panel.grid.major.y = element_line(color = "white", linewidth = 0.3),
            # Remove tick marks by setting their length to 0
            axis.ticks.length = unit(0, "mm"),
            # Remove the title for the y axis, make room for but do not show x axis label
            axis.title.y = element_text(color="white"),
            axis.title.x = element_text(color="white"),
            # Only left line of the vertical axis is painted in black
            axis.line.y.left = element_line(color = "black"),
            # Remove labels from the vertical axis
            axis.text.y = element_blank(),
            # But customize labels for the horizontal axis
            axis.text.x = element_text(family = "Helvetica", size = 8)
        )

    # Add labels back in, into the plot
    plt <- plt + 
        geom_shadowtext(
            data = subset(df, OTUs < label_cutoff),
            aes(OTUs, y = Order, label = Order),
            hjust = 0,
            nudge_x = max_x/100,
            colour = "black",
            bg.colour = "white",
            bg.r = 0.2,
            family = "Helvetica",
            size = 2.8
        ) + 
        geom_text(
            data = subset(df, OTUs >= label_cutoff),
            aes(0, y = Order, label = Order),
            hjust = 0,
            nudge_x = max_x/100,
            colour = "white",
            family = "Helvetica",
            size = 2.8
        )

    plt
}

# Get malaise data from Sweden
malaise_se <- get_iba_co1_data(data_path=path_se,
                               metadata_path=path_metadata,
                               country="SE",
                               dataset="lysate|homogenate",
                               remove_spikes=FALSE) |>  # TODO: Temporary fix because spike-ins are missing
              clean()
#lysate_se  <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="SE",dataset="lysate") |>
#              clean()
#homogenate_se  <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="SE",dataset="homogenate") |>
#                  clean()

# Get malaise data from Madagascar
malaise_mg <- get_iba_co1_data(data_path=path_mg,
                               metadata_path=path_metadata,
                               country="MG",
                               dataset="lysate",
                               remove_spikes=FALSE) |>  # TODO: Temporary fix because spike-ins are missing
              clean()

# Get soil and litter data for Sweden
soil_litter_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="SE",dataset="soil|litter") |>
                  clean()
#soil_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="SE",dataset="soil") |>
#           clean()
#litter_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="SE",dataset="litter") |>
#             clean()

# Get litter data from Madagascar
litter_mg <- get_iba_co1_data(data_path=path_mg, metadata_path=path_metadata, country="MG",dataset="litter") |>
             clean()

# Make plots
pt1 <- make_barplot(malaise_se, 10, 5, max_tick=12000, div=6)
pt2 <- make_barplot(malaise_mg, 10, 5, max_tick=30000, div=6)
pt3 <- make_barplot(soil_litter_se, 10, 5, max_tick=600, div=6)
pt4 <- make_barplot(litter_mg, 10, 5, max_tick=1200, div=6)

# Arrange plots and generate figure
figure <- ggarrange(pt1, pt2, pt3, pt4, labels=c("A","B","C","D"), ncol=2, nrow=2)
ggexport(figure,filename="../figures/composition_barcharts.pdf")

