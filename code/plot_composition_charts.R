# Include libraries 
library(ggplot2)
library(ggpubr)
library(grid)
library(tidyverse)
library(shadowtext)

# Set theme
theme_set(theme_bw() + theme(legend.position="top"))

# Set paths to processed data and metadata directories
path_se <- "~/dev/figshare-repos/iba/processed_data/SE.v2/"
path_mg <- "~/dev/figshare-repos/iba/processed_data/MG.v2/"
path_metadata <- "~/dev/figshare-repos/iba/raw_data/"

# Include function for getting data
source("get_iba_co1_data_fxn.R")

# Define a cleaning function
clean <- function(X) {
    X <- X# Set theme[X$Phylum=="Arthropoda",]
    X <- X[!grepl("_X",X$Family),]
    X <- X[!grepl("unclassified",X$Family),]
    X <- X[!grepl("unresolved",X$Family),]
    X[X$Species!="Zoarces gillii",] 
}

# Define a function that makes a bar plot
make_barplot <- function(data, ntax_insects, ntax_other, max_x, div) {

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
    print(df)

    # Draw basic plot
    plt <- ggplot(df) + geom_col(aes(OTUs,Order), fill=BLUE, width=0.7)

    # Refine plot
    plt <- plt +
        scale_x_continuous(
            limits = c(0, max_x + max_x/50),
            breaks = seq(0, max_x, by = max_x/div),
            expand = c(0, 0), # The horizontal axis does not extend to either side
            position = "top"  # Labels are located on the top
        ) +
        # The vertical axis only extends upwards
        scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
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
            axis.title.y = element_blank(),
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
            data = subset(df, OTUs < max_x/5),
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
            data = subset(df, OTUs >= max_x/5),
            aes(0, y = Order, label = Order),
            hjust = 0,
            nudge_x = max_x/100,
            colour = "white",
            family = "Helvetica",
        size = 2.8
    )

    plt
}

foo <- function() {
# Get malaise data from Sweden
malaise_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="lysate|homogenate") |>
              clean()
lysate_se  <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="lysate") |>
              clean()
homogenate_se  <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="homogenate") |>
                  clean()

# Get malaise data from Madagascar
malaise_mg <- get_iba_co1_data(data_path=path_mg, metadata_path=path_metadata, country="mg",dataset="lysate") |>
              clean()

# Get soil and litter data for Sweden
soil_litter_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="soil|litter") |>
                  clean()
soil_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="soil") |>
           clean()
litter_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="litter") |>
             clean()

# Get litter data from Madagascar
litter_mg <- get_iba_co1_data(data_path=path_mg, metadata_path=path_metadata, country="mg",dataset="litter") |>
             clean()
}

# Make plots
pt1 <- make_barplot(malaise_se, 10, 5, max_x=12000, div=6)
pt2 <- make_barplot(malaise_mg, 10, 5, max_x=15000, div=5)
pt3 <- make_barplot(soil_litter_se, 10, 5, max_x=600, div=6)
pt4 <- make_barplot(litter_mg, 10, 5, max_x=700, div=7)

# Arrange fig7e
figure <- ggarrange(pt1, pt2, pt3, pt4, labels=c("A","B","C","D"), ncol=2, nrow=2)
pdf("Fig_X.pdf")
print(figure)
dev.off()

foo <- function () {
    # Refine plot
    plt <- plt +
        scale_x_continuous(
            limits = c(0, max_x + max_x/100),
            breaks = seq(0, max_x, by = max_x/10),
            expand = c(0, 0), # The horizontal axis does not extend to either side
            position = "top"  # Labels are located on the top
        ) +
        # The vertical axis only extends upwards
        scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
        theme(
            # Set background color to white
            panel.background = element_rect(fill = "white"),
            # Set the color and the width of the grid lines for the horizontal axis
            panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
            # Remove tick marks by setting their length to 0
            axis.ticks.length = unit(0, "mm"),
            # Remove the title for both axes
            axis.title = element_blank(),
            # Only left line of the vertical axis is painted in black
            axis.line.y.left = element_line(color = "black"),
            # Remove labels from the vertical axis
            axis.text.y = element_blank(),
            # But customize labels for the horizontal axis
            axis.text.x = element_text(family = "Helvetica", size = 16)
        )


plt <- plt + 
  geom_shadowtext(
    data = subset(data, count < 8),
    aes(count, y = name, label = name),
    hjust = 0,
    nudge_x = 0.3,
    colour = BLUE,
    bg.colour = "white",
    bg.r = 0.2,
    family = "Econ Sans Cnd",
    size = 7
  ) + 
  geom_text(
    data = subset(data, count >= 8),
    aes(0, y = name, label = name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Econ Sans Cnd",
    size = 7
  )

plt <- plt +
  labs(
    title = "Escape artists",
    subtitle = "Number of laboratory-acquired infections, 1970-2021"
  ) + 
  theme(
    plot.title = element_text(
      family = "Econ Sans Cnd", 
      face = "bold",
      size = 22
    ),
    plot.subtitle = element_text(
      family = "Econ Sans Cnd",
      size = 20
    )
  )

plt
}
