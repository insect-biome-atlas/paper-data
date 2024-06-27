# Set paths to processed data and metadata directories
path_se <- "~/dev/figshare-repos/iba/processed_data/SE.v2/"
path_mg <- "~/dev/figshare-repos/iba/processed_data/MG.v2/"
path_metadata <- "~/dev/figshare-repos/iba/raw_data/"

# Include function for getting data
source("get_iba_co1_data_fxn.R")

# Define a cleaning function
clean <- function(X) {
    X <- X[X$Phylum=="Arthropoda",]
    X <- X[!grepl("_X",X$Order),]
    X <- X[!grepl("unclassified",X$Order),]
    X <- X[!grepl("unresolved",X$Order),]
    X[X$Species!="Zoarces gillii",] 
}

# Get malaise data from Sweden
malaise_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="lysate|homogenate")
malaise_se <- clean(malaise_se)

lysates_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="lysate") |>
              clean()
homogenates_se <- get_iba_co1_data(data_path=path_se, metadata_path=path_metadata, country="se",dataset="homogenate") |>
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

# Compute order counts
T1 <- table(malaise_se$Order)
T2 <- table(malaise_mg$Order)
T3 <- table(soil_litter_se$Order)
T4 <- table(litter_mg$Order)


