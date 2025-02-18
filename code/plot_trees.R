# Load required libraries
library(tidyverse)
library(data.table)
library(tidytree)
library(ape)
library(rotl)
library(taxize)
library(ggtree)
library(ggstar)
library(ggtreeExtra)
library(ggnewscale)
library(scico)
library(scales)

# Set random seed
set.seed(10)

# ---------------------------------------------------------------------------------------------
# data ----------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# All data must be 


# --------------- Spike in IDs 
spike_ins_mg  <- fread("data/Biological_spikes_MG_taxonomy.tsv")
spike_ins_se  <- fread("data/Biological_spikes_SE_taxonomy.tsv")

# --------------- Get sample IDs 
se_meta <- fread("data/co1_sequencing_metadata_se.tsv")
mg_meta <- fread("data/co1_sequencing_metadata_mg.tsv")

ID_se_malaise <- se_meta |> 
                _[str_detect(dataset, "homogenate|lysate"),] |> 
                _[lab_sample_type == "sample"] |> 
                _$sampleID_NGI

ID_mg_malaise <- mg_meta |> 
                _[str_detect(dataset, "homogenate|lysate"),] |> 
                _[lab_sample_type == "sample"] |> 
                _$sampleID_NGI

ID_se_soil <- se_meta |> 
              _[str_detect(dataset, "soil"),] |> 
              _[lab_sample_type == "sample"] |> 
              _$sampleID_NGI

ID_mg_soil <- mg_meta |> 
              _[str_detect(dataset, "litter"),] |> 
              _[lab_sample_type == "sample"] |> 
              _$sampleID_NGI


# ----------------- Get and clean cluster data
# These files need downloading from figshare

swe_clusters <- fread("data/non_numts_cleaned_SE.tsv") |>
                _[representative == 1 ,] |> 
                _[Phylum == "Arthropoda" ,] |> 
                _[!str_detect(Family , "_X*"),] |> 
                _[!str_detect(Species , paste(spike_ins_se$Species,collapse="|")),]

mg_clusters  <- fread("data/non_numts_cleaned_MG.tsv") |>
                _[representative == 1 ,] |> 
                _[Phylum == "Arthropoda" ,] |> 
                _[!str_detect(Family , "_X*"),] |> 
                _[!str_detect(Species , paste(spike_ins_mg$Species,collapse="|")),]


# ------------ Get read numbers 
swe_counts   <- fread("data/non_numts_cluster_counts_cleaned_SE.tsv")
mg_counts    <- fread("data/non_numts_cluster_counts_cleaned_MG.tsv")


# ---------------------------------------------------------------------------------------------
# tidy data ----------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# tidy counts & cluster data ------------------------------------------------------------------

# First get counts of each species in each sample
swe_counts_long <- swe_counts |> 
                     _[cluster %in% swe_clusters$cluster,] |> 
                    melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
                    _[read_count > 10,] 
                    

# Merge with cluster IDS
se_counts_DT <- swe_counts_long[swe_clusters, on = "cluster"] |> 
                _[!str_detect(Order , "unclassified"),] 


# First get counts of each species in each sample
mg_counts_long <- mg_counts |> 
                  _[cluster %in% mg_clusters$cluster,] |> 
                  melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
                  _[read_count > 10,] 

# Merge with cluster IDS
mg_counts_DT <- mg_counts_long[mg_clusters, on = "cluster"] |>
               _[!str_detect(Order , "unclassified"),] 


# Swedish malaise data -----------------------------------------------------------------------

malaise_counts_se <- se_counts_DT |>
                     _[sampleID_NGI %in% ID_se_malaise,]

# Summarise to order level
malaise_order_stats_se <- malaise_counts_se |>
                        group_by(Order) |>
                        summarise(total_reads = sum(read_count , na.rm=TRUE) , 
                                  total_OTU = n_distinct(cluster)) 

malaise_orders_se <- malaise_order_stats_se$Order

# malagasy malaise data -----------------------------------------------------------------------

malaise_counts_mg <- mg_counts_DT |>
                     _[sampleID_NGI %in% ID_mg_malaise,]

# Summarise to order level
malaise_order_stats_mg <- malaise_counts_mg |> droplevels() |> 
                           group_by(Order) |>
                           summarise(total_reads = sum(read_count , na.rm=TRUE) , 
                                     total_OTU = n_distinct(cluster)) |> 
                           filter(total_reads>0)

malaise_orders_mg <- unique(malaise_order_stats_mg$Order)

# Swedish soil data ---------------------------------------------------------------------------

soil_counts_se <- se_counts_DT |> 
                  _[sampleID_NGI %in% ID_se_soil,]

# Summarise to order level
soil_order_stats_se <- soil_counts_se |>
   group_by(Order) |>
   summarise(total_reads = sum(read_count , na.rm=TRUE) , 
             total_OTU = n_distinct(cluster)) 

soil_orders_se <- soil_order_stats_se$Order

# malagasy soil data --------------------------------------------------------------------------

soil_counts_mg <- mg_counts_DT |> 
                  _[sampleID_NGI %in% ID_mg_soil,]

# Summarise to order level
soil_order_stats_mg <- soil_counts_mg |>
   group_by(Order) |>
   summarise(total_reads = sum(read_count , na.rm=TRUE) , 
             total_OTU = n_distinct(cluster)) 

soil_orders_mg <- soil_order_stats_mg$Order

# ---------------------------------------------------------------------------------------------
# Get trees  ----------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# functions -----------------------------------------------------------------------------------
# function to fetch order IDs from ncbi
get_tree <- function(orders){
   IBA_ids        <- get_ids(orders, db = "ncbi")
   IBA_class      <- classification(IBA_ids, db = "ncbi")
   dupes          <- duplicated(names(IBA_class[[1]]))
   IBA_class[[1]] <- IBA_class[[1]][!dupes]
   tree           <- taxize::class2tree(IBA_class$ncbi, check = TRUE)
}
# ---------------------------------------------------------------------------------------------

# Get trees
malaise_tree_se <- get_tree(orders = malaise_orders_se)
malaise_tree_mg <- get_tree(orders = malaise_orders_mg)
soil_tree_se    <- get_tree(orders = soil_orders_se)
soil_tree_mg    <- get_tree(orders = soil_orders_mg)
# 
# saveRDS(malaise_tree_se , "data/malaise_tree_se.rds")
# saveRDS(malaise_tree_mg , "data/malaise_tree_mg.rds")
# saveRDS(soil_tree_se , "data/soil_tree_se.rds")
# saveRDS(soil_tree_mg , "data/soil_tree_mg.rds")


# ---------------------------------------------------------------------------------------------
# Prepare trees for plotting ------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# Sweden - malaise ----------------------------------------------------------------------------

# Read tree
#malaise_tree_se <- readRDS("data/malaise_tree_se.rds")

# Prepare tree data frame
malaise_tree_df_se <- as_data_frame(malaise_tree_se$phylo) %>% 
   mutate(Order = label) %>%
   left_join(malaise_order_stats_se, by = "Order")

# Convert tree data frame back to phylo object
malaise_tree_se$phylo <- as.phylo(malaise_tree_df_se)

# Remove nodes with no reads
malaise_noreads_se <- malaise_tree_df_se %>% filter(is.na(total_reads) , total_OTU<10) %>% pull(label)
malaise_tree_se$phylo <- drop.tip(malaise_tree_se$phylo, malaise_noreads_se)

# Define clade labels, numbers, and colors
#check tips and nodes
malaise_tree_se$phylo$node.label
malaise_tree_se$phylo$tip.label


malaise_clades_se <- c(2, 4,6, 10,11)
malaise_clade_labels_se <- malaise_tree_se$phylo$node.label[malaise_clades_se]
malaise_clade_labels_se[4] <- "Insecta" 
malaise_clade_numbers_se <- length(malaise_tree_se$phylo$tip.label) + malaise_clades_se
malaise_clade_cols_se <- scico(length(malaise_clades_se), palette = "nuuk")

malaise_highlight_df_se <- data.frame(id = malaise_clade_numbers_se, Clade = malaise_clade_labels_se)


# Madagascar - malaise ----------------------------------------------------------------------------------

# Read tree
#malaise_tree_mg <- readRDS("data/malaise_tree_mg.rds")

# Prepare tree data frame
malaise_tree_df_mg <- as_data_frame(malaise_tree_mg$phylo) %>% 
   mutate(Order = label) %>%
   left_join(malaise_order_stats_mg, by = "Order")

# Convert tree data frame back to phylo object
malaise_tree_mg$phylo <- as.phylo(malaise_tree_df_mg)

# Remove nodes with no reads & low OTUs
malaise_noreads_mg <- malaise_tree_df_mg %>%
                     filter(is.na(total_reads) , total_OTU < 10) %>%
                     pull(label)

malaise_tree_mg$phylo <- drop.tip(malaise_tree_mg$phylo, malaise_noreads_mg)

# Define clade labels, numbers, and colors
#check tips and nodes
malaise_tree_mg$phylo$node.label
malaise_tree_mg$phylo$tip.label

malaise_clades_mg <- c(2,4, 6, 9,10)
malaise_clade_labels_mg <- malaise_tree_mg$phylo$node.label[malaise_clades_mg]
malaise_clade_numbers_mg <- length(malaise_tree_mg$phylo$tip.label) + malaise_clades_mg

malaise_highlight_df_mg <- data.frame(id = malaise_clade_numbers_mg, Clade = malaise_clade_labels_mg)


# Sweden  - soil ------------------------------------------------------------------------------


#soil_tree_se <- readRDS("data/soil_tree_se.rds")

soil_tree_df_se <- as_data_frame(soil_tree_se$phylo) %>% 
   mutate(Order = label) %>%
   left_join(soil_order_stats_se, by = "Order")

# Convert tree data frame back to phylo object
soil_tree_se$phylo <- as.phylo(soil_tree_df_se)

# Remove nodes with no reads
soil_noreads_se <- soil_tree_df_se %>% filter(is.na(total_reads) , total_OTU<10) %>% pull(label)
soil_tree_se$phylo <- drop.tip(soil_tree_se$phylo, soil_noreads_se)

# Define clade labels, numbers, and colors
#check tips and nodes
soil_tree_se$phylo$node.label
soil_tree_se$phylo$tip.label

soil_clades_se <- c(2, 4,7, 11, 12)
soil_clade_labels_se <- soil_tree_se$phylo$node.label[soil_clades_se]
soil_clade_numbers_se <- length(soil_tree_se$phylo$tip.label) + soil_clades_se

soil_highlight_df_se <- data.frame(id = soil_clade_numbers_se, Clade = soil_clade_labels_se)

# Madagascar litter  --------------------------------------------------------------------------
#soil_tree_mg <- readRDS("data/soil_tree_mg.rds")

soil_tree_df_mg <- as_data_frame(soil_tree_mg$phylo) %>% 
   mutate(Order = label) %>%
   left_join(soil_order_stats_mg, by = "Order")

# Convert tree data frame back to phylo object
soil_tree_mg$phylo <- as.phylo(soil_tree_df_mg)

# Remove nodes with no reads
soil_noreads_mg <- soil_tree_df_mg %>% filter(is.na(total_reads) , total_OTU<10) %>% pull(label)
soil_tree_mg$phylo <- drop.tip(soil_tree_mg$phylo, soil_noreads_mg)

# Define clade labels, numbers, and colors
#check tips and nodes
soil_tree_mg$phylo$node.label
soil_tree_mg$phylo$tip.label

soil_clades_mg <- c(2, 6,7)
soil_clade_labels_mg <- soil_tree_mg$phylo$node.label[soil_clades_mg]
soil_clade_numbers_mg <- length(soil_tree_mg$phylo$tip.label) + soil_clades_mg

soil_highlight_df_mg <- data.frame(id = soil_clade_numbers_mg, Clade = soil_clade_labels_mg)

# ---------------------------------------------------------------------------------------------
# make the plots ------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------


 # plotting pars -------------------------------------------------------------------------------
bs                <- 5
margins           <-  margin(5,1,8,1, "cm")
tip_lab_size      <- 7
legend_title_size <- 20
legend_text_size  <- 20 
legend.position   <- c(.5 , -.2)

 nOTU_lims   <- log10(range(malaise_tree_df_mg$total_OTU , malaise_tree_df_se$total_OTU , na.rm = TRUE)+1)
 nReads_lims <- log10(range(malaise_tree_df_mg$total_reads , malaise_tree_df_se$total_reads, na.rm = TRUE)+1)

 # plots ---------------------------------------------------------------------------------------
 
 #Clade colours
 
 malaise_clades <- union(malaise_clade_labels_mg , malaise_clade_labels_se)
 malaise_clade_cols <- scico(length(malaise_clades), palette = "nuuk")
 
 
 #Sweden malaise
 p1 <-  ggtree(malaise_tree_se$phylo, layout = "circular") +
   geom_tiplab(offset = 10 ,geom="text", size = tip_lab_size)+
   geom_highlight( data = malaise_highlight_df_se, mapping = aes(fill =Clade ,  node =id ) )+
   scale_fill_manual(values = malaise_clade_cols)+
   new_scale(new_aes = "fill") +
   geom_fruit(data = malaise_tree_df_se,
              geom = geom_tile,
              mapping = aes(y = Order, fill = log10(total_reads+1)), 
              width = 4, offset = .05, colour = "white")+
   scale_fill_viridis_c(option = "mako", limits = nReads_lims,
                        guide = guide_colorbar(title.position = "top",
                                               byrow = TRUE)) +
   labs(fill = "log10(reads)") +
   new_scale(new_aes = "fill") +
   geom_fruit(data = malaise_tree_df_se,
              geom = geom_tile,
              mapping = aes(y = Order, fill = log10(total_OTU+1)), 
              width = 4, offset = .1, colour = "white") +
   scale_fill_viridis_c(option = "rocket", limits = nOTU_lims,
                        guide = guide_colorbar(title.position = "top",
                                               byrow = TRUE))+
   theme_void(base_size = bs) +
   theme(legend.text = element_text(size = legend_text_size), 
         legend.title = element_text(size = legend_title_size,face = "bold"),
         legend.title.align=0.5,
         legend.position = legend.position,
         legend.title.position = "top",
         legend.spacing.x = unit(1, "cm"),
         legend.box="horizontal", legend.margin=margin(),
         plot.margin = margins)+
   labs(fill = "log10(nOTU)") 
 
  
 # MG malaise
 p2 <- ggtree(malaise_tree_mg$phylo, layout = "circular") +
    geom_tiplab(offset = 10, geom = "text", size = tip_lab_size) +
    geom_highlight(data = malaise_highlight_df_mg, mapping = aes(fill = Clade, node = id)) +
    scale_fill_manual(values = malaise_clade_cols) +
    new_scale(new_aes = "fill") +
    geom_fruit(data = malaise_tree_df_mg,
               geom = geom_tile,
               mapping = aes(y = Order, fill = log10(total_reads)),
               width = 4, offset = .05, colour = "white") +
    scale_fill_viridis_c(option = "mako", limits = nReads_lims,
                         guide = guide_colorbar(title.position = "top",
                                                byrow = TRUE)) +
    labs(fill = "log10(reads)") +
    new_scale(new_aes = "fill") +
    geom_fruit(data = malaise_tree_df_mg,
               geom = geom_tile,
               mapping = aes(y = Order, fill = log10(total_OTU+1)),
               width = 4, offset = .1, colour = "white") +
    scale_fill_viridis_c(option = "rocket", limits = nOTU_lims,
                         guide = guide_colorbar(title.position = "top",
                                                byrow = TRUE)) +
    theme_void(base_size = bs) +
    theme(legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_title_size, face = "bold"),
          legend.title.align = 0.5,
          legend.position = legend.position,
          legend.title.position = "top",
          legend.spacing.x = unit(1, "cm"),
          legend.box = "horizontal", legend.margin = margin(),
          plot.margin = margins) +
    labs(fill = "log10(nOTU)")
 
 

# soil ----------------------------------------------------------------------------------------

 soil_clades <- union(soil_clade_labels_mg, soil_clade_labels_se)
 soil_clade_cols <- scico(length(soil_clades), palette = "nuuk")
 
 # Sweden soil
 p3 <- ggtree(soil_tree_se$phylo, layout = "circular") +
    geom_tiplab(offset = 10, geom = "text", size = tip_lab_size) +
    geom_highlight(data = soil_highlight_df_se, mapping = aes(fill = Clade, node = id)) +
    scale_fill_manual(values = soil_clade_cols) +
    new_scale(new_aes = "fill") +
    geom_fruit(data = soil_tree_df_se,
               geom = geom_tile,
               mapping = aes(y = Order, fill = log10(total_reads + 1)), 
               width = 4, offset = .05, colour = "white") +
    scale_fill_viridis_c(option = "mako", limits = nReads_lims,
                         guide = guide_colorbar(title.position = "top",
                                                byrow = TRUE)) +
    labs(fill = "log10(reads)") +
    new_scale(new_aes = "fill") +
    geom_fruit(data = soil_tree_df_se,
               geom = geom_tile,
               mapping = aes(y = Order, fill = log10(total_OTU + 1)), 
               width = 4, offset = .1, colour = "white") +
    scale_fill_viridis_c(option = "rocket", limits = nOTU_lims,
                         guide = guide_colorbar(title.position = "top",
                                                byrow = TRUE)) +
    theme_void(base_size = bs) +
    theme(legend.text = element_text(size = legend_text_size), 
          legend.title = element_text(size = legend_title_size, face = "bold"),
          legend.title.align = 0.5,
          legend.position = legend.position,
          legend.title.position = "top",
          legend.spacing.x = unit(1, "cm"),
          legend.box = "horizontal", legend.margin = margin(),
          plot.margin = margins) +
    labs(fill = "log10(nOTU)")

 
 mg_clade_cols <- soil_clade_cols[c(1,2,4)]
 p4 <- ggtree(soil_tree_mg$phylo, layout = "circular") +
    geom_tiplab(offset = 10, geom = "text", size = tip_lab_size) +
    geom_highlight(data = soil_highlight_df_mg, mapping = aes(fill = Clade, node = id)) +
    scale_fill_manual(values = mg_clade_cols) +
    new_scale(new_aes = "fill") +
    geom_fruit(data = soil_tree_df_mg,
               geom = geom_tile,
               mapping = aes(y = Order, fill = log10(total_reads + 1)), 
               width = 4, offset = .05, colour = "white") +
    scale_fill_viridis_c(option = "mako", limits = nReads_lims,
                         guide = guide_colorbar(title.position = "top",
                                                byrow = TRUE)) +
    labs(fill = "log10(reads)") +
    new_scale(new_aes = "fill") +
    geom_fruit(data = soil_tree_df_mg,
               geom = geom_tile,
               mapping = aes(y = Order, fill = log10(total_OTU + 1)), 
               width = 4, offset = .1, colour = "white") +
    scale_fill_viridis_c(option = "rocket", limits = nOTU_lims,
                         guide = guide_colorbar(title.position = "top",
                                                byrow = TRUE)) +
    theme_void(base_size = bs) +
    theme(legend.text = element_text(size = legend_text_size), 
          legend.title = element_text(size = legend_title_size, face = "bold"),
          legend.title.align = 0.5,
          legend.position = legend.position,
          legend.title.position = "top",
          legend.spacing.x = unit(1, "cm"),
          legend.box = "horizontal", legend.margin = margin(),
          plot.margin = margins) +
    labs(fill = "log10(nOTU)")
 
# plot & render -----------------------------------------------------------


# Save the plot as a TIFF file
 library(patchwork)
tiff("figures/malaise_taxonomy_se.tiff", width = 1500, height = 1000, compression = "lzw")
p1 
dev.off()


tiff("figures/malaise_taxonomy_mg.tiff", width = 1500, height = 1000, compression = "lzw")
p2 
dev.off()


tiff("figures/soil_taxonomy_se.tiff", width = 1500, height = 1000, compression = "lzw")
p3
dev.off()


tiff("figures/soil_taxonomy_mg.tiff", width = 1500, height = 1000, compression = "lzw")
p4
dev.off()


# Open the saved TIFF file
# browseURL("figures/malaise_taxonomy_se.tiff")
# browseURL("figures/malaise_taxonomy_mg.tiff")
# browseURL("figures/soil_taxonomy_se.tiff")
# browseURL("figures/soil_taxonomy_mg.tiff")
