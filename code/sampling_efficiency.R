library(tidyverse)
library(data.table)
library(vegan)
library(patchwork)

set.seed(10)
# ---------------------------------------------------------------------------------------------
# data ----------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# Spike in IDs 
spike_ins_mg  <- fread("data/Biological_spikes_MG_taxonomy.tsv")
spike_ins_se  <- fread("data/Biological_spikes_SE_taxonomy.tsv")

# --------------- Get sample IDs 
seq_meta_se <- fread("data/co1_sequencing_metadata_se.tsv")
seq_meta_mg <- fread("data/co1_sequencing_metadata_mg.tsv")

# Sample meta
trap_meta_se <- fread("data/samples_metadata_malaise_SE_2019.tsv")
trap_meta_mg <- fread("data/samples_metadata_malaise_MG_2019.tsv")

# Site meta
site_meta_se <- fread("data/sites_metadata_SE_2019.tsv")
site_meta_mg <- fread("data/sites_metadata_MG_2019.tsv")


# ----------------- Get and clean cluster data
# These files need downloading from figshare

swe_clusters <- fread("data/non_numts_cleaned_SE.tsv") |>
  _[representative == 1 ,] |> 
  _[Phylum == "Arthropoda" ,] |> 
  _[!str_detect(Family , "_X*"),] |> 
  _[!str_detect(Species , paste(spike_ins_se$Species,collapse="|")),]

mg_clusters  <- fread("data/non_numts_cleaned_SE.tsv") |>
  _[representative == 1 ,] |> 
  _[Phylum == "Arthropoda" ,] |> 
  _[!str_detect(Family , "_X*"),] |> 
  _[!str_detect(Species , paste(spike_ins_mg$Species,collapse="|")),]

# ------------ Get read numbers 
swe_counts   <- fread("data/non_numts_cluster_counts_cleaned_SE.tsv")
mg_counts    <- fread("data/non_numts_cluster_counts_cleaned_MG.tsv")


# plotting pars -------------------------------------------------------------------------------
bs    <- 35
lsize <- 2

# ---------------------------------------------------------------------------------------------
# tidy data ----------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------


# First get counts of each species in each sample
swe_counts_long <- swe_counts |> 
                   melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
                   _[read_count > 5,] 

# Merge with cluster IDS
taxa_counts_DT <- swe_counts_long[swe_clusters, on = "cluster"] 


# Assemble the metadata
full_meta_DT <-  merge(seq_meta_se , trap_meta_se , all=TRUE) 
full_meta_DT <- merge(full_meta_DT , site_meta_se,by="trapID") |> _[lab_sample_type %in% "sample",]

# Combine metadata and OTU data
OTU_DT <- merge(taxa_counts_DT , full_meta_DT, by = "sampleID_NGI") |> 
          _[lab_sample_type == "sample", ] 

# ---------------------------------------------------------------------------------------------
# madagascar data -----------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------


# First get counts of each species in each sample
mg_counts_long <- mg_counts |> 
  melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
  _[read_count > 5,] 

# Merge with cluster IDS
taxa_counts_DT_mg <- mg_counts_long[mg_clusters, on = "cluster"] 

# Assemble the metadata
full_meta_DT_mg <-  merge(seq_meta_mg , trap_meta_mg , all=TRUE) 
full_meta_DT_mg <- merge(full_meta_DT_mg , site_meta_mg,by="trapID") |> _[lab_sample_type %in% "sample",]

# Combine metadata and OTU data
OTU_DT_mg <- merge(taxa_counts_DT_mg , full_meta_DT_mg, by = "sampleID_NGI") |> 
          _[lab_sample_type == "sample", ] 


# ---------------------------------------------------------------------------------------------
# select indices of 20 random samples ---------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# Traps
site_ind <- sample(1:max(site_meta_se$siteID) , 10)
traps_se_10 <- filter(site_meta_se , siteID %in% site_ind) |> filter(trapID == first(trapID) , .by = siteID) |> pull(trapID)

site_ind_mg <- sample(1:max(site_meta_mg$siteID) , 10)
traps_mg_10 <- filter(site_meta_mg , siteID %in% site_ind_mg) |> filter(trapID == first(trapID) , .by = siteID) |> pull(trapID)

# Samples
samples_se_10 <- full_meta_DT |> filter(trapID %in% traps_se_10) |> group_by(trapID) |> sample_n(1) |> pull(sampleID_NGI)
samples_mg_10 <- full_meta_DT_mg |> filter(trapID %in% traps_mg_10) |> group_by(trapID) |> sample_n(1) |> pull(sampleID_NGI)

# ---------------------------------------------------------------------------------------------
# Sequence depth accumulation curve -----------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# Get table of 10 random samples by cluster
sample_tx_se <- OTU_DT |>
            _[ sampleID_NGI %in% samples_se_10,] |> 
            _[,c("sampleID_NGI" , "cluster" , "read_count")] |> 
            _[,lapply(.SD, sum), by = c("sampleID_NGI" , "cluster"), .SDcols = "read_count"] |> 
            dcast(sampleID_NGI ~ cluster  ,fun.aggregate = sum, value.var = "read_count") |> 
            drop_na(sampleID_NGI) |> 
            column_to_rownames("sampleID_NGI")


sample_tx_mg <- OTU_DT_mg |>
                _[ sampleID_NGI %in% samples_mg_10,] |> 
                _[,c("sampleID_NGI" , "cluster" , "read_count")] |> 
                _[,lapply(.SD, sum), by = c("sampleID_NGI" , "cluster"), .SDcols = "read_count"] |> 
                dcast(sampleID_NGI ~ cluster  ,fun.aggregate = sum, value.var = "read_count") |> 
                drop_na(sampleID_NGI) |> 
                column_to_rownames("sampleID_NGI")

#Rarefy ------------------------------------

swe_rare_samp <- rarecurve(sample_tx_se , step = 10, tidy=TRUE)
mdg_rare_samp <- rarecurve(sample_tx_mg , step = 10, tidy=TRUE)

# Sweden ------------------------------------


p1 <- ggplot(swe_rare_samp , aes(Sample , Species))+
  geom_line(size = lsize , aes(group = Site))+
  theme_linedraw(base_size=bs)+
  labs(x = "Number of sequences" , y = "Number of clsuters")+
  scale_y_continuous(limits = c(0,600))


# Madagascar ------------------------------------
  
p2 <- ggplot(mdg_rare_samp , aes(Sample , Species))+
      geom_line(size = lsize , aes(group = Site))+
       theme_linedraw(base_size=bs)+
        labs(x = "Number of sequences" , y = "Number of clusters")+
       scale_y_continuous(limits = c(0,600))


# ---------------------------------------------------------------------------------------------
# Sampling intensity --------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------


# Traps
traps_se_all <- filter(site_meta_se , trapID == first(trapID) , .by = siteID) |> pull(trapID)
traps_mg_all <- filter(site_meta_mg , trapID == first(trapID) , .by = siteID) |> pull(trapID)


#Trap
trap_tx_se <- OTU_DT |>
            _[trapID %in% traps_se_all , ] |> 
            _[,c("trapID" , "cluster" , "read_count")] |> 
            _[,lapply(.SD, function(x) (sum(x)>1)*1), by = c("trapID" , "cluster"), .SDcols = "read_count"] |> 
            dcast(trapID ~ cluster  ,fun.aggregate = sum, value.var = "read_count") |> 
            drop_na(trapID) |> 
            column_to_rownames("trapID")
  

trap_tx_mg <- OTU_DT_mg |>
              _[trapID %in% traps_mg_all , ] |> 
              _[,c("trapID" , "cluster" , "read_count")] |> 
              _[,lapply(.SD, function(x) (sum(x)>1)*1), by = c("trapID" , "cluster"), .SDcols = "read_count"] |> 
              dcast(trapID ~ cluster  ,fun.aggregate = sum, value.var = "read_count") |> 
              drop_na(trapID) |> 
              column_to_rownames("trapID")



# Species accumulation curves
spSE <- specaccum(trap_tx_se , method = "random" , permutations = 100)
spMG <- specaccum(trap_tx_mg , method = "random" , permutations = 100)

seDF <- data.frame(species = spSE$richness , sd  = spSE$sd , sites = spSE$sites)
mgDF <- data.frame(species = spMG$richness , sd  = spMG$sd , sites = spMG$sites)

# Sweden
p3 <- ggplot(seDF , aes(sites  , species))+
        geom_line(size = lsize)+
        geom_ribbon(aes(ymin = species - sd , ymax = species + sd) , alpha = .2 , fill= "blue")+
        theme_linedraw(base_size=bs)+
        labs(x = "Number of sites" , y = "Number of clusters")
  


#Madagascar 
p4 <- ggplot(mgDF , aes(sites  , species))+
        geom_line(size = lsize)+
        geom_ribbon(aes(ymin = species - sd , ymax = species + sd) , alpha = .2 , fill= "blue")+
        theme_linedraw(base_size=bs)+
        labs(x = "Number of sites" , y = "Number of clusters")


# ---------------------------------------------------------------------------------------------
# Within site sampling efficiency -------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# Traps
traps_se_multi <- filter(site_meta_se , trap_type == "Multitrap") |> mutate(n_traps = n() , .by = siteID) 
traps_mg_multi <- filter(site_meta_mg , malaise_trap_type == "Multitrap") |> mutate(n_traps = n() , .by = siteID)

# Iterate over multitrap sites for both countries. 

# Sweden ------------------------------------ 
se_list <- list()
for(i in 1:max(traps_se_multi$n_traps)){

    # Get index of traps at each site 
    se_traps_l <- traps_se_multi  |>
                  group_split(siteID) |>
                  lapply(X=_, FUN=function(x) {x %>% sample_n(ifelse(nrow(x) < i, nrow(x), i))}) |> 
                  bind_rows() |> 
                  pull(trapID)
    
    # Get species table
    trap_tx_se_l <- OTU_DT |>
                  _[trapID %in% se_traps_l , ] |> 
                  _[,c("trapID" , "cluster" , "read_count")] |> 
                  _[,lapply(.SD, function(x) (sum(x)>1)*1), by = c("trapID" , "cluster"), .SDcols = "read_count"] |> 
                  dcast(trapID ~ cluster  ,fun.aggregate = sum, value.var = "read_count") |> 
                  drop_na(trapID) |> 
                  column_to_rownames("trapID")
    
    # Species accumulation & save
    spSE         <- specaccum(trap_tx_se_l , method = "random" , permutations = 100)
    se_list[[i]] <- data.frame(species = spSE$richness , sd  = spSE$sd , sites = spSE$sites , ntraps = i)

}


# Madagascar ------------------------------------ 
mg_list <- list()
for(i in 1:max(traps_mg_multi$n_traps)){
  
    # Get index of traps at each site 
    mg_traps_l <- traps_mg_multi  |>
                    group_split(siteID) |>
                    lapply(X=_, FUN=function(x) {x %>% sample_n(ifelse(nrow(x) < i, nrow(x), i))}) |> 
                    bind_rows() |> 
                    pull(trapID)
    
    # Get species table
    trap_tx_mg_l <- OTU_DT_mg |>
                    _[trapID %in% mg_traps_l , ] |> 
                    _[,c("trapID" , "cluster" , "read_count")] |> 
                    _[,lapply(.SD, function(x) (sum(x)>1)*1), by = c("trapID" , "cluster"), .SDcols = "read_count"] |> 
                    dcast(trapID ~ cluster  ,fun.aggregate = sum, value.var = "read_count") |> 
                    drop_na(trapID) |> 
                    column_to_rownames("trapID")
    
    # Species accumulation & save
    spMG         <- specaccum(trap_tx_mg_l , method = "random" , permutations = 100)
    mg_list[[i]] <- data.frame(species = spMG$richness , sd  = spMG$sd , sites = spMG$sites , ntraps = i)
  
}


p5 <- bind_rows(se_list) |> 
  ggplot(aes(sites  , species , group=ntraps , colour = ntraps))+
  geom_line(size = lsize)+
  theme_linedraw(base_size=bs)+
  labs(x = "Number of sites" , y = "Number of clusters" , colour = "Traps per site")+
  theme(legend.position = c(0.15,.7), legend.background = element_rect(fill=NA) , legend.key.size = unit(1,"cm"))

p6 <- bind_rows(mg_list) |> 
  ggplot(aes(sites  , species , group=ntraps , colour = ntraps))+
  geom_line(size = lsize , show.legend = FALSE)+
  theme_linedraw(base_size=bs)+
  labs(x = "Number of sites" , y = "Number of clusters" , colour = "Traps per site")


# plot all together ---------------------------------------------------------------------------

tiff("figures/figure4.tiff" , width = 2000 , height = 2000 , compression = "lzw")
(p1+p2) / (p5+p6) / (p3 + p4) + plot_annotation(tag_levels = "A")
dev.off()

