raw_data_path <- "~/dev/figshare-repos/iba/raw_data/"
processed_data_se_path <- "~/dev/figshare-repos/iba/processed_data/SE.v2/"
processed_data_mg_path <- "~/dev/figshare-repos/iba/processed_data/MG.v2/"
utils_path <- "~/dev/ms-repos-iba/utils/"
neeat_res_path <- "~/dev/ms-project-repos/develop-noise-filtering/results/"

source(paste0(utils_path,"spikes_controls_fxns.R"))


# Function for getting read stats for control samples
get_read_stats <- function(C, F, T, blanks, spikein_samples) {
    
    idx <- which(colnames(C) %in% blanks)
    idx <- c(1,idx)
    control_counts <- C[,..idx]
    include_rows <- rowSums(control_counts[,-1])!=0
    control_counts <- control_counts[include_rows,]
    control_reads <- colSums(control_counts[,-1])

    cat("Summary of control sample reads:\n")
    print(summary(control_reads))
    cat("90% quantile of control sample reads:\n")
    print(quantile(control_reads,0.90))
    cat("There are",length(blanks),"control samples in total\n")
    cat("There are",length(control_reads),"samples where sequencing was successful\n")
    cat("There are",sum(control_reads>0),"samples with reads\n")
    cat("Of these:\n")
    x <- sum(control_reads<1000 & control_reads!=0)
    n <- sum(control_reads!=0)
    cat("There are",x,"samples with <1000 reads (", round(100*x/n,2), "%)\n")
    x <- sum(control_reads<10000 & control_reads!=0)
    cat("There are",x,"samples with <10000 reads (", round(100*x/n,2), "%)\n")
    x <- sum(control_reads<100000 & control_reads!=0)
    cat("There are",x,"samples with <100000 reads (", round(100*x/n,2), "%)\n")

    spikein_clusters <- identify_spikes(C, spikein_samples, T)
    include_rows <- control_counts$cluster %in% F$cluster & !(control_counts$cluster %in% spikein_clusters)
    control_counts_cleaned <- control_counts[include_rows,]
    control_reads_cleaned <- colSums(control_counts_cleaned[,-1])

    cat("Summary of control sample reads after filtering and removal of spikeins\n")
    print(summary(control_reads_cleaned))
    cat("90% quantile of control sample reads after filtering and removal of spikeins:\n")
    print(quantile(control_reads_cleaned,0.90))
    cat("There are",sum(control_reads_cleaned>0),"samples with reads\n")
}


# ------------
# Swedish data
# ------------

C <- fread(paste0(processed_data_se_path,"cluster_counts.tsv"))
T <- read.delim(paste0(processed_data_se_path,"cluster_taxonomy.tsv"))
M <- read.delim(paste0(raw_data_path,"CO1_sequencing_metadata_SE.tsv"))
F <- read.delim(paste0(neeat_res_path,"cleaned_cluster_taxonomy_SE.tsv"))

T <- T[T$representative==1,]

blanks <- M$sampleID_NGI[M$lab_sample_type %in% c("buffer_blank","buffer_blank_art_spikes","pcr_neg","extraction_neg")]
samples <- M$sampleID_NGI[M$lab_sample_type=="sample" & M$sequencing_status=="sequencing successful"]
spikein_samples <- M$sampleID_NGI[(M$lab_sample_type=="sample") & (M$sequencing_status=="sequencing successful") & (M$dataset %in% c("CO1_lysate_2019_SE","CO1_homogenate_2019_SE"))]

# Get read stats
get_read_stats(C, F, T, blanks, spikein_samples)

# Identify control clusters
remove_clusters <- identify_control_clusters(C, T, samples, blanks)


# ---------------
# Madagascar data
# ---------------

C <- fread(paste0(processed_data_mg_path,"cluster_counts.corrected.NGI_IDs.tsv"))
T <- read.delim(paste0(processed_data_mg_path,"cluster_taxonomy.tsv"))
M <- read.delim(paste0(raw_data_path,"CO1_sequencing_metadata_MG.new.tsv"))
F <- read.delim(paste0(neeat_res_path,"cleaned_cluster_taxonomy_MG.tsv"))

T <- T[T$representative==1,]

blanks <- M$sampleID_NGI[M$lab_sample_type %in% c("buffer_blank","pcr_neg","extraction_neg")]
samples <- M$sampleID_NGI[M$lab_sample_type=="sample" & grepl("successful",M$sequencing_status)]
spikein_samples <- M$sampleID_NGI[M$lab_sample_type=="sample" & grepl("successful",M$sequencing_status) & M$dataset=="CO1_lysate_2019_MG"]

# Get read stats
get_read_stats(C, F, T, blanks, spikein_samples)

# Identify control clusters
remove_clusters <- identify_control_clusters(C, T, samples, blanks)


