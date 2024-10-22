# Processed ASV data from the Insect Biome Atlas Project

## General information

Author: Fredrik Ronquist
Contact e-mail: fredrik.ronquist@nrm.se
DOI: 10.17044/scilifelab.27202368
License: CC BY 4.0
This readme file was last updated: 2024-10-12

The Insect Biome Atlas project was supported by the Knut and Alice Wallenberg Foundation (dnr 2017.0088). The project
analyzed the insect faunas of Sweden and Madagascar, and their associated microbiomes, mainly using DNA metabarcoding of
Malaise trap samples collected in 2019 (Sweden) or 2019–2020 (Madagascar).

Please cite this version of the dataset as:Miraldo A, Iwaszkiewicz-Eggebrecht E, Sundh J, Lokeshwaran M, Granqvist E,
Goodsell R, Andersson AF, Lukasik P, Roslin T, Tack A, Ronquist F. 2024. Processed ASV data from the Insect Biome Atlas
Project, version 1. doi:10.17044/scilifelab.27202368.v1 or https://doi.org/10.17044/scilifelab.27202368.v1

## Dataset description

This dataset contains the results from bioinformatic processing of version 1 of the amplicon sequence variant (ASV) data
from the Insect Biome Atlas project (Miraldo et al. 2024), that is, the cytochrome oxidase subunit 1 (CO1)
metabarcoding data from Malaise trap samples processed using the FAVIS mild lysis protocol (Iwaszkiewicz et al. 2023).
The bioinformatic processing involved: (1) taxonomic assignment of ASVs, (2) chimera removal; (3) clustering into OTUs;;
and (4) noise filtering. The clustering step involved resolution of the taxonomic annotation of the cluster and
identification of a representative ASV. The noise filtering step involved removal of ASV clusters identified as
potentially originating from nuclear mitochondrial DNA (NUMTs) or representing other types of error or noise. ASV
taxonomic assignments, ASV cluster designations, consensus taxonomies and summed counts of clusters in the sequenced
samples are provided in compressed tab-separated files. Sequences of cluster representatives are provided in compressed
FASTA format files. The bioinformatic processing pipeline is further described in Sundh et al. (2024). NB! All result
files include ASVs and clusters that represent biological spike-ins.

## Methods

### Taxonomic assignment 

ASVs were taxonomically assigned using kmer-based methods implemented in a Snakemake workflow available at
https://github.com/insect-biome-atlas/ASV-taxonomy. Specifically ASVs were assigned a taxonomy using the SINTAX
algorithm in vsearch (v2.21.2) using a CO1 database constructed from the Barcode Of Life Data System (Sundh 2022).

### Chimera removal

The workflow first identifies chimeric ASVs in the input data using the ‘uchime_denovo’ method implemented in vsearch.
This was done with a so-called ‘strict samplewise’ strategy where each sample was analysed separately (hence the
‘samplewise’ notation), only comparing ASVs present in the same sample. Further, ASVs had to be identified as chimeric
in all samples where they were present (corresponding to the ‘strict’ notation) in order to be removed as chimeric. 

### ASV clustering

Non-chimeric sequences were then split by family-level taxonomic assignments and ASVs within each family were clustered
in parallel using swarm (v3.1.0) with differences=13. Representative ASVs were selected for each generated cluster by
taking the ASV with the highest relative abundance across all samples in a cluster. Counts were generated at the cluster
level by summing over all ASVs in each cluster.

### Consensus taxonomy

A consensus taxonomy was created for each cluster by taking into account the taxonomic assignments of all ASVs in a
cluster as well as the total abundance of ASVs. For each cluster, starting at the most resolved taxonomic level, each
unique taxonomic assignment was weighted by the sum of read counts of ASVs with that assignment. If a single weighted
assignment made up 80% or more of all weighted assignments at that rank, that taxonomy was propagated to the ASV
cluster, including parent rank assignments. If no taxonomic assignment was above the 80% threshold, the algorithm
continued to the parent rank in the taxonomy. Taxonomic assignments at any available child ranks were set to the
consensus assignment prefixed with ‘unresolved’.

### Noise filtering and cleaning

ASV clusters were further analyzed to remove clusters likely to originate from nuclear mitochondrial DNA (NUMTs) or to
represent other types of errors or noise, using an algorithm based on abundance, negative controls and taxonomic
annotation uncertainty. Specifically, we removed clusters for which we only had one or two reads. For Swedish data, we
also removed clusters that were not successfully assigned to at least family level by SINTAX; according to several
tests, this procedure removed potential NUMTs with good accuracy (Sundh et al. 2024). For the Madagascar data, the
SINTAX annotation often failed at the family level, indicating that the taxonomy annotation filtering flagged too many
clusters as potential NUMTs. We therefore filtered the Madagascar data only based on the number of reads.

As a last clean-up step in the noise filtering, clusters containing at least one ASV present in more than 5% of blanks
were removed.

The chimera filtering and ASV clustering methods have been implemented in a Snakemake workflow available at
https://github.com/insect-biome-atlas/happ. This workflow takes as input:

1. The ASV sequences in FASTA format
2. A tab-delimited file of counts of ASVs (rows) in samples (columns)
3. Taxonomic assignments of ASVs

Data for 1) and 2) are available at https://doi.org/10.17044/scilifelab.25480681.v1, and 3) is available in this upload.

Noise filtering and cleaning was done with custom python scripts available at
https://github.com/insect-biome-atlas/utils-clean_asv_data. 

## Available data

### Processed ASV data files

ASV taxonomic assignments, non-chimeric ASV cluster designations, consensus taxonomies, sequences of cluster
representatives and summed counts of clusters in the sequenced samples are provided in compressed tab-separated files.
Files are organized by country (Sweden and Madagascar), marked by the suffixes SE and MG, respectively.

#### Taxonomic assignments

The files asv_taxonomy_[SE|MG].tsv.gz are tab-separated files with taxonomic assignments for all ASVs. Columns:

- ASV: The id of the ASV
- Kingdom, Phylum, Class, Order, Family, Genus, Species, BOLD_bin: Taxonomic assignment for each rank. 

If an ASV was unclassified at a particular rank, the taxonomic label is prefixed with ‘unclassified.’ followed by the
taxonomic assignment of the most resolved parent rank.

#### Cluster assignments

The files cluster_taxonomy_[SE|MG].tsv are tab-separated files containing all non-chimeric ASVs (that is, the ASVs
passing the chimera-filtering step) with their corresponding taxonomic and cluster assignments. Columns:

- ASV: ASV id
- cluster: name of designated cluster
- median: the median of normalized reads across all samples for each ASV
- Kingdom, Phylum, Class, Order, Family, Genus, Species, BOLD_bin: taxonomic assignment of each ASV
- representative: contains 1 if ASV is a representative of its cluster, otherwise 0

#### Cluster counts

The files cluster_counts_[SE|MG].tsv are tab-separated files with read counts of ASV clusters (rows) in samples
(columns). Counts have been summed for all ASVs belonging to each cluster.

#### Sequences of cluster representatives

The files cluster_reps_[SE|MG].fasta are text files in FASTA format with representative sequences for each cluster. The
fasta headers have the format “>ASV_ID CLUSTER_NAME”.

#### Consensus taxonomy

The files cluster_consensus_taxonomy_[SE|MG].tsv are tab-separated files with consensus taxonomy of each generated ASV
cluster. Columns are the same as in asv_taxonomy_[SE|MG].tsv.
