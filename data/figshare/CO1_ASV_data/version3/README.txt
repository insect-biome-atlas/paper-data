## General information 

Author: Fredrik Ronquist
Contact e-mail: fredrik.ronquist@nrm.se
DOI: 10.17044/scilifelab.25480681
License: CC BY 4.0
This readme file was last updated: 2024-06-18

Please cite as: Miraldo A, Iwaszkiewicz-Eggebrecht E, Sundh J, Lokeshwaran M, Granqvist E, Andersson AF, Lukasik P, Roslin
T, Tack A, Ronquist F. 2024. Dataset of amplicon sequence variants (ASVs) from the Insect Biome Atlas Project, version
3, https://doi.org/10.17044/scilifelab.25480681

The Insect Biome Atlas project was supported by the Knut and Alice Wallenberg Foundation (dnr 2017.0088). The project
analyzed the insect faunas of Sweden and Madagascar, and their associated microbiomes, mainly using DNA metabarcoding of
Malaise trap samples collected in 2019 (Sweden) or 2019–2020 (Madagascar).

## Dataset description 

This dataset contains amplicon sequence variants (ASVs) generated from high-throughput sequencing of the cytochrome c
oxidase subunit I (COI) gene from Malaise trap samples (lysates, homogenates and preservative ethanol) and soil and
litter samples. It includes ASV sequences and abundance information (number of reads) as well as metadata files that are
needed to interpret and analyse the data further. Future versions of the dataset will include additional data. NB! All
ASV files include ASVs that represent biological and artificial spike-ins. 

## Methods

Samples were sequenced using Illumina technology. Raw data are available at the European Nucleotide Archive (ENA) under
project PRJEB61109. The raw sequence data was preprocessed using a Snakemake workflow available at
https://github.com/biodiversitydata-se/amplicon-multi-cutadapt. Preprocessed reads were used as input to the AmpliSeq
Nextflow (v.2.1.0) (https://github.com/nf-core/ampliseq) pipeline to generate ASVs.

## Available data 
Three types of files are included: ASV files, metadata files and other complementary data files. Files marked 'SE'
contain data from Sweden while those marked 'MG' contain data from Madagascar.

The file shasum.txt contains checksums for each of the files. Run:

shasum -c shasum.txt

to check file integrity.

### ASV files 
ASV sequences in fasta format are found in files CO1_asv_seqs_SE.fasta.gz and CO1_asv_seqs_MG.fasta.gz. Counts of ASVs
in each sample are in CO1_asv_counts_SE.tsv.gz and CO1_asv_counts.MG.tsv.gz. The Swedish dataset contains 821,559 ASVs
in 6,169 samples. The Madagascar dataset contains 701,769 ASVs in 2,287 samples.

### Metadata files

Four types of metadata files are included: 
1. sequencing_metadata files with information about samples that were processed in the lab and sequenced
2. samples_metadata files with information about samples that were collected in the field.
3. sites_metadata files with information about sites where samples were collected.
4. sipke-ins metadata files with information about spike-ins added to each malaise trap sample at the time of sample processing in the lab.


#### Sequencing metadata files

The two sequencing metadata files CO1_sequencing_metadata_SE.tsv and CO1_sequencing_metadata_MG.tsv contain information
about samples that were sequenced. The columns in these files are as follows:

- sampleID_NGI: Sample id given by the sequencing facility (matching the columns in the counts file)
- sampleID_HISTORICAL: Custom user id
- sampleID_FIELD: Sample id from field sampling
- sampleID_LAB: Sample id from handling in the lab
- dataset: Dataset designation for each sample
- lab_sample_type: Type of sample, e.g. 'sample', 'buffer_blank', 'pcr_neg' etc.
- country: Country of origin for sample
- biological_spikes: True if sample has biological spike ins added
- artificial_spikes: True if sample has artificial spike ins added at the time of DNA purification
- sample_metadata_file: Corresponding metadata file for sample
- lysate_rack_ID: Identification of 96-well plate where lysate aliquot is stored in the lab (internal use only)
- lysate_well_ID: Identification of well position where lysate aliquot is stored in the lab (internal use only)
- dna_plate_ID: Identification of 96-well plate where purified DNA is stored in the lab (internal use only)
- dna_plate_well_ID: Identification of well position where lysate is stored in the lab (internal use only)
- sequencing_batch: Custom user id for sequencing batch number
- sequencing_batch_NGI: Sequencing batch number given by the sequencing facility
- notes_lab: Additional information about sample processing in the lab (only for SE file)
- sequencing_successful: True if sequencing was successful for the sample. If False then this sample will be missing from the ASV counts file.
- study_accession_ENA: Study identification at the European Nucleotide Archive
- sample_accession_ENA: Sample identification at the European Nucleotide Archive
- experiment_accession_ENA: Experiment identification at the European Nucleotide Archive
- run_accession_ENA: Run identification at the European Nucleotide Archive  

#### Samples metadata files

Four samples_metadata files are included in this dataset with information about each sample that was collected in the
field. For samples collected with malaise traps we have two files, one for each country: samples_metadata_malaise_SE.tsv
and samples_metadata_malaise_MG.tsv. The columns in these files are as follows:

- sampleID_FIELD: Sample id from field sampling
- trapID: Malaise trap id from field sampling
- biomass_grams: Wet weight of each bulk sample
- placing_time: Time when sampling started
- placing_date: Date when sampling started
- collecting_time: Time when sampling ended
- collecting_date: Date when sampling ended
- duration_min: Total number of minutes the sample was collecting
- trap_condition_collection: Condition of the malaise trap at the time of collecting the sample from the trap (good; acceptable; poor)
- sample_ethanol_conc: Concentration of preservative ethanol at the time of DNA extraction (only for SE file)
- processing_group: Processing batch id (for internal use only)
- sample_accession_ENA: Sample identification at the European Nucleotide Archive
- sample_status: Additional information about sample processing status in the lab

For arthropod samples collected from litter and soil we have two files, one for each country:
samples_metadata_soil_litter_SE.tsv and samples_metadata_litter_MG.tsv. Note that for Madagascar we did not collect
arthropod samples from soil. Also note that for Madagascar we collected four leaf litter samples at each trap location,
one sample in each direction of the Malaise trap (front, back, left and right); whilst for Sweden we collected only one
sample at each trap location. Columns in these files are as follows:

- sampleID_FIELD: Sample id from field sampling.
- trapID: Malaise trap id from field sampling.
- sample_type: Identifies if a sample is a soil or litter sample (only for SE file).
- transectID: Identifies the transect where the sample was collected: front of Malaise trap (transectID=1), right hand side of Malaise trap (transectID=2), back of Malaise trap (transectID=3), left hand side of Malaise trap (transectID=4) (only for MG file).
- biomass_grams: Wet weight of each bulk sample (only for MG file).
- date: Date when sample was collected.
- time: Time when sample was collected.
- sample_accession_ENA: Sample identification at the European Nucleotide Archive.
- sample_status: Additional information about sample processing status in the lab.

#### Sites metadata files

There are two files that contain information about sampling sites, one for each country: sites_metadata_SE.tsv and
sites_metadata_MG.tsv. Columns in these files are as follows:

- siteID: Sampling site id number. Note that for some sites there can be several Malaise traps assembled (malaise_trap_type=Multitrap)
- trapID: Malaise trap id from field sampling
- latitude_WGS84: Latitude in WGS84 coordinate system. This info specifies the Malaise trap location at the sampling site
- longitude_WGS84: Longitude in WGS84 coordinate system. This info specifies the Malaise trap location at the site
- trap_habitat: Habitat where the Malaise trap was located
- malaise_trap_type: Identifies if there are multiple traps assembled at the sampling site (Multitrap) or only one (Single_trap)
- parkID: Name of national park (for MG only)
- provinceID: Name of province (for MG only)
- NILS_mhabitat: Habitat for nearest plot of the National Inventory of Landscapes in Sweden (NILS) from the malaise trap location (only for SE file). For more information about NILS sampling design, check: https://www.slu.se/centrumbildningar-och-projekt/nils_old/Datainsamling/bakgrund-och-mal/
- NILS_square: Identification of nearest NILS square for sampling site (only for SE file)
- NILS_plot: Identification of nearest NILS plot to the Malaise trap location (only for SE file)
- trap_orientation_degrees_S: Orientation in degrees of the collection head of the Malaise trap
- notes: notes associated with the Malaise trap (only for SE file)

#### Spike-ins metadata files

We provide three files with information about spike-ins used when processing samples in the lab: biological_spikes_taxonomy_SE.tsv and biological_spikes_taxonomy_MG.tsv contain taxonomic information on biological spike ins while the file synthetic_spikes_info.tsv has information on synthetic spike ins. See README.txt for more information.

Columns in the files biological_spikes_taxonomy_[SE|MG].tsv are as follows:

- Kingdom, Phylum, Class, Order, Family, Genus, Species: Taxonomic assignment for each spike-in.
- n: number of specimens added to each malaise trap sample.
- Comment: comments on alternative species nomenclature.

Columns in the file synthetic_spikes_info.tsv are as follows:

- Cluster_ID: identification synthetic spike-in cluster.
- Spike_name: name of synthetic spike-in.
- Amount_added: number of copies added (M=million).
- Sequence: nucleotide sequence.

### Other complementary data files

We present complementary data on soil chemistry collected at each sampling location in both Sweden and Madagascar, stand
characteristics collected at each sampling location in Madagascar and biomass/count data for a selected number of
malaise trap samples from the Insect Biome Atlas project (n=24) and the Swedish Insect Inventory Project (n=224).

#### Soil chemistry data

We provide two datasets, one for each country, on soil chemistry (soil_chemistry_SE.tsv and soil_chemistry_MG.tsv) that
store information on soil nutrients from soil samples collected at the same sampling sites as the arthropod communities.
Topsoil (0-20cm) was sampled at 5 sites around each Malaise trap in both Sweden and Madagascar: one soil core (6 cm
diameter) at the center of trap and one soil core on each of the four “sides” of the trap five meters away from the
trap. Soil samples at each site were taken as composite samples from the five locations. Soil samples collected in
Sweden were analysed at Eurofins in Sweden and the ones collected in Madagascar were analysed at the Laboratoire des
Radioisotopes in Madagascar. As samples from each country were analysed at different laboratories the variables on soil
nutrients presented in each dataset differ slightly. 

soil_chemistry_SE.tsv
- ​trapID: Malaise trap id from field sampling that identifies the location where soil samples for chemistry analyses were taken.
- date: date (yyyy-mm-dd) when soil sample was collected. 
- Raining: true if it was raining when the soil sample was taken; false if it was not raining when the soil sample was taken. 
- litter_depth_center (mm): litter depth in millimeters measured at the center of the Malaise trap.
- litter_depth_1 (mm): litter depth in millimeters measured five meters away from the collecting head of the Malaise trap (front of Malaise trap). 
- litter_depth_2 (mm): litter depth in millimeters measured five meters away from the right hand side of the Malaise trap at 90 degrees angle from the collecting head of the Malaise trap.
- litter_depth_3 (mm): litter depth in millimeters measured five meters away from the back of the Malaise trap.
- litter_depth_4 (mm): litter depth in millimeters measured five meters away from the right hand side of the Malaise trap at 270 degrees angle from the collecting head of the Malaise trap.
- soil_humidity_center (%): soil humidity in percentage measured at the center of the Malaise trap.
- soil_humidity_1 (%): soil humidity in percentage measured five meters away from the collecting head of the Malaise trap (front of Malaise trap). 
- soil_humidity_2 (%): soil humidity in percentage measured five meters away from the right hand side of the Malaise trap at 90 degrees angle from the collecting head of the Malaise trap.
- soil_humidity_3 (%): soil humidity in percentage measured five meters away from the back of the Malaise trap.
- soil_humidity_4 (%): soil humidity in percentage measured five meters away from the right hand side of the Malaise trap at 270 degrees angle from the collecting head of the Malaise trap.
- pH: pH in water of the composite soil sample. ISO:10390: 2005-12.
- NH4-N (mg/100g): nitrogen content of ammonium ion present in the composite soil sample. Measured using ADAS method 53.
- NO3-N (mg/100g): nitrogen content of ammonia present in composite soil sample. ADAS method 53.
- P_AL (mg/100g): phosphorus present in the composite soil sample measured using Amonium lactate method (AL method as in Egnér et al. 1960). ISO 118855-2009:09.
- K_AL (mg/100g): potassium present in the composite soil sample measured using Amonium lactate method (AL method as in Egnér et al. 1960). ISO 118855-2009:09.
- Mg_AL (mg/100g): magnesium present in the composite soil sample measured using Amonium lactate method (AL method as in Egnér et al. 1960). ISO 118855-2009:09.
- Ca_AL (mg/100g): calcium present in the composite soil sample measured using Amonium lactate method (AL method as in Egnér et al. 1960). ISO 118855-2009:09.
- K/Mg_ratio: potassium to magnesium ratio.
- Nmin: minimum total nitrogen.

soil_chemistry_MG.tsv
- trapID: Malaise trap id from field sampling that identifies the location where soil samples for chemistry analyses were taken.
- date: date (yyyy-mm-dd) when soil sample was collected. 
- litter_depth_center (mm): litter depth in millimeters measured at the center of the Malaise trap.
- litter_depth_1 (mm): litter depth in millimeters measured five meters away from the collecting head of the Malaise trap (front of Malaise trap). 
- litter_depth_2 (mm): litter depth in millimeters measured five meters away from the right hand side of the Malaise trap at 90 degrees angle from the collecting head of the Malaise trap.
- litter_depth_3 (mm): litter depth in millimeters measured five meters away from the back of the Malaise trap.
- litter_depth_4 (mm): litter depth in millimeters measured five meters away from the right hand side of the Malaise trap at 270 degrees angle from the collecting head of the Malaise trap.
- soil_humidity_center (%): soil humidity in percentage measured at the center of the Malaise trap.
- soil_humidity_1 (%): soil humidity in percentage measured five meters away from the collecting head of the Malaise trap (front of Malaise trap). 
- soil_humidity_2 (%): soil humidity in percentage measured five meters away from the right hand side of the Malaise trap at 90 degrees angle from the collecting head of the Malaise trap.
- soil_humidity_3 (%): soil humidity in percentage measured five meters away from the back of the Malaise trap.
- soil_humidity_4 (%): soil humidity in percentage measured five meters away from the right hand side of the Malaise trap at 270 degrees angle from the collecting head of the Malaise trap.
- pH: pH in water of the composite soil sample.
- organic_carbon (g/kg): organic carbon from the composite soil sample.
- total_nitrogen (g/kg): total nitrogen from the composite soil sample.
- total_phosphorus (g/kg): total phosphorus from the composite soil sample.
- potassium_exchange (cmol/kg): exchangeable potassium from the composite soil sample.
- calcium_exchange (cmol/kg): exchangeable calcium from the composite soil sample.
- magnesium_exchange (cmol/kg): exchangeable magnesium from the composite soil sample.

#### Stand characteristics data

Standing characteristics were only measured in Madagascar as extensive data on landscape composition and vegetation
structure at the sampling sites in Sweden had already been compiled as part of the National Inventory of Landscapes in
Sweden (NILS) and data are publicly available at
https://www.slu.se/en/Collaborative-Centres-and-Projects/nils/nils-inventory-2003-2020/. 

The file stand_characteristics_MG.tsv contains information on a set of standing characteristics from Madagascar related
to tree density (DBH, shading, etc). Columns in this file are as follows:

- trapID: Malaise trap id from field sampling that identifies the location where standing characteristics were measured.
- ​​canopy_cover_center (%): percentage of canopy cover measured at the center of the malaise trap.
- canopy_cover_1 (%): percentage of canopy cover measured five meters away from the collecting head of the Malaise trap (front of Malaise trap).
- canopy_cover_2 (%): percentage of canopy cover measured to the left hand side of the Malaise trap (when one is facing the collecting head), five meters away from the trap. 
- canopy_cover_3 (%): percentage of canopy cover measured five meters away from the back of the Malaise trap.
- canopy_cover_4 (%): percentage of canopy cover measured to the right hand side of the Malaise trap (when one is facing the collecting head), five meters away from the trap. 
- T1_DBH_1 to T1_DBH_n: diameter at breast height of tree 1 (tree 2, tree 3, tree n) that have its roots or stem within 10 cm of transect 1 (T1). Transect 1 is 10 m long starting at the collecting head of the Malaise trap. (L) denotes lianas, (P) denotes palms and (F) denotes ferns. Units in cm.
- T2_DBH_1 to T2_DBH_n: diameter at breast height of tree 1 (tree 2, tree 3, tree n) that have its roots or stem within 10 cm of transect 2 (T2). Transect 2 is 10 m long starting at the extremity of the Malaise trap at a direction of 90 degrees angle from transect 1. (L) denotes lianas, (P) denotes palms and (F) denotes ferns. Units in cm.
- T3_DBH_1 to T3_DBH_n: diameter at breast height of tree 1 (tree 2, tree 3, tree n) that have its roots or stem within 10 cm of transect 3 (T3). Transect 3 is 10 m long starting at the back of the Malaise trap. (L) denotes lianas, (P) denotes palms and (F) denotes ferns. Units in cm.
- T4_DBH_1 to T4_DBH_n: diameter at breast height of tree 1 (tree 2, tree 3, tree n) that have its roots or stem within 10 cm of transect 4 (T4). Transect 4 is 10 m long starting at the extremity of the Malaise trap at a direction of 270 degrees angle from transect 1. (L) denotes lianas, (P) denotes palms and (F) denotes ferns. Units in cm.
- BIG TREE_1 to BIG_TREE_n: circumference of each tree with a DBH larger than 30 cm and located within a 10 m radius of the Malaise trap. Units in cm.

#### Biomass and count data

To allow an assessment of how the biomass of a Malaise trap sample translates to the number of specimens, we provide two
files describing samples from Sweden, for which we measured the biomass and also counted all the specimens in the
sample. The first set comprises 24 samples from the IBA field campaign (biomass_count_IBA.tsv), and the second set
comprises 224 samples from a separate Swedish Malaise trapping campaign (Swedish Insect Inventory Project) in 2018–2019
(biomass_count_SIIP.tsv). For the latter dataset, we provide the site and sample metadata in the same file.

biomass_count_IBA.tsv
- sampleID_FIELD: Sample id from Insect Biome Atlas field sampling campaign.
- trapID: Malaise trap id from field sampling that identifies the location where sample was collected.
- biomass_grams: Wet weight of each bulk sample biomass in grams.
- count: Total number of specimens present in the bulk sample.

biomass_count_SIIP.tsv
- siteID: Sampling site id number.
- sampleID: Sample id from the Swedish Insect Inventory project field sampling campaign.
- start_date: Date (yyyy-mm-dd) when sampling started.
- end_date: Date (yyyy-mm-dd) when sampling ended.
- region: Region where sample was collected.
- province: Province where samples was collected.
- municipality: Municipality where sample was collected.
- locality: Locality where sample was collected.
- latitude_WGS84: Latitude in WGS84 coordinate system.
- longitude_WGS84: Longitude in WGS84 coordinate system.
- habitat: Habitat where the Malaise trap was located.
- biomass_grams: Wet weight of each bulk sample biomass in grams.
- count: Total number of specimens present in the bulk sample.

References:
Egnér, H., Riehm, H., & Domingo, W. (1960). Untersuchungen über die chemische Bodenanalyse als Grundlage für die
Beurteilung des Nährstoffzustandes der Böden. II. Chemische Extraktionsmethoden zur Phosphor-und Kaliumbestimmung.
Kungliga Lantbrukshögskolans Annaler, 26, 199–215.
