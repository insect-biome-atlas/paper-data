# Set data path
path <- "~/dev/figshare-repos/iba/raw_data/"

# Read in samples metadata
samples_se <- read.delim(paste0(path,"samples_metadata_malaise_SE.tsv"))
samples_mg <- read.delim(paste0(path,"samples_metadata_malaise_MG.tsv"))

# Remove NAs and zeros for biomass
samples_se <- samples_se[!is.na(samples_se$biomass_grams) & samples_se$biomass_grams!=0,]
samples_mg <- samples_mg[!is.na(samples_mg$biomass_grams) & samples_mg$biomass_grams!=0,]

# Compute proportional estimates from lin-lin regression equation
# y = 338.5 x
prop_spec_se <- sum(samples_se$biomass) * 338.5
prop_spec_mg <- sum(samples_mg$biomass) * 338.5

# Compute non-proportional estimates from log-log regression equation
# equation unknown...

# nonprop_spec_se <- sum(10^(k * log10(samples_se$biomass) + b))
# nonprop_spec_mg <- sum(10^(k * log10(samples_mg$biomass) + b))

# Print results

cat("Specimen estimates for the proportional (lin-lin) model:\n")
cat("Sweden:", prop_spec_se, " specimens\n")
cat("Madagascar:", prop_spec_mg, " specimens\n")

# cat("Specimen estimates for the nonproportional (log-log) model:\n")
# cat("Sweden:", nonprop_spec_se, " specimens\n")
# cat("Madagascar:", nonprop_spec_mg, " specimens\n")

