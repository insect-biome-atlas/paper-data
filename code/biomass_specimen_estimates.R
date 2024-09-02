# Set data path
path <- "~/dev/figshare-repos/iba/raw_data/"

# Compute regression models
biomass_iba <- read.delim("../biomass_count_IBA.tsv", dec=",")
biomass_siip <- read.delim("../biomass_count_SIIP.tsv", dec=",")

count <- c(biomass_iba$count, biomass_siip$count)
biomass <- c(biomass_iba$biomass_grams, biomass_siip$biomass_grams)

fit_lin <- lm(count~biomass+0)
fit_log <- lm(log10(count)~log10(biomass))
sum_lin <- summary(fit_lin)
sum_log <- summary(fit_log)
k1 <- sum_lin$coefficients[1]
k2 <- sum_log$coefficients[2]
b2 <- sum_log$coefficients[1]
cat("Coefficients:",k1, k2, b2, "\n")

# Read in samples metadata
samples_se <- read.delim(paste0(path,"samples_metadata_malaise_SE.tsv"))
samples_mg <- read.delim(paste0(path,"samples_metadata_malaise_MG.tsv"))

# Remove NAs and zeros for biomass
samples_se <- samples_se[!is.na(samples_se$biomass_grams) & samples_se$biomass_grams>0,]
samples_mg <- samples_mg[!is.na(samples_mg$biomass_grams) & samples_mg$biomass_grams>0,]

# Compute proportional estimates from lin-lin regression equation
# y = k1 x
prop_spec_se <- sum(samples_se$biomass_grams) * k1
prop_spec_mg <- sum(samples_mg$biomass_grams) * k1

# Compute non-proportional estimates from log-log regression equation
# log10 y = k2 x + b2
nonprop_spec_se <- sum(10^(k2 * log10(samples_se$biomass_grams) + b2))
nonprop_spec_mg <- sum(10^(k2 * log10(samples_mg$biomass_grams) + b2))

# Print results
cat("Specimen estimates for the proportional (lin-lin) model:\n")
cat("Sweden:", prop_spec_se, " specimens\n")
cat("Madagascar:", prop_spec_mg, " specimens\n")

cat("Specimen estimates for the nonproportional (log-log) model:\n")
cat("Sweden:", nonprop_spec_se, " specimens\n")
cat("Madagascar:", nonprop_spec_mg, " specimens\n")

