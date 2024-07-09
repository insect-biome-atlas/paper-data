############################################
# Required packages: 
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid) 

###########################
# Set data path
##########################
# TODO: set to figshare repo before publication
# path <- "~/dev/figshare-repos/iba/raw_data/"

###########################
# Weights - read in data
##########################
weights_IBA <- read.delim("../biomass_count_IBA.tsv", dec=",")
weights_SIIP <- read.delim("../biomass_count_SIIP.tsv", dec=",")
colnames(weights_IBA)[3:4] <- c("Wei","Ind")
colnames(weights_SIIP)[12:13] <- c("Wei","Ind")

###################################
# Weight curves - linear regression
###################################
full_set <- rbind(weights_SIIP[,12:13],weights_IBA[,3:4])

#MODEL 1 (log)
model1 <- lm(log10(Ind) ~ log10(Wei), data = full_set)
model_summary_1 <- summary(model1)

coefficients_1 <- model_summary_1$coefficients
r_squared_1 <- model_summary_1$r.squared
adj_r_squared_1 <- model_summary_1$adj.r.squared
f_statistic_1 <- model_summary_1$fstatistic
p_value_f_1 <- pf(f_statistic_1[1], f_statistic_1[2], f_statistic_1[3], lower.tail = FALSE)
residual_se_1 <- model_summary_1$sigma
conf_intervals_1 <- confint(model1, level = 0.95)

# Print Results
cat("Model 1 (log-log): \n")
cat("Model Equation: y =", coefficients_1[1, 1], "+", coefficients_1[2, 1], "* x\n")
cat("R-squared:", r_squared_1, "\n")
cat("Adjusted R-squared:", adj_r_squared_1, "\n")
cat("P-value for the slope:", coefficients_1[2, 4], "\n")
cat("Standard Error for the slope:", coefficients_1[2, 2], "\n")
cat("F-statistic:", f_statistic_1[1], "\n")
cat("P-value for the model:", p_value_f_1, "\n")
cat("Residual Standard Error:", residual_se_1, "\n")
cat("95% Confidence Intervals:\n")
print(conf_intervals_1)

#MODEL 2
model2 <- lm(Ind ~ 0 + Wei, data = full_set)
model_summary_2 <- summary(model2)

coefficients_2 <- model_summary_2$coefficients
r_squared_2 <- model_summary_2$r.squared
adj_r_squared_2 <- model_summary_2$adj.r.squared
f_statistic_2 <- model_summary_2$fstatistic
p_value_f_2 <- pf(f_statistic_1[1], f_statistic_1[2], f_statistic_1[3], lower.tail = FALSE)
residual_se_2 <- model_summary_2$sigma
conf_intervals_2 <- confint(model2, level = 0.95)

# Print Results
cat("Model 2 (lin-lin): \n")
cat("Model Equation: y =", coefficients_2[1, 1], "* x\n")
cat("R-squared:", r_squared_2, "\n")
cat("Adjusted R-squared:", adj_r_squared_2, "\n")
cat("F-statistic:", f_statistic_2[1], "\n")
cat("P-value for the model:", p_value_f_2, "\n")
cat("Residual Standard Error:", residual_se_2, "\n")
cat("95% Confidence Intervals:\n")
print(conf_intervals_2)

###################################
# Plots
###################################
plotA <- ggplot(full_set, aes(Wei, Ind)) +
  geom_point() +  # Add points for each data pair
  geom_point(data = weights_IBA, color = "red") +  # Overlay points from IBA
  geom_smooth(data=full_set, method = "lm", formula= y ~ 0 + x, se = TRUE) +
  annotate("text", x = 80, y = 52000, label = sprintf(paste0("y = ", round(coefficients_2[1], digits=0)," * x")), vjust = 3, hjust = 1, color = "blue") +
  annotate("text", x = 140, y = 49000, label = sprintf(paste0("95%% CI of slope: [", round(conf_intervals_2[1,1],digits=0),", ", round(conf_intervals_2[1,2],digits=0),"]")), vjust = 3, hjust = 1, color = "blue") +
  labs(title = "A.",
       y = "No. of individuals",
       x = "Biomass (g)") +
  theme_minimal()+
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
          scale_y_continuous(
    breaks = (c(10000, 20000, 30000, 40000, 50000)),  
    labels = c("10'000", "20'000", "30'000", "40'000","50'000")  
  ) 

plotB <- ggplot(full_set, aes(log10(Wei), log10(Ind))) +
  geom_point() +  # Add points for each data pair
  geom_point(data = weights_IBA, aes(x = log10(Wei), y = log10(Ind)), color = "red") +  # Overlay points from IBA
  geom_smooth(data=full_set, method = "lm", formula= y ~ x , se = TRUE) +
  annotate("text", x = 1.35, y = log10(100000), label = sprintf(paste0("log10 y = ", round(coefficients_1[2], digits=2)," * log10 x + ", round(coefficients_1[1], digits=2))), vjust = 3, hjust = 1, color = "blue") +
  annotate("text", x = 1.2, y = log10(60000), label = sprintf(paste0("95%% CI of slope: [", round(conf_intervals_1[2,1],digits=2),", ", round(conf_intervals_1[2,2],digits=2),"]")), vjust = 3, hjust = 1, color = "blue") +
  labs(title = "B.",
       y = "No. of individuals",
       x = "Biomass (g)") +
  theme_minimal()+
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
      scale_y_continuous(
    breaks = log10(c(1, 10, 100, 1000,10000)),  
    labels = c("1", "10", "100", "1000","10'000")  
  ) +
  scale_x_continuous(
    breaks = log10(c(1, 10, 100)),  
    labels = c("1", "10", "100")  
  )


###################################
# Output figures
###################################
combined_plot <- arrangeGrob(plotA, plotB, ncol = 2)
jpeg("plot_biomass_specimens.jpg", width = 2600, height = 1400, res = 300)
grid.draw(combined_plot)
dev.off()
###################################
