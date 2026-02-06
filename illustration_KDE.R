## Generate example KDE

rm(list=ls())

## Path and libraries
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")
# Set seed for reproducibility
set.seed(123)

# Generate 100 normally distributed random numbers
data <- rnorm(100, mean = 0, sd = 1)

# Create a data frame for the histogram
hist_data <- data.frame(value = data)

# Perform KDE with different bandwidths
kde_005 <- density(data, bw = 0.05)
kde_01 <- density(data, bw = 0.1)
kde_03 <- density(data, bw = 0.3)

# Create data frames for the KDEs
kde_005_df <- data.frame(x = kde_005$x, y = kde_005$y, bandwidth = "0.05")
kde_01_df <- data.frame(x = kde_01$x, y = kde_01$y, bandwidth = "0.1")
kde_03_df <- data.frame(x = kde_03$x, y = kde_03$y, bandwidth = "0.3")

# Combine KDE data frames
kde_df <- bind_rows(kde_005_df, kde_01_df, kde_03_df)

# Plot histogram and KDEs
kde_example <- ggplot(hist_data, aes(x = value)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#9ed6e3", color = "black", alpha = 0.7) +
  geom_line(data = kde_df, aes(x = x, y = y, color = bandwidth, linetype = bandwidth), size = 1.2) + theme_bw() +
  labs(x = "Value", y = "Density") +
  scale_color_manual(name = "Bandwidth h", values = c("0.05" = "#ed6e6c", "0.1" = "#9dc544", "0.3" = "#423089")) +
  scale_linetype_manual(name = "Bandwidth h", values = c("0.05" = "solid", "0.1" = "solid", "0.3" = "solid"))

##########################
###    SAVE PLOTS
## to tikZ
setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
tikz('kde_example.tex', standAlone = TRUE, width=3.5, height=2.0)
kde_example
dev.off()
tools::texi2dvi('kde_example.tex',pdf=T) ## LUALATEX

