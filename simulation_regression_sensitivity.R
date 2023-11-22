library(broom)
library(dplyr)

# Set the seed for reproducibility
set.seed(42)

# Generate a normal distribution of ages between 55 to 104
ages <- rnorm(1000, mean = 79.5, sd = 10) # Mean and standard deviation adjusted for the desired range

# Ensure ages are within the range of 55 to 104
ages <- pmin(pmax(ages, 55), 104)

# Sample 1000 individuals from the generated ages
sampled_ages <- sample(ages, 1000, replace = FALSE)

# Print the sampled ages
# print(sampled_ages)


num_individuals <- 1000
num_loci <- 100


af_matrix <- matrix(runif(num_individuals * num_loci, min = 0, max = 0.5), nrow = num_individuals, ncol = num_loci)
af_dataFrame <- as.data.frame(af_matrix)

rownames(af_dataFrame) <- paste0("Ind.", seq(1:1000))
colnames(af_dataFrame) <- paste0("snp.", seq(1:100))
af_dataFrame$age <- c(sampled_ages)

# Calculate correlations between age and specific SNPs
correlations <- cor(af_dataFrame[c("age", "snp.10", "snp.20", "snp.30", "snp.40", "snp.50",
                                   "snp.60", "snp.70", "snp.80", "snp.90", "snp.100")])

# Mutate SNPs to be positively correlated with age
mutated_dataFrame <- af_dataFrame %>%
  mutate(
    snp.10 = snp.10 + 0.5 * correlations["age", "snp.10"] * (age - mean(age)),
    snp.20 = snp.20 + 0.5 * correlations["age", "snp.20"] * (age - mean(age)),
    snp.30 = snp.30 + 0.5 * correlations["age", "snp.30"] * (age - mean(age)),
    snp.40 = snp.40 + 0.5 * correlations["age", "snp.40"] * (age - mean(age)),
    snp.50 = snp.50 + 0.5 * correlations["age", "snp.50"] * (age - mean(age)),
    snp.60 = snp.60 + 0.5 * correlations["age", "snp.60"] * (age - mean(age)),
    snp.70 = snp.70 + 0.5 * correlations["age", "snp.70"] * (age - mean(age)),
    snp.80 = snp.80 + 0.5 * correlations["age", "snp.80"] * (age - mean(age)),
    snp.90 = snp.90 + 0.5 * correlations["age", "snp.90"] * (age - mean(age)),
    snp.100 = snp.100 + 0.5 * correlations["age", "snp.100"] * (age - mean(age))
  )



# Order the dataframe by the age column
af_dataFrame <- mutated_dataFrame[order(mutated_dataFrame$age), ]

# Define the bin size and overlap
# bin_size <- 30
# overlap <- 15

# Define the range of window sizes and overlap sizes
bin_sizes <- c(20, 30, 40, 50)
overlaps <- c(5, 10, 15)

out.df = data.frame()
rows = list()

# Perform sensitivity analysis for different window and overlap sizes
for (bin_size in bin_sizes) {
  for (overlap in overlaps) {
    
    # Calculate the number of bins
    num_bins <- ceiling((nrow(af_dataFrame) - overlap) / (bin_size - overlap))
    
    # Initialize an empty list to store bin means
    bin_means_list <- list()
    
    # Split the data frame into bins and calculate the mean of each bin
    for (i in 1:num_bins) {
      start <- (i - 1) * (bin_size - overlap) + 1
      end <- start + bin_size - 1
      if (end > nrow(af_dataFrame)) {
        end <- nrow(af_dataFrame)
      }
      bin <- af_dataFrame[start:end, ]
      bin_means <- colMeans(bin)  # Include all columns in the calculation of means
      bin_means_list[[i]] <- bin_means
    }
    
    # Combine the bin means into a new data frame
    bin_means_df <- as.data.frame(do.call(rbind, bin_means_list))
    
    # Set column names for the new data frame
    colnames(bin_means_df) <- colnames(af_dataFrame)
    
    
    # Initialize a list to store regression results
    regression_results <- list()
    
    # Perform linear regression for each column against the 'age' column
    for (i in 1:(ncol(bin_means_df) - 1)) {  # Exclude the last column (age) for regression
      formula <- as.formula(paste0(names(bin_means_df)[i], " ~ age"))
      regression <- lm(formula, data = bin_means_df)
      p_value <- summary(regression)$coefficients[2, 4]
      beta_coefficient <- summary(regression)$coefficients[2, 1]
      snp_name <- names(bin_means_df)[i]
      regression_results[[snp_name]] <- c(beta_coefficient, p_value)
    }
    
    # Convert the list of regression results to a data frame
    regression_results_df <- as.data.frame(t(do.call(rbind, regression_results)))
    
    # Set column names for the new data frame
    colnames(regression_results_df) <- c("Beta Coefficient", "P-Value")
    
    # Print the data frame with regression results
    # print(regression_results_df)
    
    
    # Initialize a list to store regression results
    regression_results <- list()
    
    # Perform linear regression for each column against the 'age' column
    for (i in 1:(ncol(bin_means_df) - 1)) {  # Exclude the last column (age) for regression
      formula <- as.formula(paste0(names(bin_means_df)[i], " ~ age"))
      regression <- lm(formula, data = bin_means_df)
      p_value <- tidy(summary(regression))$p.value[2]
      beta_coefficient <- tidy(summary(regression))$estimate[2]
      snp_name <- names(bin_means_df)[i]
      regression_results[[snp_name]] <- c(beta_coefficient, p_value)
    }
    
    # Convert the list of regression results to a data frame
    regression_results_df <- as.data.frame(do.call(rbind, regression_results))
    
    # Set column names for the new data frame
    colnames(regression_results_df) <- c("Beta Coefficient", "P.value")
    
    
    print(paste0("Window_", bin_size, "_Overlap_", overlap))
    print(regression_results_df %>% filter(P.value < 5e-8))
    
    out_p_value <- t(regression_results_df %>% filter(P.value < 5e-8))[2,]
    r.name <- c(as.character(paste0("Window_", bin_size, "_Overlap_", overlap)))
    
    rows <- append(rows, r.name)
    out.df <- rbind(out.df, out_p_value)
  }
}

colnames(out.df) <- c("snp.10", "snp.20", "snp.30", "snp.40", "snp.50","snp.60", "snp.70", "snp.80", "snp.100")
# rownames(out.df) <- rows
out.df$window <- rows
# row.names()

View(out.df)


melt_out_df <- melt(out.df)
colnames(melt_out_df) <- c("window", "snp", "p_value")

ggplot(melt_out_df, aes(snp, window)) +
  geom_tile(aes(fill = -log10(p_value)), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")+
  geom_text(aes(label = p_value), size= 2)+
  theme_bw()



