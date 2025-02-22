## Load packages
library(tidyverse)
library(pracma)
library(ggplot2)
library(gridExtra)

## Import the LabChip output files
data <- read.csv("mock_SizeTable.csv")
freq_df <- read.csv("Mix_data.csv")

row <- data.frame(pair = 'baseline', pairwise_distance = 0.24164, dH = 0.62, frxfr = 0)
freq_df  <- rbind(freq_df, row)

## Data cleaning
data <- subset(data, Time >= 20 & Time <= 64.13) # filter by running time 72
data[data < 0] <- 0 # replace negative values with 0

ggplot() +
  geom_line(data = data, aes(x = Time, y = Mix, colour = "heteroduplex")) +  
  labs(title = "Heteroduplex frequencies", 
       x = "Time",
       y = "Fluorescence",
       color = "Legend") 

####### Calculation of Heteroduplex dH ######
## Lower marker time
LM_data <- subset(data, data$Size >= 0 & data$Size < 0.2)
all <- findpeaks(LM_data$Mix, minpeakheight = 50)
LM <- LM_data$Time[LM_data$Mix == all[1, 1]]

## Custom marker (CO1) time
CO1_data <- subset(data, data$Size > 0.7 & data$Size < 0.8)
all <- findpeaks(CO1_data$Mix, minpeakheight = 5)
CM <- CO1_data$Time[CO1_data$Mix == all[1, 1]]

## Homoduplex peak time
homo_data <- subset(data, data$Size > 1.4 & data$Size < 1.7)
all <- findpeaks(homo_data$Mix, minpeakheight = 50)
homo <- homo_data$Time[homo_data$Mix == all[1, 1]]

## Align the migration time using two markers and the homoduplex
## gdH = (heteroduplex - homoduplex) / (custom marker - lower marker)
data$gdH_Mix = (data$Time - homo) / (CM - LM) 
data <- subset(data, Time >= homo)

## Normalisation
data$Mix <- data$Mix / sum(data$Mix)

## Remove homoduplex
data <- subset(data, Time >= 46.75)
#############################################

###### Binning #######
# Divide dH values into 300 bins
data$bin <- cut(data$gdH_Mix, breaks = 300, labels = FALSE)
freq_df$bin <- cut(freq_df$pairwise_distance, breaks = 300, labels = FALSE)

# Calculate min and max for each dataframe
gdH_Mix_min <- min(data$gdH_Mix, na.rm = TRUE)
gdH_Mix_max <- max(data$gdH_Mix, na.rm = TRUE)
pd_min <- min(freq_df$pairwise_distance, na.rm = TRUE)
pd_max <- max(freq_df$pairwise_distance, na.rm = TRUE)

# Calculate bin boundaries
data_bins <- seq(gdH_Mix_min, gdH_Mix_max, length.out = 301)
freq_bins <- seq(pd_min, pd_max, length.out = 301)

# Initialize Heteroduplex_binned dataframe
Heteroduplex_binned <- data.frame(
  gdH_mean = numeric(300),
  Mix = numeric(300)
)

# Initialize Nucleotide_binned dataframe
Nucleotide_binned <- data.frame(
  pd_mean = numeric(300),
  frxfr = numeric(300)
)

# Calculate mean heteroduplex frequency for each bin
for (b in 1:300) {
  bin_data <- data[data$bin == b, ]
  Heteroduplex_binned$Mix[b] <- if (nrow(bin_data) > 0) {
    round(mean(bin_data$Mix, na.rm = TRUE), 4)
  } else {
    0.000001
  }
  
  # Calculate gdH_mean for each bin using bin boundaries
  Heteroduplex_binned$gdH_mean[b] <- mean(data_bins[b:(b + 1)], na.rm = TRUE)
}

# Calculate mean frxfr for each bin
for (b in 1:300) {
  bin_freq_data <- freq_df[freq_df$bin == b, ]
  Nucleotide_binned$frxfr[b] <- if (nrow(bin_freq_data) > 0) {
    round(mean(bin_freq_data$frxfr, na.rm = TRUE), 4)
  } else {
    NA
  }
  
  # Calculate pd_mean for each bin using bin boundaries
  Nucleotide_binned$pd_mean[b] <- mean(freq_bins[b:(b + 1)], na.rm = TRUE)
}

# Exclude NA values from the data
valid_data <- na.omit(Nucleotide_binned)

# Perform interpolation
# Use the 'approx' function for linear interpolation of missing frxfr values
Nucleotide_binned$interpolated <- Nucleotide_binned$frxfr  # Create a new column

## Linear interpolation
Nucleotide_binned$interpolated[is.na(Nucleotide_binned$frxfr)] <-
  approx(valid_data$pd_mean, valid_data$frxfr,
         xout = Nucleotide_binned$pd_mean[is.na(Nucleotide_binned$frxfr)])$y

## Spline interpolation
# Calculate interpolated values
# interpolated_values <- spline(valid_data$dH_mean, valid_data$frxfr, 
#                               xout = Nucleotide_binned$dH_mean[is.na(Nucleotide_binned$frxfr)])$y

# Set interpolated values to 0 if they are less than or equal to 0
# Nucleotide_binned$interpolated[is.na(Nucleotide_binned$frxfr)] <- 
#   pmax(interpolated_values,  0.000001)

# Graph after interpolation
ggplot(Nucleotide_binned, aes(x = pd_mean)) +
  geom_point(aes(y = frxfr, color = "original"), size = 3) +  # Original values (NA values are not displayed)
  geom_point(aes(y = interpolated, color = "interpolation")) +  # Interpolated values
  labs(title = "Nucleotide Binned Data with Interpolation", x = "pd_mean", y = "frxfr") +
  geom_line(aes(y = interpolated, color = "interpolation"), linetype = "dashed") 

## Graph before optimisation
ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = gdH_mean, y = Mix, colour = "heteroduplex")) + 
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), 
           stat = "identity") +   
  labs(title = "Heteroduplex frequencies and Nucleotide frequencies", 
       x = "relative time (gdH or pd)",
       y = "freq*freq",
       color = "Legend") 

## Graph before optimisation
ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = gdH_mean, y = Mix, colour = "heteroduplex")) +  
  geom_line(data = Nucleotide_binned, aes(x = pd_mean, y = interpolated, colour = "nucleotide")) +  
  labs(title = "Heteroduplex frequencies and Nucleotide frequencies", 
       x = "relative time (gdH or pd)",
       y = "freq*freq",
       color = "Legend") 
######################

###### Optimisation ######
# Set variable names
npd <- Nucleotide_binned$pd_mean
gdH <- Heteroduplex_binned$gdH_mean
dna_freq <- Nucleotide_binned$interpolated
hma_freq <- Heteroduplex_binned$Mix

## Find the best model for the fluorescence adjustment 
## 1. Optimise the linear regression model 
obj_lm <- function(param) {
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  
  estimated_npd <- a + b * gdH
  estimated_nfreq <- hma_freq * (c + d * estimated_npd)
  
  c1 <- sqrt(sum((dna_freq - estimated_nfreq)^2) + sum((npd - estimated_npd)^2))
  return(c1)
}
opti_lm <- optim(par = c(1, 1, 1, 1), fn = obj_lm, method = "L-BFGS-B", control = list(maxit = 1e9, fnscale = 1),
                 lower = c(-Inf, 1e-10, 1e-10, -Inf))

## 2. Optimise the exponential distribution model
obj_exp <- function(param) {
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  
  estimated_npd <- a + b * gdH
  estimated_nfreq <- hma_freq * (c * exp(d *  estimated_npd))
  
  c2 <- sqrt(sum((dna_freq - estimated_nfreq)^2) + sum((npd - estimated_npd)^2))
  return(c2)
}
opti_exp <- optim(par = c(1, 1, 1, 1), fn = obj_exp, method = "L-BFGS-B", control = list(maxit = 1e9, fnscale = 1),
                  lower = c(-Inf, 1e-10, 1e-10, 1e-10))

## 3. Optimise the log-normal distribution model
obj_logn <- function(param) {
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  e <- param[5]
  
  estimated_npd <- a + b * gdH
  estimated_nfreq <- e * hma_freq * (dlnorm(x = estimated_npd, meanlog = c, sdlog = d))
  
  c3 <- sqrt(sum((dna_freq - estimated_nfreq)^2) + sum((npd - estimated_npd)^2))
  return(c3)
}
opti_logn <- optim(par = c(1, 1, 1, 1, 1), fn = obj_logn, method = "L-BFGS-B", control = list(maxit = 1e9, fnscale = 1), 
                   lower = c(-Inf, 1e-10, 1e-10, -Inf ,1e-10))

## 4. Optimise the gamma distribution model
obj_gamma <- function(param) {
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  e <- param[5]
  
  estimated_npd <- a + b * gdH
  estimated_nfreq <- e * hma_freq * (dgamma(x = estimated_npd, shape = c, scale = d))
  
  c4 <- sqrt(sum((dna_freq -  estimated_nfreq)^2) + sum((npd - estimated_npd)^2))
  return(c4)
}
opti_gamma <- optim(par = c(1, 1, 1, 1, 1), fn = obj_gamma, method = "L-BFGS-B", control = list(maxit = 1e9, fnscale = 1),
                    lower = c(1e-10, 1e-10, 1, 1e-10, -Inf))

## Sort by value
df_opti <- data.frame(matrix(ncol = 2, nrow = 4))
colnames(df_opti) <- c("Methods", "Optimize_value")
df_opti[,1] <- c("Linear", "Exponential", "Log-normal", "Gamma")
df_opti[,2] <- c(opti_lm$value, opti_exp$value, opti_logn$value, opti_gamma$value)
df_opti[order(df_opti$Optimize_value),] #Top one is the best fit model
###############################################################
## After optimisation
# Parameter from the linear interpolation - bin 300
RMSE3 <- 0.2933
Heteroduplex_binned$estimated_npd <- (-0.0074458) + 0.2382800 * gdH
Heteroduplex_binned$estimated_nfreq <- hma_freq * (11.5949748 * exp(5.0376538 * Heteroduplex_binned$estimated_npd))

p3 <- ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = gdH_mean, y = Mix, colour = "heteroduplex")) +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = estimated_nfreq, colour = "optimised_heteroduplex")) +
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  annotate("text", x = Inf, y = Inf, 
           label = paste("RMSE =", RMSE3), 
           hjust = 1, vjust = 1, size = 4, color = "black") +
  theme(legend.position = "none") +
  labs(title = "SA+ZB parameter", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

q3 <- ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = Mix, colour = "heteroduplex")) +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = estimated_nfreq, colour = "optimised_heteroduplex")) +
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  annotate("text", x = -Inf, y = Inf, 
           label = paste("RMSE =", RMSE3), 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme(legend.position = "none") +
  labs(title = "(c) SA+ZB parameter", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

data$opti_gdH_Mix <- (-0.0074458) + 0.2382800 * data$gdH_Mix
data$opti_Mix <- data$Mix * (11.5949748 * exp(5.0376538 * data$opti_gdH_Mix))

ggplot() +
  geom_line(data = data, aes(x = opti_gdH_Mix , y = opti_Mix, colour = "heteroduplex")) + 
  geom_bar(data = freq_df, aes(x = pairwise_distance, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  geom_text(aes(x = 0, y = 0.06, label = paste("RMSE =", RMSE3)), size = 4, color = "blue", hjust = 0) +
  labs(title = "without bin", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 
######
# SA model & parameter
Heteroduplex_binned$estimated_npd <- (-0.01248065) + 0.27736854 * gdH
Heteroduplex_binned$estimated_nfreq <- hma_freq * (0.34416230 * exp(18.80257762 * Heteroduplex_binned$estimated_npd))
RMSE1 <- sqrt(sum((dna_freq - Heteroduplex_binned$estimated_nfreq)^2) + sum((npd - Heteroduplex_binned$estimated_npd)^2)) %>%
  round(., digits = 4)

# Create the plot
p1 <- ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = gdH_mean, y = Mix, colour = "heteroduplex")) +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = estimated_nfreq, colour = "optimised_heteroduplex")) +
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), stat = "identity") +
  annotate("text", x = Inf, y = Inf, 
           label = paste("RMSE =", RMSE1), 
           hjust = 1, vjust = 1, size = 4, color = "black") +
  theme(legend.position = "none") +
  labs(title = "SA parameter", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

q1 <- ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = Mix, colour = "heteroduplex")) +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = estimated_nfreq, colour = "optimised_heteroduplex")) +
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  annotate("text", x = -Inf, y = Inf, 
           label = paste("RMSE =", RMSE1), 
           hjust = 0, vjust = 1, size = 4, color = "black") +  theme(legend.position = "none") +
  labs(title = "(a) SA parameter", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

######
# ZB model & parameter
Heteroduplex_binned$estimated_npd <- (-0.007776805) + 0.229360702 * gdH
Heteroduplex_binned$estimated_nfreq <- hma_freq * (4.608453533 * exp(11.271586081 * Heteroduplex_binned$estimated_npd))
RMSE2 <- sqrt(sum((dna_freq - Heteroduplex_binned$estimated_nfreq)^2) + sum((npd - Heteroduplex_binned$estimated_npd)^2)) %>%
  round(., digits = 4)

# Create the plot
p2 <- ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = gdH_mean, y = Mix, colour = "heteroduplex")) +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = estimated_nfreq, colour = "optimised_heteroduplex")) +
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  annotate("text", x = Inf, y = Inf, 
           label = paste("RMSE =", RMSE2), 
           hjust = 1, vjust = 1, size = 4, color = "black") +
  theme(legend.position = "none") +
  labs(title = "ZB parameter", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

q2 <- ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = Mix, colour = "heteroduplex")) +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = estimated_nfreq, colour = "optimised_heteroduplex")) +
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  annotate("text", x = -Inf, y = Inf, 
           label = paste("RMSE =", RMSE2), 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme(legend.position = "none") +
  labs(title = "(b) ZB parameter", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

######
# all3 model & parameter
Heteroduplex_binned$estimated_npd <- (-0.007638021) + 0.244851730 * gdH
Heteroduplex_binned$estimated_nfreq <- hma_freq * (4.720633550 * exp(8.577815568 * Heteroduplex_binned$estimated_npd))
RMSE4 <- sqrt(sum((dna_freq - Heteroduplex_binned$estimated_nfreq)^2) + sum((npd - Heteroduplex_binned$estimated_npd)^2)) %>%
  round(., digits = 4)

# Create the plot
p4 <- ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = gdH_mean, y = Mix, colour = "heteroduplex")) +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = estimated_nfreq, colour = "optimised_heteroduplex")) +
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  annotate("text", x = Inf, y = Inf, 
           label = paste("RMSE =", RMSE4), 
           hjust = 1, vjust = 1, size = 4, color = "black") +
  theme(legend.position = "none") +
  labs(title = "All datasets parameter", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend")  

q4 <- ggplot() +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = Mix, colour = "heteroduplex")) +
  geom_line(data = Heteroduplex_binned, aes(x = estimated_npd, y = estimated_nfreq, colour = "optimised_heteroduplex")) +
  geom_bar(data = Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  annotate("text", x = -Inf, y = Inf, 
           label = paste("RMSE =", RMSE4), 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme(legend.position = "none") +
  labs(title = "(d) All datasets parameter", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

#######
grid.arrange(p1, p2, p3, p4, ncol = 2)
grid.arrange(q1, q2, q3, q4, ncol = 2)

