## Load packages
library(tidyverse)  
library(pracma)
library(ggplot2)

## Import the LabChip output files
data <- read.csv("mock_SizeTable.csv")

freq_SA <- read.csv("SA_data.csv")
new_row <- data.frame(pair = 'baseline', pairwise_distance = 0.24164, dH = 0.62, frxfr = 0)
freq_SA  <- rbind(freq_SA, new_row)

freq_ZB <- read.csv("ZB_data.csv")
new_row <- data.frame(pair = 'baseline', pairwise_distance = 0.2156, dH = 0.55, frxfr = 0)
freq_ZB  <- rbind(freq_ZB, new_row)

freq_Mix <- read.csv("Mix_data.csv")
new_row <- data.frame(pair = 'baseline', pairwise_distance = 0.24164, dH = 0.62, frxfr = 0)
freq_Mix  <- rbind(freq_Mix, new_row)

## Data cleaning
data[data < 0] <- 0 # replace negative values with 0

data_SA <- subset(data, Time >= 20 & Time <= 62.27) # filter by running time 
data_ZB <- subset(data, Time >= 20 & Time <= 63.10) # filter by running time 
data_Mix <- subset(data, Time >= 20 & Time <= 64.13) # filter by running time 

####### Calculation of Heteroduplex dH ######
## Lower marker time
LM_data <- subset(data, data$Size >= 0 & data$Size < 0.2)
all <- findpeaks(LM_data$SA, minpeakheight = 50)
LM_SA <- LM_data$Time[LM_data$SA == all[1, 1]]
all <- findpeaks(LM_data$ZB, minpeakheight = 50)
LM_ZB <- LM_data$Time[LM_data$ZB == all[1, 1]]
all <- findpeaks(LM_data$Mix, minpeakheight = 50)
LM_Mix <- LM_data$Time[LM_data$Mix == all[1, 1]]

## Custom marker (CO1) time
CO1_data <- subset(data, data$Size > 0.7 & data$Size < 0.8)
all <- findpeaks(CO1_data$SA, minpeakheight = 5)
CM_SA <- CO1_data$Time[CO1_data$SA == all[1, 1]]
all <- findpeaks(CO1_data$ZB, minpeakheight = 5)
CM_ZB <- CO1_data$Time[CO1_data$ZB == all[1, 1]]
all <- findpeaks(CO1_data$Mix, minpeakheight = 5)
CM_Mix <- CO1_data$Time[CO1_data$Mix == all[1, 1]]

## Homoduplex peak time
homo_data <- subset(data, data$Size > 1.4 & data$Size < 1.7)
all <- findpeaks(homo_data$SA, minpeakheight = 50)
homo_SA <- homo_data$Time[homo_data$SA == all[1, 1]]
all <- findpeaks(homo_data$ZB, minpeakheight = 50)
homo_ZB <- homo_data$Time[homo_data$ZB == all[1, 1]]
all <- findpeaks(homo_data$Mix, minpeakheight = 50)
homo_Mix <- homo_data$Time[homo_data$Mix == all[1, 1]]

## Align the migration time using two markers and the homoduplex
## gdH = (heteroduplex - homoduplex) / (custom marker - lower marker)
data_SA$gdH_SA = (data_SA$Time - homo_SA) / (CM_SA - LM_SA) 
data_SA <- subset(data_SA, Time >= homo_SA)
data_SA$SA <- data_SA$SA / sum(data_SA$SA) # Normalisation
data_SA <- subset(data_SA, Time >= 47.88) # Remove homoduplex

data_ZB$gdH_ZB = (data_ZB$Time - homo_ZB) / (CM_ZB - LM_ZB) 
data_ZB <- subset(data_ZB, Time >= homo_ZB)
data_ZB$ZB <- data_ZB$ZB / sum(data_ZB$ZB) # Normalisation
data_ZB <- subset(data_ZB, Time >= 46.95) # Remove homoduplex

data_Mix$gdH_Mix = (data_Mix$Time - homo_Mix) / (CM_Mix - LM_Mix) 
data_Mix <- subset(data_Mix, Time >= homo_Mix)
data_Mix$Mix <- data_Mix$Mix / sum(data_Mix$Mix) # Normalisation
data_Mix <- subset(data_Mix, Time >= 46.75) # Remove homoduplex
#############################################

###### Binning #######
# Divide dH values into 300 bins
data_SA$bin <- cut(data_SA$gdH_SA, breaks = 300, labels = FALSE)
freq_SA$bin <- cut(freq_SA$pairwise_distance, breaks = 300, labels = FALSE)
data_ZB$bin <- cut(data_ZB$gdH_ZB, breaks = 300, labels = FALSE)
freq_ZB$bin <- cut(freq_ZB$pairwise_distance, breaks = 300, labels = FALSE)
data_Mix$bin <- cut(data_Mix$gdH_Mix, breaks = 300, labels = FALSE)
freq_Mix$bin <- cut(freq_Mix$pairwise_distance, breaks = 300, labels = FALSE)

# Calculate min and max for each dataframe
gdH_SA_min <- min(data_SA$gdH_SA, na.rm = TRUE)
gdH_SA_max <- max(data_SA$gdH_SA, na.rm = TRUE)
gdH_ZB_min <- min(data_ZB$gdH_ZB, na.rm = TRUE)
gdH_ZB_max <- max(data_ZB$gdH_ZB, na.rm = TRUE)
gdH_Mix_min <- min(data_Mix$gdH_Mix, na.rm = TRUE)
gdH_Mix_max <- max(data_Mix$gdH_Mix, na.rm = TRUE)

pd_SA_min <- min(freq_SA$pairwise_distance, na.rm = TRUE)
pd_SA_max <- max(freq_SA$pairwise_distance, na.rm = TRUE)
pd_ZB_min <- min(freq_ZB$pairwise_distance, na.rm = TRUE)
pd_ZB_max <- max(freq_ZB$pairwise_distance, na.rm = TRUE)
pd_Mix_min <- min(freq_Mix$pairwise_distance, na.rm = TRUE)
pd_Mix_max <- max(freq_Mix$pairwise_distance, na.rm = TRUE)

# Calculate bin boundaries
SA_bins <- seq(gdH_SA_min, gdH_SA_max, length.out = 301)
ZB_bins <- seq(gdH_ZB_min, gdH_ZB_max, length.out = 301)
Mix_bins <- seq(gdH_Mix_min, gdH_Mix_max, length.out = 301)
SAfreq_bins <- seq(pd_SA_min, pd_SA_max, length.out = 301)
ZBfreq_bins <- seq(pd_ZB_min, pd_ZB_max, length.out = 301)
Mixfreq_bins <- seq(pd_Mix_min, pd_Mix_max, length.out = 301)

# Initialize Heteroduplex_binned dataframe
SA_Heteroduplex_binned <- data.frame(
  gdH_mean = numeric(300),
  SA = numeric(300)
)

ZB_Heteroduplex_binned <- data.frame(
  gdH_mean = numeric(300),
  ZB = numeric(300)
)

Mix_Heteroduplex_binned <- data.frame(
  gdH_mean = numeric(300),
  Mix = numeric(300)
)

# Initialize Nucleotide_binned dataframe
SA_Nucleotide_binned <- data.frame(
  pd_mean = numeric(300),
  frxfr = numeric(300)
)

ZB_Nucleotide_binned <- data.frame(
  pd_mean = numeric(300),
  frxfr = numeric(300)
)

Mix_Nucleotide_binned <- data.frame(
  pd_mean = numeric(300),
  frxfr = numeric(300)
)

# Calculate mean heteroduplex frequency for each bin
for (b in 1:300) {
  bin_data <- data_SA[data_SA$bin == b, ]
  SA_Heteroduplex_binned$SA[b] <- if (nrow(bin_data) > 0) {
    round(mean(bin_data$SA, na.rm = TRUE), 4)
  } else {
    0.000001
  }
  
  # Calculate gdH_mean for each bin using bin boundaries
  SA_Heteroduplex_binned$gdH_mean[b] <- mean(SA_bins[b:(b + 1)], na.rm = TRUE)
}

for (b in 1:300) {
  bin_data <- data_ZB[data_ZB$bin == b, ]
  ZB_Heteroduplex_binned$ZB[b] <- if (nrow(bin_data) > 0) {
    round(mean(bin_data$ZB, na.rm = TRUE), 4)
  } else {
    0.000001
  }
  
  # Calculate gdH_mean for each bin using bin boundaries
  ZB_Heteroduplex_binned$gdH_mean[b] <- mean(ZB_bins[b:(b + 1)], na.rm = TRUE)
}

for (b in 1:300) {
  bin_data <- data_Mix[data_Mix$bin == b, ]
  Mix_Heteroduplex_binned$Mix[b] <- if (nrow(bin_data) > 0) {
    round(mean(bin_data$Mix, na.rm = TRUE), 4)
  } else {
    0.000001
  }
  
  # Calculate gdH_mean for each bin using bin boundaries
  Mix_Heteroduplex_binned$gdH_mean[b] <- mean(Mix_bins[b:(b + 1)], na.rm = TRUE)
}

## Calculate mean frxfr for each bin
for (b in 1:300) {
  bin_freq_data <- freq_SA[freq_SA$bin == b, ]
  SA_Nucleotide_binned$frxfr[b] <- if (nrow(bin_freq_data) > 0) {
    round(mean(bin_freq_data$frxfr, na.rm = TRUE), 4)
  } else {
    NA
  }
  
  # Calculate pd_mean for each bin using bin boundaries
  SA_Nucleotide_binned$pd_mean[b] <- mean(SAfreq_bins[b:(b + 1)], na.rm = TRUE)
}

for (b in 1:300) {
  bin_freq_data <- freq_ZB[freq_ZB$bin == b, ]
  ZB_Nucleotide_binned$frxfr[b] <- if (nrow(bin_freq_data) > 0) {
    round(mean(bin_freq_data$frxfr, na.rm = TRUE), 4)
  } else {
    NA
  }
  
  # Calculate pd_mean for each bin using bin boundaries
  ZB_Nucleotide_binned$pd_mean[b] <- mean(ZBfreq_bins[b:(b + 1)], na.rm = TRUE)
}

for (b in 1:300) {
  bin_freq_data <- freq_Mix[freq_Mix$bin == b, ]
  Mix_Nucleotide_binned$frxfr[b] <- if (nrow(bin_freq_data) > 0) {
    round(mean(bin_freq_data$frxfr, na.rm = TRUE), 4)
  } else {
    NA
  }
  
  # Calculate pd_mean for each bin using bin boundaries
  Mix_Nucleotide_binned$pd_mean[b] <- mean(Mixfreq_bins[b:(b + 1)], na.rm = TRUE)
}

# Exclude NA values from the data
SA_valid_data <- na.omit(SA_Nucleotide_binned)
ZB_valid_data <- na.omit(ZB_Nucleotide_binned)
Mix_valid_data <- na.omit(Mix_Nucleotide_binned)

## Linear interpolation
# Use the 'approx' function for linear interpolation of missing frxfr values
SA_Nucleotide_binned$interpolated <- SA_Nucleotide_binned$frxfr  # Create a new column
ZB_Nucleotide_binned$interpolated <- ZB_Nucleotide_binned$frxfr  
Mix_Nucleotide_binned$interpolated <- Mix_Nucleotide_binned$frxfr  

# Perform interpolation
SA_Nucleotide_binned$interpolated[is.na(SA_Nucleotide_binned$frxfr)] <-
  approx(SA_valid_data$pd_mean, SA_valid_data$frxfr,
         xout = SA_Nucleotide_binned$pd_mean[is.na(SA_Nucleotide_binned$frxfr)])$y

ZB_Nucleotide_binned$interpolated[is.na(ZB_Nucleotide_binned$frxfr)] <-
  approx(ZB_valid_data$pd_mean, ZB_valid_data$frxfr,
         xout = ZB_Nucleotide_binned$pd_mean[is.na(ZB_Nucleotide_binned$frxfr)])$y

Mix_Nucleotide_binned$interpolated[is.na(Mix_Nucleotide_binned$frxfr)] <-
  approx(Mix_valid_data$pd_mean, Mix_valid_data$frxfr,
         xout = Mix_Nucleotide_binned$pd_mean[is.na(Mix_Nucleotide_binned$frxfr)])$y

###### Optimisation ######
# Set variable names
gdH_SA <- SA_Heteroduplex_binned$gdH_mean
gdH_ZB <- ZB_Heteroduplex_binned$gdH_mean
gdH_Mix <- Mix_Heteroduplex_binned$gdH_mean
npd_SA <- SA_Nucleotide_binned$pd_mean
npd_ZB <- ZB_Nucleotide_binned$pd_mean
npd_Mix <- Mix_Nucleotide_binned$pd_mean

hfreq_SA <- (SA_Heteroduplex_binned$SA) 
hfreq_ZB <- (ZB_Heteroduplex_binned$ZB) 
hfreq_Mix <- (Mix_Heteroduplex_binned$Mix) 
nfreq_SA <- SA_Nucleotide_binned$interpolated
nfreq_ZB <- ZB_Nucleotide_binned$interpolated
nfreq_Mix <- Mix_Nucleotide_binned$interpolated

## Find the best model for the fluorescence adjustment (i.e. mismatch optimisation)
## 1. Optimise the linear regression model 
obj_lm <- function(param) {
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]

  estimated_npd_SA <- a + b * gdH_SA
  estimated_npd_ZB <- a + b * gdH_ZB
  estimated_npd_Mix <- a + b * gdH_Mix
  
  estimated_nfreq_SA <- hfreq_SA * (c + d * estimated_npd_SA)
  estimated_nfreq_ZB <- hfreq_ZB * (c + d * estimated_npd_ZB)
  estimated_nfreq_Mix <- hfreq_Mix * (c + d * estimated_npd_Mix)
  
  resid_SA <- sum((nfreq_SA - estimated_nfreq_SA)^2) + sum((npd_SA - estimated_npd_SA)^2)
  resid_ZB <- sum((nfreq_ZB - estimated_nfreq_ZB)^2) + sum((npd_ZB - estimated_npd_ZB)^2)
  resid_Mix <- sum((nfreq_Mix - estimated_nfreq_Mix)^2) + sum((npd_Mix - estimated_npd_Mix)^2)
  
  c1 <- sqrt(resid_SA + resid_ZB + resid_Mix)
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
  
  estimated_npd_SA <- a + b * gdH_SA
  estimated_npd_ZB <- a + b * gdH_ZB
  estimated_npd_Mix <- a + b * gdH_Mix
  
  estimated_nfreq_SA <- hfreq_SA * (c * exp(d * estimated_npd_SA))
  estimated_nfreq_ZB <- hfreq_ZB * (c * exp(d * estimated_npd_ZB))
  estimated_nfreq_Mix <- hfreq_Mix * (c * exp(d * estimated_npd_Mix))
  
  resid_SA <- sum((nfreq_SA - estimated_nfreq_SA)^2) + sum((npd_SA - estimated_npd_SA)^2)
  resid_ZB <- sum((nfreq_ZB - estimated_nfreq_ZB)^2) + sum((npd_ZB - estimated_npd_ZB)^2)
  resid_Mix <- sum((nfreq_Mix - estimated_nfreq_Mix)^2) + sum((npd_Mix - estimated_npd_Mix)^2)
  
  c2 <- sqrt(resid_SA + resid_ZB + resid_Mix)
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
  
  estimated_npd_SA <- a + b * gdH_SA
  estimated_npd_ZB <- a + b * gdH_ZB
  estimated_npd_Mix <- a + b * gdH_Mix
  
  estimated_nfreq_SA <- e * hfreq_SA * (dlnorm(x = estimated_npd_SA, meanlog = d, sdlog = c))
  estimated_nfreq_ZB <- e * hfreq_ZB * (dlnorm(x = estimated_npd_ZB, meanlog = d, sdlog = c))
  estimated_nfreq_Mix <- e * hfreq_Mix * (dlnorm(x = estimated_npd_Mix, meanlog = d, sdlog = c))
  
  resid_SA <- sum((nfreq_SA - estimated_nfreq_SA)^2) + sum((npd_SA - estimated_npd_SA)^2)
  resid_ZB <- sum((nfreq_ZB - estimated_nfreq_ZB)^2) + sum((npd_ZB - estimated_npd_ZB)^2)
  resid_Mix <- sum((nfreq_Mix - estimated_nfreq_Mix)^2) + sum((npd_Mix - estimated_npd_Mix)^2)
  
  c3 <- sqrt(resid_SA + resid_ZB + resid_Mix)
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
  
  estimated_npd_SA <- a + b * gdH_SA
  estimated_npd_ZB <- a + b * gdH_ZB
  estimated_npd_Mix <- a + b * gdH_Mix
  
  estimated_nfreq_SA <- e * hfreq_SA * (dgamma(x = estimated_npd_SA, shape = c, scale = d))
  estimated_nfreq_ZB <- e * hfreq_ZB * (dgamma(x = estimated_npd_ZB, shape = c, scale = d))
  estimated_nfreq_Mix <- e * hfreq_Mix * (dgamma(x = estimated_npd_Mix, shape = c, scale = d))
  
  resid_SA <- sum((nfreq_SA - estimated_nfreq_SA)^2) + sum((npd_SA - estimated_npd_SA)^2)
  resid_ZB <- sum((nfreq_ZB - estimated_nfreq_ZB)^2) + sum((npd_ZB - estimated_npd_ZB)^2)
  resid_Mix <- sum((nfreq_Mix - estimated_nfreq_Mix)^2) + sum((npd_Mix - estimated_npd_Mix)^2)
  
  c4 <- sqrt(resid_SA + resid_ZB + resid_Mix)
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
RMSE <- 0.5038972
SA_Heteroduplex_binned$estimated_npd <- (-0.007638021) + 0.244851730 * gdH_SA
SA_Heteroduplex_binned$estimated_nfreq <- hfreq_SA * (4.720633550 * exp(8.577815568 * SA_Heteroduplex_binned$estimated_npd))

# Create the plot
ggplot() +
  geom_line(data = SA_Heteroduplex_binned, aes(x = estimated_npd , y = estimated_nfreq, colour = "heteroduplex")) + 
  geom_bar(data = SA_Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), 
           stat = "identity") +  
  labs(title = "(a)", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

######
ZB_Heteroduplex_binned$estimated_npd <- (-0.007638021) + 0.244851730 * gdH_ZB
ZB_Heteroduplex_binned$estimated_nfreq <- hfreq_ZB * (4.720633550 * exp(8.577815568 * ZB_Heteroduplex_binned$estimated_npd))

# Create the plot
ggplot() +
  geom_line(data = ZB_Heteroduplex_binned, aes(x = estimated_npd , y = estimated_nfreq, colour = "heteroduplex")) + 
  geom_bar(data = ZB_Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), 
           stat = "identity") +  
  labs(title = "(b)", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

######

Mix_Heteroduplex_binned$estimated_npd <- (-0.007638021) + 0.244851730 * gdH_Mix
Mix_Heteroduplex_binned$estimated_nfreq <- hfreq_Mix * (4.720633550 * exp(8.577815568 * Mix_Heteroduplex_binned$estimated_npd))

# Create the plot
ggplot() +
  geom_line(data = Mix_Heteroduplex_binned, aes(x = estimated_npd , y = estimated_nfreq, colour = "heteroduplex")) + 
  geom_bar(data = Mix_Nucleotide_binned, aes(x = pd_mean, y = frxfr, colour = "nucleotide"), 
           stat = "identity") +  
  labs(title = "(c)", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 


