## Load packages
library(tidyverse)
library(pracma)
library(ggplot2)

## Import the LabChip output files
data <- read.csv("main-run_SizeTable.csv")
data <- data[, c("Time", "Size", "Alpha4")] %>%
  rename(Sample = Alpha4)

## Import the Sequence data
freq_df <- read.csv("freq_df/freq_df_Alpha4.csv")

## Data cleaning
data[data < 0] <- 0 # replace negative values with 0
data <- subset(data, Time >= 20 & Time <= 72) # filter by running time 

ggplot() +
  geom_line(data = data, aes(x = Time, y = Sample, colour = "heteroduplex")) +  
  labs(title = "Heteroduplex frequencies", 
       x = "Time",
       y = "Fluorescence",
       color = "Legend") 

ggplot() +
  geom_bar(data = freq_df, aes(x = pairwise_distance, y = frxfr, colour = "freq*freq"), stat = "identity") +   
  labs(title = "Seqeunce data of Zymobiomics", 
       x = "Pairwise distance",
       y = "Frequency * frequency",
       color = "Legend") 


####### Normalisation ######
## Lower marker time
LM_data <- subset(data, data$Size >= 0 & data$Size < 0.2)
all <- findpeaks(LM_data$Sample, minpeakheight = 50)
LM <- LM_data$Time[LM_data$Sample == all[1, 1]]

## Custom marker (CO1) time
CO1_data <- subset(data, data$Size > 0.7 & data$Size < 1.0)
all <- findpeaks(CO1_data$Sample, minpeakheight = 5)
CM <- CO1_data$Time[CO1_data$Sample == all[1, 1]]
# CM <- 

## Homoduplex peak time
homo_data <- subset(data, data$Size > 0.8 & data$Size < 12.0)
# Manually find the start and endpoint of homoduplex
homo_start <- 45.86
homo_end <- 47.55
homo <- (homo_start + homo_end) / 2 

## Align the migration time using two markers and the homoduplex
## gdH = (heteroduplex - homoduplex) / (custom marker - lower marker)
data$gdH <- (data$Time - homo) / (CM - LM) 

# Normalisation of fluorescence intensity 
data <- subset(data, Time >= homo_start)
data$Sample <- data$Sample / sum(data$Sample) 

# remove homoduplex
data <- subset(data, Time >= homo_end)

ggplot() +
  geom_line(data = data, aes(x = gdH, y = Sample, colour = "heteroduplex")) +  
  labs(title = "Heteroduplex frequencies", 
       x = "Time",
       y = "Fluorescence",
       color = "Legend") 

## Optimisation
data$opti_gdH <- (-0.004366542) + 0.242374649 * data$gdH
data$opti_flu <- data$Sample * (5.484994353 * exp(8.542979393 * data$opti_gdH))

## Calculate HMAa
HMAa <- sum(data$opti_gdH * data$opti_flu)
HMAa

## plot
ggplot() +
  geom_line(data = data, aes(x = opti_gdH, y = Sample, colour = "heteroduplex")) +
  geom_line(data = data, aes(x = opti_gdH, y = opti_flu, colour = "optimised_heteroduplex")) +
  geom_bar(data = freq_df, aes(x = pairwise_distance, y = frxfr, colour = "nucleotide"), stat = "identity") + 
  theme(legend.position = "none") +
  labs(title = "Alpha3", 
       x = "pairwise_distance",
       y = "freq*freq",
       color = "Legend") 

