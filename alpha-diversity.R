###### Load packages ######
library(tidyverse)
library(ape)
library(seqinr)
library(msa) # If failed to install the package, run the following code:
########################################
## if (!requireNamespace("BiocManager", quietly = TRUE))
## install.packages("BiocManager")
## BiocManager::install("msa")
## library(msa)
########################################
library(phyloseq) # If failed to install the package, run the following code:
########################################
## if (!requireNamespace("BiocManager", quietly = TRUE))
## install.packages("BiocManager")
## BiocManager::install("phyloseq")
## library(phyloseq)
########################################
library(VGAM)
library(phyloseq.extended) # If this package is not installed, run the following code:
########################################
##install.packages("remotes")
##remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")
########################################
## Manually define the function if 'phyloseq.extended' fails to install.
##' Compute incidence matrix of a tree
##'
##' @title incidenceMatrix
##' @param phy Required. A \code{phylo} class object
##' @return An incidence matrix M of size nedges(tree) x ntaxa(tree) where
##'         M[i,j] is set to 1 if taxa derives from edge i and 0 otherwise.
##' @note incidenceMatrix assumes that the tree is rooted. If not, it roots it
##'       it at an arbitrary taxa (taxa 1). 
incidenceMatrix <- function(phy) {
  if (!is.rooted(phy)) {
    warning("Tree is not rooted, incidence matrix may be meaningless")        
  }
  ## Construct incidence matrix of the tree (taxa x edge matrix)
  ## All taxa descending from an edge are set to 1, all others to -1
  ntaxa <- length(phy$tip.label)
  nedges <- nrow(phy$edge)
  incidence <- matrix (0,
                       nrow = ntaxa,
                       ncol = nedges,
                       dimnames = list(phy$tip.label, phy$edge[, 2]))
  ## Incidence of internal edges
  phy.part <- prop.part(phy) ## clade composition indexed by (shifted) node number
  for (i in 2:length(phy.part)) { ## first clade corresponds to root node
    edge <- which(phy$edge[, 2] == i + ntaxa) ## node numbers are shifted by ntaxa 
    incidence[phy.part[[i]] , edge] <- 1
  }
  ## Incidence of pendant edges
  ## pendant.edges[i] is the edge leading to tip i. 
  pendant.edges <- match(seq_len(ntaxa), phy$edge[ , 2])
  for (i in seq_len(ntaxa)) {
    incidence[i, pendant.edges[i]] <- 1
  }
  attr(incidence, "pendant.edges") <- pendant.edges
  return(incidence)
}
########################################
##############################

## Import the All-Species Living Tree Project tree file
tree_link <- "LTP_all_08_2023.ntree"
LTP_tree <- read.tree(tree_link) # ape

## Import the All-Species Living Tree Project alignment file
alignment_link <- "LTP_08_2023_aligned.fasta"
LTP_alignment <- read.alignment(alignment_link, format="fasta") # seqinr

## Remove Archaea
Bacteria_node <- which(LTP_tree$node.label == "Bacteria")
tree_no_archaea <- extract.clade(LTP_tree, node = (length(LTP_tree$tip.label) + Bacteria_node))

## Function to simulate alpha diversity
sim_alpha <- function(seed, N, LTP_tree){
  ## Set a seed number for reproducibility
  set.seed(seed)
  
  ## Create a new directory to save the result
  if (!dir.exists("alpha_outputs")) {
  dir.create("alpha_outputs")
    }
  dir_name <- paste0("alpha_outputs/output_", seed)
  dir.create(dir_name, showWarnings = FALSE)
  
  ###### Random sampling taxa from the phylogenetic tree ######
  ## Conduct random sampling of N tips from the tree
  sampled_tips <- sample(LTP_tree$tip.label, N)
  
  ## Create a tree containing only the sampled tips
  sample_tree <- drop.tip(LTP_tree, LTP_tree$tip.label[!LTP_tree$tip.label %in% sampled_tips]) # ape
  
  ## Extract GenBank ID of sampled taxa
  # Regex pattern to extract codes inside parentheses (focusing on the format similar to GenBank ID)
  pattern <- "\\(([^)]+)\\)"
  
  # Extract the codes and remove the parentheses from the extracted results
  GenBank <- str_extract(sampled_tips, pattern)
  GenBank <- gsub("[()]", "", GenBank)
  #############################################################
  
  ###### Simulate frequencies ######
  ###### referring to the Harrison et al. (2019)
  ## Set the number of OTUs
  notus <-length(sampled_tips)
  
  ## 1. All taxa are equally abundant (1/n)
  df_frequency_eq <- data.frame(otu = sampled_tips,
                                GenBank = GenBank,
                                freq = 1 / notus)
  
  ## Write dataframe to a CSV file
  freq_file <- paste0(N, "_df_frequency_eq.csv")
  write.csv(df_frequency_eq, file = file.path(dir_name, freq_file), row.names = FALSE)
  
  ## 2. Only a few taxa are abundant, some are intermediate, and most are rare
  ## Pareto (Shape = 4, location = 1)
  repeat {
    alphas_4 <- rpareto(n = notus, scale = 1, shape = 4)

    # Separate by abundance class
    abundant_4 <- which(alphas_4 >= 5)
    medium_4 <- which(alphas_4 < 5 & alphas_4 > 2)
    rare_4 <- which(alphas_4 <= 2)

    # Repeat the above steps until all abundance classes have at least one taxon
    if (length(abundant_4) > 0 & length(medium_4) > 0 & length(rare_4) > 0) {
      freq_4 <- alphas_4 / sum(alphas_4) # sum to 1 # Get frequencies proportional to alphas
      break
    }
  }
  
  ## Make a frequency data frame  
  df_frequency_4 <- data.frame(otu = sampled_tips,
                               GenBank = GenBank,
                               freq = as.vector(freq_4))
  
  ## Write dataframe to a CSV file
  freq_file <- paste0(N, "_df_frequency_4.csv")
  write.csv(df_frequency_eq, file = file.path(dir_name, freq_file), row.names = FALSE)
  
  ## 3. Only a few taxa are abundant, most of the taxa are rare
  ## Pareto (Shape = 0.7, location = 1)
  repeat {
    alphas_0.7 <- rpareto(notus, scale = 1, shape = 0.7)   
    
    # Separate by abundance class
    abundant_0.7 <- which(alphas_0.7 >= 1000) 
    medium_0.7 <- which(alphas_0.7 < 1000 & alphas_0.7 > 100)
    rare_0.7 <- which(alphas_0.7 <= 100) 
  
    # Repeat the above steps until all abundance classes have at least one taxon
    if(length(abundant_0.7) > 0 & length(medium_0.7) > 0 & length(rare_0.7) > 0) {
      freq_0.7 <- alphas_0.7 / sum(alphas_0.7) # sum to 1
      break
    }
  }
  
  ## Make a frequency data frame  
  df_frequency_0.7 <- data.frame(otu = sampled_tips,
                                 GenBank = GenBank,
                                 freq = as.vector(freq_0.7))

  ## Write dataframe to a CSV file
  freq_file <- paste0(N, "_df_frequency_0.7.csv")
  write.csv(df_frequency_eq, file = file.path(dir_name, freq_file), row.names = FALSE)
  ###################################
  
  ###### Calculating BWPD ######
  ## Function to calculate BWPD-theta
  ###### Code modified from 'phylodiv' function of 'mahendra-mariadassou/phyloseq-extended' ######
  calc_BWPD <- function(tree, freq_data, theta) {
    ## Convert the frequency matrix into the community table format of 'phyloseq'
    com_mat <- rbind(freq = freq_data$freq)
    colnames(com_mat) <- sampled_tips
    otu_mat <- otu_table(com_mat, taxa_are_rows = FALSE)
    
    ## Construct incidence matrix of the tree
    incidence <- incidenceMatrix(tree)
    
    ## Order incidence matrix according to community tables
    incidence <- incidence[colnames(otu_mat), ]
    
    ## Create community phylogeny matrix by multiplying (community x edge matrix)
    ## where cpm_{ij} gives the abundance of OTUs originating from branch j in community i.
    cpm <- otu_mat %*% incidence
    ## Convert to incidence matrix (0/1) and multiply by edge length to obtain PD per community.
    if (theta == 0) {
      cpm[cpm > 0] <- 1
    } else {
      cpm <- (2 * pmin(cpm, 1-cpm))^theta
    }
    BWPD <- cpm %*% tree$edge.length
    
    return(BWPD)
  }

  ## Calculate BWPD-theta (theta = 1)
  BWPD_eq <- calc_BWPD(sample_tree, df_frequency_eq, 1)  
  BWPD_4 <- calc_BWPD(sample_tree, df_frequency_4, 1)
  BWPD_0.7 <- calc_BWPD(sample_tree, df_frequency_0.7, 1)
  
  df_BWPD <- data.frame(frequency = c("equally abundant", "intermediate", "extremely skewed"),
                        BWPD = c(BWPD_eq, BWPD_4, BWPD_0.7))
  ##############################

  ###### Calculating hmaPD-alpha metrics ######
  ## Filter the alignment of sampled taxa from the entire alignment file
  filter_alignment <- function(LTP_alignment, GenBank) {
    nam <- LTP_alignment$nam[LTP_alignment$nam %in% GenBank]
    seq  <- LTP_alignment$seq[LTP_alignment$nam %in% GenBank]
    com <- LTP_alignment$com[LTP_alignment$nam %in% GenBank]
    nb <- length(nam)

    return(seqinr::as.alignment(nb = nb, nam = nam, seq = seq, com = com))
  }
  alignment <- filter_alignment(LTP_alignment, GenBank)
  
  ## Compute a pairwise distance matrix 
  dist_matrix <- as.matrix(dist.alignment(alignment, matrix = "identity", gap = TRUE))
                           
  ## Extract the species/strain names from the distance matrix
  otu <- rownames(dist_matrix)

  ## Create a matrix containing the extracted species/strain names
  ## to establish comparison names in the format of x---y
  names <- outer(otu, otu, function(x, y) paste(x, "---", y, sep = ""))

  ## Extract the comparison names from the 'names' matrix
  pair <- names[upper.tri(names)] %>% as.character()

  ## Extract pairwise distances from the distance matrix 
  sqrt_pairwise_distance <-dist_matrix[upper.tri(dist_matrix)] %>% as.numeric()
  pairwise_distance <- (sqrt_pairwise_distance)^2
  
  ## Calculate the relative heteroduplex mobilities referring to Li et al. (2021)
  dH <- (pairwise_distance - 0.011) / 0.372

  ## Create a data frame to store comparison names and pairwise distances
  df_distance <- data.frame(pair = pair,
                            pairwise_distance = pairwise_distance,
                            dH = dH)

  ## Initialise a dataframe to store the HMAa results
  df_HMAa <- data.frame(frequency = character(), HMAa = numeric(), stringsAsFactors = FALSE)

  ## Save frequency dataframes to a list
  df_frequency_list <- list(
    "equally abundant" = df_frequency_eq,
    "intermediate" = df_frequency_4,
    "extremely skewed" = df_frequency_0.7
    )

  ## Calculate HMAa for each df_frequency
  for (freq_name in names(df_frequency_list)) {
    df_frequency <- df_frequency_list[[freq_name]]
  
    ## Calculate the pairwise distance multiplied by frequencies for each pair
    df_distance$distXfreq <- sapply(df_distance$pair, function(pair) {
      ## Split the pair into two OTUs
      otus <- strsplit(pair, "---") [[1]]

      ## Extract frequencies of corresponding OTUs from the 'df_frequency' data frame
      freq_1 <- df_frequency$freq[df_frequency$GenBank == otus[1]]
      freq_2 <- df_frequency$freq[df_frequency$GenBank == otus[2]]
    
      ## Calculate the relative heteroduplex mobility multiplied by frequencies
      distXfreq <- freq_1 * freq_2 * df_distance$dH[df_distance$pair == pair]
    
      return(distXfreq)
    })
  
    ## Calculate the hmaPD-alpha value by summing up  
    HMAa <- 2 * sum(df_distance$distXfreq)

    ## Add to the result dataframe
    df_HMAa <- rbind(df_HMAa, data.frame(frequency = freq_name, HMAa = HMAa))
    }
  #############################################

  ## Create a local dataframe for the results
  local_result <- data.frame(N = rep(N, nrow(df_HMAa)),
                             frequency = df_HMAa$frequency,
                             BWPD = df_BWPD$BWPD,
                             HMAa = df_HMAa$HMAa,
                             stringsAsFactors = FALSE)

  return(local_result)
}

## Define seed numbers and N values
seeds <- 1001:1030
N_values <- c(10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 180, 200)  # Add more values as needed

## Initialise an empty dataframe to store results from all seeds and N values
df_all_results <- data.frame(N = integer(),
                             frequency = character(),
                             BWPD = numeric(),
                             HMAa = numeric(),
                             seed = integer(),
                             stringsAsFactors = FALSE)

## Run simulation
for (seed in seeds) {
  for (N in N_values) {
    # Perform simulation for the current N
    result <- sim_alpha(seed, N, LTP_tree) 
    result$seed <- seed  # Add an identifier for the current seed to the result
    
    # Append results to the main dataframe
    df_all_results <- rbind(df_all_results, result)
  }
}

## Save the combined results to a CSV file
write.csv(df_all_results, file = file.path("alpha_outputs", "alpha_simulation_results.csv"), row.names = FALSE)
