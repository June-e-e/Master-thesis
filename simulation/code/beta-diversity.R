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

# Import the All-Species Living Tree Project tree file
tree_link <- "LTP_all_08_2023.ntree"
LTP_tree <- read.tree(tree_link)

## Import the All-Species Living Tree Project alignment file
alignment_link <- "LTP_08_2023_aligned.fasta"
LTP_alignment <- read.alignment(alignment_link, format="fasta")

## Remove Archaea
Bacteria_node <- which(LTP_tree$node.label == "Bacteria")
tree_no_archaea <- extract.clade(LTP_tree, node = (length(LTP_tree$tip.label) + Bacteria_node))

## Function to simulate beta diversity
sim_beta <- function(seed, M_N_pairs, tree) {
  ## Set a seed number for reproducibility
  set.seed(seed)
  
  ## Create a new directory to save the result
  if (!dir.exists("beta_outputs")) {
    dir.create("beta_outputs")
  }
  dir_name <- paste0("beta_outputs/output_", seed)
  dir.create(dir_name, showWarnings = FALSE)

  for (pair in M_N_pairs) {
    M <- pair[1]
    N <- pair[2]
    
    ###### Random sampling taxa from the phylogenetic tree ######
    comA <- sample(tree$tip.label, M)
    comB <- sample(tree$tip.label, N)
    comTotal <- union(comA, comB)
    
    # Create a tree containing only the sampled tips
    tree_comA <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% comA])
    tree_comB <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% comB])
    tree_comTotal <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% comTotal])
    
    # Extract GenBank ID of sampled taxa
    pattern <- "\\(([^)]+)\\)"
    
    GenBank_A <- str_extract(comA, pattern)
    GenBank_A <- gsub("[()]", "", GenBank_A)
    
    GenBank_B <- str_extract(comB, pattern)
    GenBank_B <- gsub("[()]", "", GenBank_B)
    
    GenBank_Total <- str_extract(comTotal, pattern)
    GenBank_Total <- gsub("[()]", "", GenBank_Total)
    
    ## Check for shared otus
    shared_otu <- length(intersect(GenBank_A, GenBank_B))
    #############################################################
    
    ###### Simulate frequencies using a Pareto distribution (Shape = 0.7, location = 1) ######
    ###### Only a few taxa are abundant, most of the taxa are rare
    ###### referring to the Harrison et al. (2019)
    freq_A <- NULL
    freq_B <- NULL
    
    repeat {
      abun_A <- rpareto(n = M, scale = 1, shape = 0.7)   
      abundant_A <- which(abun_A >= 1000) 
      medium_A <- which(abun_A < 1000 & abun_A > 100)
      rare_A <- which(abun_A <= 100) 
      
      if(length(abundant_A) > 0 & length(medium_A) > 0 & length(rare_A) > 0) {
        freq_A <- abun_A / sum(abun_A) # sum to 1
        break
      }
    }
    
    repeat {
      abun_B <- rpareto(n = N, scale = 1, shape = 0.7)   
      abundant_B <- which(abun_B >= 1000) 
      medium_B <- which(abun_B < 1000 & abun_B > 100)
      rare_B <- which(abun_B <= 100) 
      
      if(length(abundant_B) > 0 & length(medium_B) > 0 & length(rare_B) > 0) {
        freq_B <- abun_B / sum(abun_B) # sum to 1
        break
      }
    }
    
    # Make a frequency data frame for communities A and B
    df_freq_A <- data.frame(otu = comA,
                            GenBank = GenBank_A,
                            freq = as.numeric(freq_A))
    
    df_freq_B <- data.frame(otu = comB,
                            GenBank = GenBank_B,
                            freq = as.numeric(freq_B))
    
    ## Combine two frequency data frames for the combined community
    df_freq_Total <- rbind(
      transform(df_freq_A, freq = freq / 2),
      transform(df_freq_B, freq = freq / 2)
    )
    
    ## Remove duplicates when a shared taxa is present in two communities
    ## Sum the frequencies for duplicates
    df_freq_Total <- df_freq_Total %>%
      group_by(otu, GenBank) %>%
      summarise(freq = sum(freq), .groups = 'drop')
    
    freq_fileA <- paste0(M, "_", N, "_df_frequency_A.csv")
    write.csv(df_freq_A, file = file.path(dir_name, freq_fileA), row.names = FALSE)
    freq_fileB <- paste0(M, "_", N, "_df_frequency_B.csv")
    write.csv(df_freq_B, file = file.path(dir_name, freq_fileB), row.names = FALSE)
    freq_fileTotal <- paste0(M, "_", N, "_df_frequency_Total.csv")
    write.csv(df_freq_Total, file = file.path(dir_name, freq_fileTotal), row.names = FALSE)
    ###################################
    
    ###### Calculating hmaPD-alpha ######
    ## function to calculate HMAa 
    calc_HMAa <- function(data) {
      ## Filter the alignment of sampled taxa from the entire alignment file
      nam <- LTP_alignment$nam[LTP_alignment$nam %in% data$GenBank]
      seq  <- LTP_alignment$seq[LTP_alignment$nam %in% data$GenBank]
      com <- LTP_alignment$com[LTP_alignment$nam %in% data$GenBank]
      nb <- length(nam)
      alignment <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq, com = com)
      
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
      
      ## Calculate the pairwise distance multiplied by frequencies for each pair
      df_distance$distXfreq <- sapply(df_distance$pair, function(pair) {
        ## Split the pair into two OTUs
        otus <- strsplit(pair, "---") [[1]]
        
        ## Extract frequencies of corresponding OTUs from the 'df_frequency' data frame
        freq_1 <- data$freq[data$GenBank == otus[1]]
        freq_2 <- data$freq[data$GenBank == otus[2]]
        
        ## Calculate the relative heteroduplex mobility multiplied by frequencies
        distXfreq <- freq_1 * freq_2 * df_distance$dH[df_distance$pair == pair]
        
        return(distXfreq)
      })
      
      ## Calculate the hmaPD-alpha value by summing up  
      HMAa <- 2 * sum(df_distance$distXfreq)
      
      return(HMAa)
    }
    
    ## Calculate alpha diversity of each community  
    comA_HMAa <- calc_HMAa(df_freq_A)
    comB_HMAa <- calc_HMAa(df_freq_B)
    #############################
    
    ###### Calculating weighted UniFrac ######
    ## Create a community table
    # Initialise a matrix for a community table
    com_mat <- matrix(0, nrow = 2, ncol = length(comTotal))
    
    # Assign column names and row names to the community table
    colnames(com_mat) <- comTotal
    rownames(com_mat) <- c("comA", "comB")
    
    # Fill in the frequencies
    com_mat["comA", comA] <- df_freq_A$freq
    com_mat["comB", comB] <- df_freq_B$freq
    
    ## Create phyloseq-class object
    physeq <- phyloseq(otu_table(com_mat, taxa_are_rows = FALSE), tree_comTotal)
    
    ## Calcualte weighted UniFrac
    UniFrac <- phyloseq::UniFrac(physeq, weighted = TRUE)
    ##########################################
    
    ###### Calculating hmaPD-beta ######
    ## Calculate alpha diversity of the combined communities
    All_HMAa <- calc_HMAa(df_freq_Total)
    
    ## Calculate the hmaPD-beta 
    HMAb <- 2 * All_HMAa - (comA_HMAa + comB_HMAa)
    ####################################
    
    ###### Calculating Rao's dissimilarity coefficient ######
    ## Function to calculate between-community diversity Hij
    calc_Hij <- function(data1, data2, data_all){
      ## Filter the alignment of sampled taxa from the entire alignment file
      nam <- LTP_alignment$nam[LTP_alignment$nam %in% data_all$GenBank]
      seq  <- LTP_alignment$seq[LTP_alignment$nam %in% data_all$GenBank]
      com <- LTP_alignment$com[LTP_alignment$nam %in% data_all$GenBank]
      nb <- length(nam)
      alignment <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq, com = com)
      
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
      
      ## Abundance matrix for between-community diversity Hij
      abundance_Hij <- outer(data1$freq, data2$freq, FUN = function(x, y) x * y)
      rownames(abundance_Hij) <- data1$GenBank
      colnames(abundance_Hij) <- data2$GenBank
      
      # Initialize distance_values matrix
      dist_matrix_Hij <- matrix(0, nrow = length(data1$GenBank), ncol = length(data2$GenBank))
      
      # Set row and column names for distance_values matrix
      rownames(dist_matrix_Hij) <- data1$GenBank
      colnames(dist_matrix_Hij) <- data2$GenBank
      
      # Match and extract distances
      for (i in seq_along(data1$GenBank)) {
        for (j in seq_along(data2$GenBank)) {
          # Find the corresponding distance in the dist_matrix
          dist_matrix_Hij[i, j] <- dist_matrix[match(data1$GenBank[i], otu), match(data2$GenBank[j], otu)]
        }
      }
      
      # Calculate the element-wise product of abundance_Hij and the distance values
      Hij_matrix <- abundance_Hij * dist_matrix_Hij
      Hij <- sum(Hij_matrix)
      
      return(Hij)
    }
    
    ## Calculate Rao's dissimilarity coefficient
    Hij <- calc_Hij(df_freq_A, df_freq_B, df_freq_Total)
    Rao_beta <- Hij - 0.5 * (comA_HMAa + comB_HMAa)
    
    #########################################################
    
    size_difference <- abs(M - N) # Calculate the absolute difference between M and N
    
    ## Create a local dataframe for the results
    local_result <- data.frame(size_difference = size_difference,
                               shared_otu = shared_otu,
                               comA_HMAa = comA_HMAa,
                               comB_HMAa = comB_HMAa,
                               All_HMAa = All_HMAa,
                               HMAb = HMAb,
                               Rao_beta = Rao_beta,
                               UniFrac = as.numeric(UniFrac),
                               seed = seed,
                               stringsAsFactors = FALSE)
    
  return(local_result)
  }
}

## Define seed numbers and (M, N) pairs
seeds <- 1001:1030
M_N_pairs <- list(c(20, 20), c(20, 40), c(20, 60), c(20, 80),
                  c(40, 40), c(40, 60), c(40, 80), c(40, 100),
                  c(60, 60), c(60, 80), c(60, 100), c(60, 140),
                  c(80, 80), c(80, 100), c(80, 120), c(80, 160))

## Initialise an empty dataframe to store results from all seeds and (M, N) pairs
df_all_results <- data.frame(size_difference = integer(),
                             shared_otu = integer(),
                             comA_HMAa = numeric(),
                             comB_HMAa = numeric(),
                             All_HMAa = numeric(),
                             HMAb = numeric(),
                             Rao_beta = numeric(),
                             UniFrac = numeric(),
                             seed = integer(),
                             stringsAsFactors = FALSE)

## Run simulation
for (seed in seeds) {
  for (pair in M_N_pairs) {
  result <- sim_beta(seed, list(pair), tree_no_archaea) 
  df_all_results <- rbind(df_all_results, result)
  }
}

## Save the combined results to a CSV file
write.csv(df_all_results, file = file.path("beta_outputs", "beta_simulation_results.csv"), row.names = FALSE)

