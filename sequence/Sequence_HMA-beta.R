## Load packages
library(tidyverse)
library(Biostrings)
library(seqinr)
library(ape)
library(phyloseq)
library(phangorn)

## set seed
set.seed(1000)

## Import files
alignment1 <- read.alignment("consensus_alignment/W4_consensus_alignment.fasta", format = "fasta")
alignment2 <- read.alignment("consensus_alignment/B3_consensus_alignment.fasta", format = "fasta")
alignment_all <- read.alignment("Beta_alignments/Beta21_alignment.fasta", format = "fasta")
abundance1 <- read_tsv("EMU_output/W4/W4_rel-abundance_nondoubletons.tsv")
abundance2 <- read_tsv("EMU_output/B3/B3_rel-abundance_nondoubletons.tsv")

## Abundance matrix for between-community diversity Hij
abundance_Hij <- outer(abundance1$abundance, abundance2$abundance, FUN = function(x, y) x * y)
rownames(abundance_Hij) <- abundance1$tax_id
colnames(abundance_Hij) <- abundance2$tax_id

## Check for duplicated sequences
any(duplicated(unlist(alignment_all$seq)))

if (any(duplicated(unlist(alignment_all$seq)))) {
  # Identify sequence headers to remove
  nam_rm <- alignment_all$nam[duplicated(unlist(alignment_all$seq))]
  
  nam <- alignment_all$nam[!alignment_all$nam %in% nam_rm]
  seq <- alignment_all$seq[!alignment_all$nam %in% nam_rm]
  com <- NA
  nb <- length(nam)
  
  alignment_all <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq, com = com)
} 

## Check for common taxa
taxa1 <- gsub(".*tax_id=(\\d+).*", "\\1", alignment1$nam)
taxa2 <- gsub(".*tax_id=(\\d+).*", "\\1", alignment2$nam)
common_taxa <- intersect(taxa1, taxa2)
rm_taxa <-alignment2$nam[taxa2 %in% common_taxa]
rm_taxa <- c(rm_taxa, paste0("_R_", rm_taxa), sub("^_R_", "", rm_taxa))

## Remove common taxa from the alignment file 
if (length(common_taxa) > 0) {
  nam <- alignment_all$nam[!alignment_all$nam %in% rm_taxa]
  seq <- alignment_all$seq[!alignment_all$nam %in% rm_taxa]
  com <- NA
  nb <- length(nam)
  
  alignment_all <- seqinr::as.alignment(nb = nb, nam = nam, seq = seq, com = com)
} 

## Combine the abundances
abundance_all <- bind_rows(abundance1, abundance2) %>%
  group_by(tax_id) %>%  
  summarise(
    abundance = sum(abundance) / 2
  )

## Create a community table
# Initialise a matrix for a community table
com_mat <- matrix(0, nrow = 2, ncol = length(abundance_all$tax_id))

# Assign column names and row names to the community table
colnames(com_mat) <- abundance_all$tax_id
rownames(com_mat) <- c("com1", "com2")

# Fill in the frequencies
com_mat["com1", match(abundance1$tax_id, colnames(com_mat))] <- abundance1$abundance
com_mat["com2", match(abundance2$tax_id, colnames(com_mat))] <- abundance2$abundance 

## Compute a pairwise distance matrix 
dist_matrix <- as.matrix(dist.alignment(alignment_all, matrix = "identity", gap = TRUE))

## Construct the Neighbour-joining tree
nj_tree <- nj(dist_matrix)

## Change the tip labels to only include the tax_id
nj_tree$tip.label <- gsub(".*tax_id=(\\d+).*", "\\1", nj_tree$tip.label)

## Set the root using the midpoint method
rooted_tree <- midpoint(nj_tree)

## Visualise the tree
plot(rooted_tree)

## Create phyloseq-class object
physeq <- phyloseq(otu_table(com_mat, taxa_are_rows = FALSE), rooted_tree)

## Calcualte weighted UniFrac
UniFrac <- phyloseq::UniFrac(physeq, weighted = TRUE, normalized = TRUE)

## Function to calculate between-community diversity Hij
calc_Hij <- function(alignment, abundance){
  ## Compute a pairwise distance matrix 
  dist_matrix <- as.matrix(dist.alignment(alignment, matrix = "identity", gap = TRUE))
  
  ## Extract the species/strain names from the distance matrix
  otu <- rownames(dist_matrix) %>%
    gsub(".*tax_id=(\\d+).*", "\\1", .)
  
  ## Create a matrix containing the extracted species/strain names
  ## to establish comparison names in the format of x---y
  names <- outer(otu, otu, function(x, y) paste(x, "---", y, sep = ""))
  
  ## Extract the comparison names from the 'names' matrix
  pair <- names[upper.tri(names)] %>% as.character()
  
  ## Extract pairwise distances from the distance matrix 
  sqrt_pairwise_distance <- dist_matrix[upper.tri(dist_matrix)] %>% as.numeric()
  pairwise_distance <- (sqrt_pairwise_distance)^2
  
  ## Create a data frame to store comparison names and pairwise distances
  df_pd <- data.frame(pair = pair,
                      pairwise_distance = pairwise_distance)
  
  # Get row and column taxa
  row_taxa <- rownames(abundance)
  col_taxa <- colnames(abundance)
  
  # Initialize distance_values matrix
  dist_matrix_Hij <- matrix(0, nrow = length(row_taxa), ncol = length(col_taxa))
  
  # Set row and column names for distance_values matrix
  rownames(dist_matrix_Hij) <- row_taxa
  colnames(dist_matrix_Hij) <- col_taxa
  
  # Match and extract distances
  for (i in seq_along(row_taxa)) {
    for (j in seq_along(col_taxa)) {
      # Find the corresponding distance in the dist_matrix
      dist_matrix_Hij[i, j] <- dist_matrix[match(row_taxa[i], otu), match(col_taxa[j], otu)]
    }
  }
  
  # Calculate the element-wise product of abundance_Hij and the distance values
  Hij_matrix <- abundance * dist_matrix_Hij
  Hij <- sum(Hij_matrix)
  
  return(Hij)
}

## Function to calculate HMA-alpha diversity 
calc_alpha<- function(alignment, abundance){
  ## Compute a pairwise distance matrix 
  dist_matrix <- as.matrix(dist.alignment(alignment, matrix = "identity", gap = TRUE))
  
  ## Extract the species/strain names from the distance matrix
  otu <- rownames(dist_matrix) %>%
    gsub(".*tax_id=(\\d+).*", "\\1", .)
  
  ## Create a matrix containing the extracted species/strain names
  ## to establish comparison names in the format of x---y
  names <- outer(otu, otu, function(x, y) paste(x, "---", y, sep = ""))
  
  ## Extract the comparison names from the 'names' matrix
  pair <- names[upper.tri(names)] %>% as.character()
  
  ## Extract pairwise distances from the distance matrix 
  sqrt_pairwise_distance <- dist_matrix[upper.tri(dist_matrix)] %>% as.numeric()
  pairwise_distance <- (sqrt_pairwise_distance)^2
  
  ## Create a data frame to store comparison names and pairwise distances
  df_pd <- data.frame(pair = pair,
                      pairwise_distance = pairwise_distance)
  
  ## Add the product of abundances in each pair to the data frame
  df_pd <- df_pd %>%
    rowwise() %>%
    mutate(
      tax_ids = list(strsplit(pair, "---")[[1]]),
      frxfr = prod(abundance$abundance[abundance$tax_id %in% tax_ids])
    ) %>%
    ungroup() %>%
    select(-tax_ids)  
  
  ## Calculate HMA-alpha
  hma_alpha <- 2 * sum(df_pd$frxfr * df_pd$pairwise_distance) 
  
  return(hma_alpha)
}

## Calculate HMA-beta diversity
alpha_all <- calc_alpha(alignment_all, abundance_all)
alpha1 <- calc_alpha(alignment1, abundance1)
alpha2 <- calc_alpha(alignment2, abundance2)

hma_beta <- 2 * alpha_all - (alpha1 + alpha2)

## Calculate Rao's dissimilarity coefficient
Hij <- calc_Hij(alignment_all, abundance_Hij)

Rao_beta <- Hij - 0.5 * (alpha1 + alpha2)

UniFrac
Rao_beta
hma_beta


