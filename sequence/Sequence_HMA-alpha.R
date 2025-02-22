## Load packages
library(tidyverse)
library(seqinr)
library(ape)
library(phyloseq)
library(phangorn)
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
## set seed
set.seed(1000)

## Import files
alignment <- read.alignment("consensus_alignment/MW1_consensus_alignment.fasta", format = "fasta")
abundance <- read_tsv("EMU_output/MW1/MW1_rel-abundance_nondoubletons.tsv")

## Create a community table
com_mat <- matrix(0, nrow = 1, ncol = length(abundance$tax_id),
                  dimnames = list("sample", abundance$tax_id))

## Fill in the frequencies
com_mat["sample", ] <- abundance$abundance

## Compute a pairwise distance matrix 
dist_matrix <- as.matrix(dist.alignment(alignment, matrix = "identity", gap = TRUE))

## Construct the Neighbour-joining tree
nj_tree <- nj(dist_matrix)

## Change the tip labels to only include the tax_id
nj_tree$tip.label <- gsub(".*tax_id=(\\d+).*", "\\1", nj_tree$tip.label)

## Set the root using the midpoint method
rooted_tree <- midpoint(nj_tree)

## Visualise the tree
plot(rooted_tree)

## Function to calculate BWPD-theta
###### Code modified from 'phylodiv' function of 'mahendra-mariadassou/phyloseq-extended' ######
calc_BWPD <- function(tree, com_mat, theta) {
  ## Create phyloseq-class object
  otu_mat <- otu_table(com_mat, taxa_are_rows = FALSE)
  physeq <- phyloseq(otu_mat, tree)
  
  ## Construct incidence matrix of the tree
  incidence <- incidenceMatrix(tree)
  
  ## Order incidence matrix according to community tables
  incidence <- incidence[colnames(otu_mat), ]
  
  ## Create community phylogeny matrix by multiplying (community x edge matrix)
  ## where cpm_{ij} gives the abundance of OTUs originating from branch j in community i.
  cpm <- otu_mat %*% incidence
  
  ## Convert to incidence matrix (0/1) and multiply by edge length to obtain PD per community.
  cpm <- (2 * pmin(cpm, 1-cpm)) ^ theta
  
  BWPD <- cpm %*% tree$edge.length
  
  return(BWPD)
}

## Function to calculate HMAa
calc_HMAa <- function(dist_matrix, freq_df){
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
      frxfr = prod(freq_df$abundance[freq_df$tax_id %in% as.numeric(tax_ids)])
    ) %>%
    ungroup() %>%
    select(-tax_ids)  
 
  ## Calculate HMA-alpha
  hma_alpha <- 2 * sum(df_pd$frxfr * df_pd$pairwise_distance) 
  
  return(hma_alpha)
}

BWPD <- calc_BWPD(rooted_tree, com_mat, 1)
hma_alpha <- calc_HMAa(dist_matrix, abundance)

BWPD
hma_alpha

