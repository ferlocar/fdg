#setwd("C:/Users/ferlo/OneDrive/Documents/NYU/Foster (Research)/network data")

# Load Node features
feat <- read.csv("107.feat", header = FALSE, sep = " ")
colnames(feat) <- gsub("V", "F", colnames(feat))
colnames(feat) <- c("ID", colnames(feat)[1:ncol(feat)-1])
feat$ID <- paste0("N", feat$ID)

# Load edges
edges <- read.csv("107.edges", header = FALSE, sep = " ")
edges$V1 <- paste0("N", edges$V1)
edges$V2 <- paste0("N", edges$V2)

# Set up adjacency matrix
library(Matrix)
adj_mat <- Matrix(data = 0, ncol=nrow(feat), nrow = nrow(feat), sparse = TRUE)
colnames(adj_mat) <- feat$ID
rownames(adj_mat) <- feat$ID
adj_mat[cbind(match(edges$V1, rownames(adj_mat)), match(edges$V2, colnames(adj_mat)))] <- 1

# Set up igraph
library(igraph)
net <- graph_from_adjacency_matrix(adj_mat, mode="undirected")
# Evaluate force-directed graph method
source("fdg.R")
l <- get_l(net)
results <- matrix(0, ncol=2, nrow=ncol(feat)-1)
rownames(results) <- colnames(feat)[-1]

# Test different features as dependent variables
for (vn in rownames(results)){
  print(vn)
  net <- set_vertex_attr(net, "y", index=V(net), feat[[vn]])
  results[vn,] <- eval_alg(net, l, adj_mat)
}