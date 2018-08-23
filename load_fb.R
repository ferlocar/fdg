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
get_adj_mat <- function(default_val){
  adj_mat <- Matrix(data = default_val, ncol=nrow(feat), nrow = nrow(feat), sparse = TRUE)
  colnames(adj_mat) <- feat$ID
  rownames(adj_mat) <- feat$ID
  adj_mat[cbind(match(edges$V1, rownames(adj_mat)), match(edges$V2, colnames(adj_mat)))] <- 1
  return(adj_mat)
}


# Use Floyd-Warshall algorithm
library(e1071)
paths <- allShortestPaths(get_adj_mat(NaN))$length


# Set up igraph
library(igraph)
net <- graph_from_adjacency_matrix(get_adj_mat(0), mode="undirected")
# Evaluate force-directed graph method
source("fdg.R")
coords <- get_coords(net)
results <- matrix(0, ncol=4, nrow=ncol(feat)-1)
rownames(results) <- colnames(feat)[-1]

# Test different features as dependent variables
set.seed(0)
for (vn in rownames(results)){
  print(vn)
  net <- set_vertex_attr(net, "y", index=V(net), feat[[vn]])
  target_var <- as.numeric(V(net)$y)
  results[vn,] <- eval_alg(target_var, coords, get_adj_mat(0), paths)
}

colnames(results) <- c("fdg_auc", "bench_auc", "both_auc", "floyd_auc")
write.csv(results, file = "results.csv")