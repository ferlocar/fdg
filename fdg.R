library(igraph)
library(igraphdata)
library(e1071)
library(AUC)

get_coords <- function(net, plot = FALSE){
  l <- layout_with_fr(net)
  if(plot){
    plot(net, edge.arrow.mode = 0,layout=l)
  }
  # Normalize layout
  l[,1] = (l[,1]-mean(l[,1]))/sd(l[,1])
  l[,2] = (l[,2]-mean(l[,2]))/sd(l[,2])
  return(l)
}

get_strat_folds <- function(data, n_folds){
  pos_ixs <- which(data$y==1,arr.ind=TRUE)
  neg_ixs <- which(data$y==0,arr.ind=TRUE)
  pos_folds <- cut(1:length(pos_ixs), breaks=n_folds, labels=FALSE)
  neg_folds <- cut(1:length(neg_ixs), breaks=n_folds, labels=FALSE)
  folds <- numeric(nrow(data))
  folds[pos_ixs] <- pos_folds[sample(length(pos_ixs))]
  folds[neg_ixs] <- neg_folds[sample(length(neg_ixs))]
  return(folds)
}

eval_alg <- function(y, l, adj_mat, paths){
  # Set up cross-validation
  data <- data.frame(y=y,l1=l[,1],l2=l[,2], row.names = rownames(adj_mat))
  n_folds <- min(10, sum(data$y))
  if(n_folds < 2){
    print("Not enough positive observations to perform cross validation")
    return (c(0, 0))
  }
  folds <- get_strat_folds(data, n_folds)
  avg_auc <- 0
  bench_auc <- 0
  both_auc <-0
  floyd_auc <-0
  #Perform Cross Validation
  for(i in 1:n_folds){
    print(paste0("Fold ", i))
    test_ixs <- which(folds==i,arr.ind=TRUE)
    train_data <- data[-test_ixs,]
    test_data <- data[test_ixs,]
    # Use Floyd's stuff
    f_scores <- rowSums(paths[,(data$y>0)[-test_ixs], drop=FALSE]) * -1
    # The min_val is to do something when none of the nodes can reach a "positive node"
    min_val <- -1024*1024
    f_scores[is.na(f_scores)] <- min(min_val, f_scores, na.rm=TRUE)
    floyd_auc <- floyd_auc + auc(roc(f_scores[test_ixs], factor(test_data$y)))
    # Count edges to positives (benchmark)
    counts <- (data$y[-test_ixs] %*% adj_mat[-test_ixs,])[1,]
    bench_auc <- bench_auc + auc(roc(counts[test_ixs], factor(test_data$y)))
    # Add counts as part of the data
    test_data$counts <- counts[test_ixs]
    train_data$counts <- counts[-test_ixs]
    # Force-directed graph
    model <- svm(y~l1+l2, train_data, type="C",kernel='radial', probability=TRUE)
    y_pred <- predict(model, test_data, probability=TRUE)
    scores <- attr(y_pred, "probabilities")[,2]
    avg_auc <- avg_auc + auc(roc(scores, factor(test_data$y)))
    # Both 
    model <- svm(y~l1+l2+counts, train_data, type="C",kernel='radial', probability=TRUE)
    y_pred <- predict(model, test_data, probability=TRUE)
    scores <- attr(y_pred, "probabilities")[,2]
    both_auc <- both_auc + auc(roc(scores, factor(test_data$y)))
  }
  floyd_auc <- floyd_auc/n_folds
  avg_auc <- avg_auc/n_folds
  bench_auc <- bench_auc/n_folds
  both_auc <- both_auc/n_folds
  print("Force-Directed Graph:")
  print(paste0("AVG. AUC (", n_folds, " folds): ", avg_auc))
  print("Counting edges:")
  print(paste0("Benchmark AUC (", n_folds, " folds): ", bench_auc))
  print("Edges and FD Graph:")
  print(paste0("Both AUC (", n_folds, " folds): ", both_auc))
  print("Floyd Algorithm:")
  print(paste0("Floyd AUC (", n_folds, " folds): ", floyd_auc))
  return(c(avg_auc, bench_auc, both_auc, floyd_auc))
}

# Karate example
#data(karate)
#net <- karate
#net <- set_vertex_attr(net, "y", index=V(net), vertex_attr(net, "Faction") - 1)