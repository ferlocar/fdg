library(igraph)
library(igraphdata)
library(e1071)
library(AUC)

get_l <- function(net, plot = FALSE){
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

eval_alg <- function(net, l, adj_mat){
  # Set up cross-validation
  data <- data.frame(y=as.numeric(V(net)$y),l1=l[,1],l2=l[,2], row.names = rownames(adj_mat))
  n_folds <- min(10, sum(data$y))
  if(n_folds < 2){
    print("Not enough positive observations to perform cross validation")
    return (c(0, 0))
  }
  folds <- get_strat_folds(data, n_folds)
  avg_auc <- 0
  bench_auc <- 0
  #Perform Cross Validation
  for(i in 1:n_folds){
    print(paste0("Fold ", i))
    fold_ixs <- which(folds==i,arr.ind=TRUE)
    train_data <- data[-fold_ixs,]
    test_data <- data[fold_ixs,]
    # Force-directed graph
    model <- svm(y~l1+l2, train_data, type="C",kernel='radial', probability=TRUE)
    y_pred <- predict(model, test_data, probability=TRUE)
    scores <- attr(y_pred, "probabilities")[,2]
    avg_auc <- avg_auc + auc(roc(scores, factor(test_data$y)))
    # Count edges to positives (benchmark)
    bench_scores <- (data$y[-fold_ixs] %*% adj_mat[-fold_ixs, fold_ixs])[1,]
    bench_auc <- bench_auc + auc(roc(bench_scores, factor(test_data$y)))
  }
  avg_auc <- avg_auc/n_folds
  bench_auc <- bench_auc/n_folds
  print("Force-Directed Graph:")
  print(paste0("AVG. AUC (", n_folds, " folds): ", avg_auc))
  print("Counting edges:")
  print(paste0("Benchmark AUC (", n_folds, " folds): ", bench_auc))
  return(c(avg_auc, bench_auc))
}

# Karate example
#data(karate)
#net <- karate
#net <- set_vertex_attr(net, "y", index=V(net), vertex_attr(net, "Faction") - 1)