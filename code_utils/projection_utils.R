
#' Get a knn matrix given SVD matrix for reference and query
#' 
#' @param data cell x dims matrix for reference dataset
#' @param query cell x dims matrix for query dataset (optional)
#' @param k number of neighbors to calculate
#' @param threads Number of threads to use. Note that result is non-deterministic
#'          if threads > 1
#' @param ef ef parameter for RccppHNSW::hnsw_search. Increase for slower search but
#'          improved accuracy
#' @param verbose whether to print progress information during search
#' @return List of 2 matrices -- idx for cell x K neighbor indices,
#'         dist for cell x K neighbor distances.
#'         If no query is given, nearest neighbors are found mapping
#'         the data matrix to itself, prohibiting self-neighbors
get_knn <- function(data, query=NULL, k = 10, metric="euclidean", verbose=TRUE, threads=1, ef=100) {
  index <- RcppHNSW::hnsw_build(
    data, 
    distance=metric, 
    verbose=verbose,
    n_threads=threads
  )
  if (is.null(query)) {
      res <- RcppHNSW::hnsw_search(
        data,
        index,
        k+1,
        ef=ef,
        verbose=verbose,
        n_threads = threads
      )
      missed_searches <- sum(res$idx[,1] != seq_len(nrow(data)))
      if (any(res$idx[,1] != seq_len(nrow(data)))) {
        stop(sprintf("KNN search didn't find self-neighbor for %d datapoints", missed_searches))
      }
      return(list(
        idx = res$idx[,-1,drop=FALSE],
        dist = res$dist[,-1,drop=FALSE]
      ))
  } else {
      res <- RcppHNSW::hnsw_search(
        query,
        index,
        k,
        ef=ef,
        verbose=verbose,
        n_threads = threads
      )
      return(res)
  }
}


#' Internal worker function for label transfer and quality estimation
#' @param knn KNN object as returned from `get_knn`.
#' @param labels Vector of label values to transfer. Should be type integer with minimum >= 1
#'               (convert from strings or factors before calling this function)
#' @return Matrix of cells x labels
get_label_scores <- function(knn, labels) {
  stopifnot(is.integer(labels))
  stopifnot(min(labels) >= 1)

  nearest_labels <- matrix(labels[knn$idx], nrow=nrow(knn$idx))
  
  weights <- 1 - knn$dist / knn$dist[,ncol(knn$dist)]
  weights <- weights / rowSums(weights)
 
  label_scores <- map(
    seq_len(max(labels)),
    ~ rowSums(weights * (nearest_labels == .x))                    
  ) %>% do.call(cbind, .)
  
  return(label_scores)
}

#' Transfer labels using the most common value from the K nearest neighbors
#' Neighbors are weighted similarly to Stuart et al. (2019), by 1 - (dist_i/dist_k)
#' @param knn KNN object as returned from `get_knn`.
#' @param labels Vector of label values to transfer
transfer_labels <- function(knn, labels) {
  label_values <- unique(labels)
  label_indices <- match(labels, label_values)
  
  label_scores <- get_label_scores(knn, label_indices)
  best_label <- max.col(label_scores, ties.method = "first")
  
  label_values[best_label]
}