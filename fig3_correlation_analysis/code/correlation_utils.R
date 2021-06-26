
#### Correlations ####
#' Pearson or spearman correlation on sparse matrix
#' Operates like cor function, finding all pairwise correlations
sparse.cor <- function(x, y, method=c("pearson", "spearman")) {
    stopifnot(nrow(x) == nrow(y))
    method <- match.arg(method)
    if (method == "spearman") {
        if (is(x, "dgCMatrix")) x <- column_rank(x)
        else                    x <- apply(x, 2, rank)
        
        if (is(y, "dgCMatrix")) y <- column_rank(y)
        else                    y <- apply(y, 2, rank)
    }

    n <- nrow(x)
    
    sd_x <- sqrt(n * colSums(x*x) - colSums(x) ^ 2)
    sd_y <- sqrt(n * colSums(y*y) - colSums(y) ^ 2)
    
    cor_mat <- (n * (t(x) %*% y) - tcrossprod(colSums(x), colSums(y))) /
        tcrossprod(sd_x, sd_y)
    return(as.matrix(cor_mat))
}

#' Pearson or spearman correlation on sparse matrix
#' Finds only a subset of pairwise correlations based on the indices given in
#' x_idx and y_idx
sparse.cor.partial <- function(x, y, x_idx, y_idx, method=c("pearson", "spearman")) {
    stopifnot(nrow(x) == nrow(y))
    stopifnot(length(x_idx) == length(y_idx))
    
    stopifnot(ncol(x) >= max(x_idx) && 1 <= min(x_idx))
    stopifnot(ncol(y) >= max(y_idx) && 1 <= min(y_idx))

    stopifnot(is(x, "dgCMatrix") || is(x, "matrix"))
    stopifnot(is(y, "dgCMatrix") || is(y, "matrix"))
    method <- match.arg(method)

    if (method == "spearman") {
        if (is(x, "dgCMatrix")) x <- column_rank(x)
        else                    x <- apply(x, 2, rank)
        
        if (is(y, "dgCMatrix")) y <- column_rank(y)
        else                    y <- apply(y, 2, rank)
    }
    
    n <- nrow(x)
    
    sd_x <- sqrt(n * colSums(x*x) - colSums(x) ^ 2)
    sd_y <- sqrt(n * colSums(y*y) - colSums(y) ^ 2)

    # Call appropriate C++ method to do partial dot product
    if(is(x, "matrix") && is(y, "matrix")) {
        dotprod <- partialDot_dense_dense_cpp(x, y, x_idx - 1, y_idx - 1)
    } else if (is(x, "dgCMatrix") && is(x, "dgCMatrix")) {
        dotprod <- partialDot_sparse_sparse_cpp(x, y, x_idx - 1, y_idx - 1)
    } else {
        if (is(x("dgCMatrix"))) {
            dotprod <- partialDot_sparse_dense_cpp(x, y, x_idx - 1, y_idx - 1)
        } else {
            dotprod <- partialDot_sparse_dense_cpp(y, x, y_idx - 1, x_idx - 1)
        }
    }

    cor <- (n * dotprod - colSums(x)[x_idx]*colSums(y)[y_idx]) /
        (sd_x[x_idx]*sd_y[y_idx])

    return(cor)
}


#' Helper function that computes column-wise shifted rank transform on either
#' sparse or dense matrices
col_rank_transform <- function(x) {
    if (is(x, "dgCMatrix"))  x <- column_rank_cpp(x)
    else if(is(x, "matrix")) x <- apply(x, 2, rank)
    else stop("Matrix must be dgCMatrix or matrix")
}

#' Calculate permutation-based background mean + standard deviation for correlation
#' results. Numerically stable algorithm (Welford)
#' @param x Matrix of observations x features
#' @param y Matrix of observations x features (This matrix gets row-permuted, so should have smaller dimensions for best performance)
#' @param n Number of permutation iterations to run
#' @param replace Whether to use replacement when sampling rows
#' @param seed Seed to use for sampling
#' @param calculate_means Whether to assume mean correlations are 0, or calculate using re-sampling
#' @param method What kind of correlation to use, either "pearson" or "spearman"
#' @param pseudobulk_mat If non-null a matrix of dimension n-pseudobulks x n-observations, which should
#'      left-multiply both X and Y prior to running correlation. Note that if a pseudobulk_mat is given,
#'      both X and Y will have their rows separately permuted
#' @param save_iterations Boolean of whether to save correlation values from all iterations in an array 
#'      (dimensions ncol(x), ncol(y), n) 
sparse.cor.null <- function(x, y, n=50, replace=FALSE, seed=125125, calculate_means=FALSE, method=c("pearson", "spearman"),
                            pseudobulk_mat=NULL, save_iterations=FALSE) {
    method <- match.arg(method)
    if(!is.null(seed)) {
        set.seed(seed)
    }

    if (save_iterations) iterations <- list()
    else iterations <- NULL

    means <- matrix(0, nrow=ncol(x), ncol=ncol(y), dimnames = list(colnames(x), colnames(y)))
    m2 <- matrix(0, nrow=ncol(x), ncol=ncol(y))

    # Convert to ranks just once if needed for spearman and we aren't pseudobulking
    if (method == "spearman" && is.null(pseudobulk_mat)) {
        x <- col_rank_transform(x)
        y <- col_rank_transform(y)   
    }
    
    for (i in seq_len(n)) {
        shuffle <- sample.int(nrow(y), replace=replace)
        y_shuffle <- y[shuffle,]
        if(!is.null(pseudobulk_mat)) {
            y_shuffle <- as.matrix(pseudobulk_mat %*% y_shuffle)
            if(method == "spearman") y_shuffle <- col_rank_transform(y_shuffle)

            # We also have to shuffle x if we have a pseudobulk mat
            shuffle2 <- sample.int(nrow(x), replace=replace)
            x_shuffle <- x[shuffle2,]
            if(!is.null(pseudobulk_mat)) {
                x_shuffle <- as.matrix(pseudobulk_mat %*% x_shuffle)
                if(method == "spearman") x_shuffle <- col_rank_transform(x_shuffle)
            }
        } else {
            # No secondary shuffling needed if there isn't a pseudobulk mat
            x_shuffle <- x
        }

        res <- sparse.cor(x_shuffle, y_shuffle)
        
        # Adapted from https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
        if (calculate_means) {
            delta <- res - means
            means <- means + delta / i
            delta2 <- res - means
            m2 <- m2 + delta * delta2
        } else {
            m2 <- m2 + res*res
        }
        if (save_iterations) iterations[[i]] <- res
    }
    return(list(
        mean=means,
        sd=sqrt(m2/(n-1)),
        iterations=iterations
    ))
}


#' Returns a p-value for correlation based on T-test
#' https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Testing_using_Student's_t-distribution
#' I believe this is the same test that R's cor.test uses
#' @param correlations matrix/vector/numeric pearson (or spearman) correlation coefficients
#' @param npoints (single number) Number of points used for all the correlations (must be the same)
#' @return p-value for two-sided test with null hypothesis of zero correlation
cor.pval <- function(correlations, npoints) {
    t <- abs(correlations)*sqrt((npoints-2)/(1-correlations^2))
    return(
        pt(t, npoints-2, lower.tail=FALSE)*2
    )
}

#' @param x Matrix of observations x features
#' @param y Matrix of observations x features 
#' @param use_permutations whether to calculate (semi-)empirical p values with sparse.cor.null
#' @param permutation_args named list of additional arguemnts to sparse.cor.null 
#' @param method Correlation method to use, can be pearson or spearman
#' @param pseudobulk_mat If non-null a matrix of dimension n-pseudobulks x n-observations, which should
#'      left-multiply both X and Y prior to running correlation
#' @return Data table with columns feature_x, feature_y, pearson, p.ttest, (optionally p.permutation)
correlation_test <- function(x, y, use_permutations=TRUE, permutation_args=NULL, method=c("pearson", "spearman"),
                             pseudobulk_mat=NULL) {
  method <- match.arg(method)

  if (is.null(pseudobulk_mat)) {
    cor <- sparse.cor(x, y, method=method)
  } else {
    cor <- sparse.cor(as.matrix(pseudobulk_mat %*% x), as.matrix(pseudobulk_mat %*% y), method=method)
  }
  
  p.ttest <- cor.pval(cor, nrow(x))
  
  data_tibbles <- list(
    cor = as_tibble(cor, rownames="feature_x") %>% rename_with(~paste0(method, ".", .x), !feature_x),
    p.ttest = as_tibble(p.ttest) %>% rename_with(~paste0("p.ttest.", .x))
  )

  if(use_permutations) {
    if (is.null(permutation_args)) {
        permutation_args <- list()
    }
    permutation_args[["x"]] <- x
    permutation_args[["y"]] <- y
    permutation_args[["method"]] <- method
    permutation_args[["pseudobulk_mat"]] <- pseudobulk_mat
    background <- do.call(sparse.cor.null, permutation_args)
    data_tibbles[["means"]] <- as_tibble(background$mean) %>% rename_with(~paste0("mean.", .x))
    data_tibbles[["sds"]] <- as_tibble(background$sd) %>% rename_with(~paste0("sd.", .x))
  }

  res <- bind_cols(data_tibbles) %>%
    pivot_longer(!feature_x, 
        names_to=c(".value", "feature_y"), 
        names_pattern=c("(pearson|spearman|p\\.ttest|mean|sd)\\.(.*)"))
    
  if(use_permutations) {
    res <- mutate(res, p.permutation=pnorm(abs((.data[[method]]-mean)/sd), lower.tail=FALSE)*2) %>%
        select(!c("mean", "sd"))
  }
  return(res)
}




#' Get a MAGIC-style imputation matrix
#' @param pca cells x dimensions linear dimensionality reduction
#' @param k Number of nearest neighbors considered
#' @param k_a Index of nearest neighbor normalized to distance 1 (generally 1/3 k)
#' @param distance_metric Distance metric to use for nearest neighbor search
#' @return Sparse matrix of cellsxcells that can left-multiply a cellxfeature matrix
#'   in order to perform one diffusion iteration
#' Notes: So the MAGIC paper and 4 "official" implementations disagree on algorithm
#' details. I'm going with the version in dpeerlab/magic R implementation, which
#' is also what ArchR uses. Compared to the manuscript 
#' (methods section "Constructing MAGIC’s Markov Affinity Matrix")
#' this implemetation has two changes:
#' 1. Steps 2+3 are swapped (now symmetry is done before exponential transform)
#' 2. Exponential transform is just e^-x instead of e^(-x^2) (only the manuscript
#'    lists it this way, so I'm comfortable making the change)
get_impute_matrix <- function(pca, k_a=4, k=15, distance_metric=c("cosine", "euclidean")) {
    distance_metric <- match.arg(distance_metric)
    knn <- RcppHNSW::hnsw_knn(pca, k = k, distance = distance_metric)

    W <- sparseMatrix(
        x = as.numeric(knn$dist / knn$dist[,k_a]),
        i = rep(seq_len(nrow(pca)), k),
        j = as.integer(knn$idx),
        dims = c(nrow(pca), nrow(pca))
    )

    W <- (W + t(W))
    W@x <- exp(-(W@x))
    W <- W / rowSums(W)
    return(W)
}

#' Use a MAGIC-style imputation matrix
#' @param data cells x features matrix
#' @param impute_matrix Imputation matrix returned from get_impute_matrix
#' @param t Number of diffusion steps to run
impute_data <- function(data, impute_matrix, t=3) {
    rownames <- rownames(data)
    for (i in seq_len(t)) {
        data <- impute_matrix %*% data
    }
    data <- as.matrix(data)
    rownames(data) <- rownames
    return(as.matrix(data))
}

#' Load a text-based pwm file containing counts data (used for loading Helios motifs)
load_pwm <- function(pfm_path, name) {
  #Try following the normalization setup from chromVARmotifs
  read_tsv(pfm_path, col_types = "ddddc") %>%
    select(A, C, G, T) %>% as.matrix() %>%
    {(. / rowSums(.)) + .008 } %>%
    {log(. / rowSums(.) / .25)} %>%
    t() %>%
    PWMatrix(ID=name, name="Helios", strand="*", profileMatrix = .)
}


#' Create pseudobulks in the ArchR/Cicero style
#' @param pca Cells x dimensions matrix for neighbor finding
#' @param k Number of cells per pseudobulk
#' @param starter_cells Number of starter cells to use for generating pseudobulks prior to overlap filtering.
#'    Alternatively, if the length is >1 starter_cells is interpreted as a vector of row indices for the starter cells in the PCA
#' @param max_overlap Maximum percent overlaping cells in two pseudobulks
#' @param hnsw_build_args Additional arguments to pass to RcppHNSW::hnsw_build
#' @param hnsw_search_args Additional arguments to pass to RcppHNSW::hnsw_search
#' @return  pseudobulks x cells sparse matrix, where a 1 means that a cell should be counted
#'   as part of a pseudobulk
#' @details Algorithm approach: 1. Sample 'iterations' cells to use as the pseudobulk
#' starters. 2. Find nearest neighbors for each sampled cell (including self) and use as the pseudobulk
#' candidates. 3. Iteratively filter out pseudobulks by arbitrarily discarding whenever
#' we detect a pseudobulk with greater than the maximum overlap.
#' Note that a minor difference between this and the ArchR approach is that ArchR
#' appears to not include the initial sampled cell in the pseudobulk, whereas I will
#' include both the cell and its neighbors
#' https://github.com/GreenleafLab/ArchR/blob/968e4421ce7187a8ac7ea1cf6077412126876d5f/R/IntegrativeAnalysis.R#L736-L750
knn_pseudobulk <- function(pca, k=100, starter_cells=500, max_overlap=0.8, hnsw_build_args=list(), hnsw_search_args=list()) {
    if (length(starter_cells) == 1) 
        starter_cells <- sample.int(nrow(pca), starter_cells, replace= !nrow(pca) > starter_cells)
    
    hnsw_build_args[["X"]] <- pca
    knn_index <- do.call(RcppHNSW::hnsw_build, hnsw_build_args)
    hnsw_search_args[["X"]] <- pca[starter_cells,]
    hnsw_search_args[["ann"]] <- knn_index
    hnsw_search_args[["k"]] <- k
    knn <- do.call(RcppHNSW::hnsw_search, hnsw_search_args)

    groups_mat <- sparseMatrix(
        x = 1,
        i = as.integer(knn$idx),
        j = rep(seq_along(starter_cells), k),
        dims = c(nrow(pca), length(starter_cells))
    )

    overlaps <- as.matrix(t(groups_mat) %*% groups_mat) > floor(max_overlap * k)
    diag(overlaps) <- 0
    # Iteratively discard the pseudobulk with the most overlaps to other pseudobulks
    keep <- filter_pseudobulk_helper_cpp(overlaps)

    return(t(groups_mat[,keep]))
}