read_10x_gene_h5 <- function(path) {
    h5 <- hdf5r::H5File$new(path, mode="r")
    rna <- Matrix::sparseMatrix(
        i=h5[["matrix/indices"]][]+1, x=h5[["matrix/data"]][], p=h5[["matrix/indptr"]][], dims = h5[["matrix/shape"]][]
    )
    rownames(rna) <- h5[["matrix/features/id"]][]
    colnames(rna) <- h5[["matrix/barcodes"]][]
    rna <- rna["Gene Expression" == h5[["matrix/features/feature_type"]][], ]
    h5$close()
    rna
}

read_10x_gene_names_h5 <- function(path) {
    h5 <- hdf5r::H5File$new(path, mode="r")
    res <- tibble::tibble(
        gene_id = h5[["matrix/features/id"]][],
        gene_name = h5[["matrix/features/name"]][]
    )
    h5$close()
    res
}



#' @param mat raw counts matrix, features x cells
#' @return log10 of CLR-normalized counts
clr_normalize <- function(mat, base=2) {
    t(t(log(mat + 1, base)) - colMeans(log(mat + 1, base)))
}

#' Compute cluster means of single cell matrix
#' @param mat input matrix (1 column per cell)
#' @param clusters cluster identites (length = ncol(mat))
cluster_means <- function(mat, clusters) {
    stopifnot(length(clusters) == ncol(mat))
    cluster_ids <- levels(as.factor(clusters))
    output <- lapply(cluster_ids, function(id) {
        Matrix::rowMeans(mat[,clusters == id])
    })
    output <- do.call(cbind, output)
    colnames(output) <- cluster_ids
    rownames(output) <- rownames(mat)
    output
}

#' Compute cluster standard deviations of single cell matrix
#' @param mat input matrix (1 column per cell)
#' @param clusters cluster identites (length = ncol(mat))
cluster_sds <- function(mat, clusters) {
    stopifnot(length(clusters) == ncol(mat))
    cluster_ids <- levels(as.factor(clusters))
    output <- lapply(cluster_ids, function(id) {
        matrixStats::rowSds(mat[,clusters == id])
    })
    output <- do.call(cbind, output)
    colnames(output) <- cluster_ids
    rownames(output) <- rownames(mat)
    output
}



#' Demultiplex samples based on data from cell hashtag oligos. 
#' 
#' Observation: on my existing dataset, the kmedoids clustering is kinda meh.
#' For that dataset I should definitely use n and not n+1 clusters
#' 
#' Based on Seurat's HTODemux approach, but with two main differences
#' 1. Rather than fitting the cutoff based on a negative binomial of the raw counts
#'    data, I use either a normal distribution or empirical distribution of the
#'    normalized counts. I justify this because contaminating read counts scale with
#'    on-target read counts, so our cutoffs should use the normalized data, not raw
#'    read counts
#' 2. Instead of using the minimum kmedoids cluster for the background distribution, I use all 
#'    but the top 2 clusters, assuming the top 2 clusters will be the positive and multiplet clusters. 
#'    This gave better results on my dataset
#' 
#' @param hto_matrix A filtered+normalized counts matrix of hashtag oligos x cells
#' @param cutoff.quantile What quantile of the null distribution to set as a positive cutoff
#' @param cutoff.model What model to use of the null distribution? Normal assumes normal distribution;
#'  empirical takes the quantile directly on the observed data
#' @param kmedoids.clusters The number of clusters to use in k medoids clustering. In most cases this should be the number
#'   of HTOs, however in cases with heavily imbalanced proportions decreasing this number may help compensate
#' @param background.clusters The number of clusters to use as a background distribution. In most cases this should be
#'   the kmedoids.clusters - 1, however in cases with heavily imbalanced proportions decreasing this number may help compensate
#' @param sample.labels Optional vector of sample labels for each cell. If not NULL, cells will be split by sample
#'   and cutoffs for each HTO will be calculated separately per sample
#' 
#' @return A list with two elements: 
#'  cutoffs - a tibble with positive cutoff values, one for each hashtag (and optionally sample label).
#'  hto_count - a vector of how many HTOs match as positive for each cell.
#'  hto_assignment - a vector of HTO assignments, NA for cells where hto_count != 1
#'  kmedoids_cluster - a vector of kmedoid cluster ID assignments for debugging purposes
#' 
#' Singlets will have hto_count == 1, and their sample will be given by top_hto.
#' 
#' @importFrom cluster clara
#' @importFrom matrixStats rowMins
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' 
#' @export
#'
HTO_demultiplex <- function(hto_matrix, 
    cutoff.quantile=0.99, cutoff.model=c("normal", "empirical"), 
    kmedoids.clusters=nrow(hto_matrix), background.clusters=kmedoids.clusters-1,
    sample.labels=NULL) 
{
    cutoff.model <- match.arg(cutoff.model)
    if(is.null(rownames(hto_matrix))) {
        hto_names <- seq_len(nrow(hto_matrix))
    } else {
        hto_names <- rownames(hto_matrix)
    }

    if (is.null(sample.labels)) {
        # 1. Cluster based on hashtag expression
        kmedoids <- cluster::clara(
            t(hto_matrix),
            kmedoids.clusters,
            samples=100
        )$clustering
        names(kmedoids) <- colnames(hto_matrix)

        # 2. Calculate mean of each hto in each cluster
        means <- cluster_means(hto_matrix, kmedoids)
        stopifnot(sum(means == 0) == 0) #Ensure every cluster has *some* expression of each HTO

        # 3. Calculate positive cutoff
        # Take all but the top two clusters as the negative distribution
        cutoff_list <- lapply(seq_len(nrow(hto_matrix)), function(i) {
            min_clusters <- names(sort(means[i,]))[seq_len(background.clusters)]
            values <- hto_matrix[i,kmedoids %in% min_clusters]
            if(cutoff.model == "normal") {
                cutoff <- mean(values) + qnorm(cutoff.quantile) * sd(values)
            } else if(cutoff.model == "empirical") {
                cutoff <- quantile(values, cutoff.quantile, type=2)
            }
            cutoff
        })
        cutoffs <- tibble(
            hto = hto_names,
            cutoff = unlist(cutoff_list)
        )
        
        # 4. Classify cells based on cutoffs
        hto_count <- colSums(hto_matrix > cutoffs$cutoff)
        hto_max_id <- apply(hto_matrix > cutoffs$cutoff, 2, which.max)
        hto_max_id[hto_count != 1] <- NA

        res <- list(
            cutoffs = cutoffs,
            hto_count = hto_count,
            hto_assignment = hto_names[hto_max_id],
            kmedoids_cluster = kmedoids
        )
        return(res)
    } else { # Split cells by sample, run separately, then combine results
        results <- lapply(unique(sample.labels), function(sample) {
            HTO_demultiplex(
                hto_matrix[,sample.labels==sample],
                cutoff.quantile=cutoff.quantile, cutoff.model=cutoff.model,
                kmedoids.clusters=kmedoids.clusters, background.clusters=background.clusters,
                sample.labels=NULL
            )
        })
        names(results) <- unique(sample.labels)
        cutoffs <- lapply(results, function(res) {res$cutoffs})
        hto_count <- rep_len(NA_integer_, ncol(hto_matrix))
        hto_assignment <- rep_len(NA_character_, ncol(hto_matrix))
        kmedoids_cluster <- rep_len(NA_character_, ncol(hto_matrix))
        for (sample in unique(sample.labels)) {
            hto_count[sample.labels == sample] <- results[[sample]]$hto_count
            hto_assignment[sample.labels == sample] <- results[[sample]]$hto_assignment
            kmedoids_cluster[sample.labels == sample] <- paste0(sample, ".", results[[sample]]$hto_count)
        }
    
        res <- list(
            cutoffs = bind_rows(cutoffs, .id="sample"),
            hto_count = hto_count,
            hto_assignment = hto_assignment,
            kmedoids_cluster = kmedoids_cluster
        )
        return(res)
    }
}