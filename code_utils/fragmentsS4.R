
setClass("Insertions",
    slots = c(
        cell_levels = "character",
        cell_ids = "list",
        iranges = "list"
    ),
    prototype = c(
        cell_levels = NA_character_,
        cell_ids = list(),
        iranges = list()
    )
)

# To implement: length; pretty printing, accessors

Insertions <- function(fragments_granges, cell_ids, cell_levels=NULL) {
    stopifnot(is(fragments_granges, "GRanges"))
    stopifnot(length(fragments_granges) == length(cell_ids))
    cell_ids <- as.character(cell_ids)
    if (is.null(cell_levels)) {
        cell_levels <- unique(cell_ids)
    } else {
        cell_levels <- as.character(cell_levels)
    }
    
    iranges <- lapply(seqlevels(fragments_granges), function (x) {
        ranges(fragments_granges[seqnames(fragments_granges) == x])
    })
    cell_ids <- lapply(seqlevels(fragments_granges), function (x) {
        ids <- cell_ids[as.vector(seqnames(fragments_granges) == x)]
        as(factor(ids, cell_levels), "Rle")
    })
    names(iranges) <- seqlevels(fragments_granges)
    names(cell_ids) <- seqlevels(fragments_granges)
    new("Insertions", cell_levels=cell_levels, cell_ids=cell_ids, iranges=iranges)
}

setValidity("Insertions", function(object) {
    if (!identical(names(object@cell_ids), names(object@iranges))) {
        return("Names of @cell_ids and @iranges must match")
    }
    problems <- NULL
    for (seq in names(object@iranges)) {
        iranges <- object@iranges[[seq]]
        cell_ids <- object@cell_ids[[seq]]

        if(!is(ranges(iranges), "IRanges")) {
            problems <- c(problems, sprintf("@iranges[[\"%s\"]] is not of class IRanges", seq))
        }
        if(!is(cell_ids, "Rle") || !is(runValue(cell_ids), "factor")) {
            problems <- c(problems, sprintf("@cell_ids[[\"%s\"]] is not of class factor-Rle", seq))
        } else if (!identical(levels(runValue(cell_ids)), object@cell_levels)) {
            problems <- c(problems, sprintf("@cell_ids[[\"%s\"]] does not have levels @cell_levels", seq))
        }
        if(length(iranges) != length(cell_ids)) {
            problems <- c(problems, 
                sprintf("@cell_ids[[\"%s\"]] and @cell_ids[[\"%s\"]] have unequal lengths", seq, seq))
        }
    }
    if (is.null(problems)) {
        return(TRUE)
    } else {
        return(problems)
    }
})

setMethod("[", "Insertions", function(x, i) {
    if(length(setdiff(i, x@cell_levels)) > 0)
        stop("i must be a subset of x@cell_levels")
    if(length(unique(i)) != length(i))
        stop("All elements of i must be unique")
    
    new_levels <- match(x@cell_levels, i)

    iranges <- list()
    cell_ids <- list()

    for (chr in names(x@iranges)) {
        keeper_fragments <- x@cell_ids[[chr]] %in% i

        raw_cell_ids <- as.integer(x@cell_ids[[chr]][keeper_fragments])
        cell_ids[[chr]] <- as(
            factor(i[new_levels[raw_cell_ids]], levels=i),
            "Rle"
        )
        iranges[[chr]] <- x@iranges[[chr]][keeper_fragments]
    }
    new("Insertions", cell_levels=i, cell_ids=cell_ids, iranges=iranges)
})

# Testing code for subsetting
# sub1 <- insertions[proj_all$cellNames[2:3]]
# sub2 <- insertions[proj_all$cellNames[3:2]]

# subset_peaks1 <- peakMatrix(sub1, getPeakSet(proj_all))
# subset_peaks2 <- peakMatrix(sub2, getPeakSet(proj_all))
# full_peaks <- peakMatrix(insertions, getPeakSet(proj_all))

# all.equal(full_peaks[,2:3], subset_peaks)
# all.equal(subset_peaks[,2:1], subset_peaks2)


peakMatrix <-  function(insertions, peakset, threads=1) {
    if(length(setdiff(as.character(seqnames(peakset)), names(insertions@iranges))) != 0) {
        stop("peakset contains seqnames not present in x")
    }
    triplets_per_seqlevel <- function(seq) {
        peak_indices <- which(as.vector(seqnames(peakset) == seq))
        peaks <- ranges(peakset[peak_indices])
        overlaps_start <- findOverlaps(start(insertions@iranges[[seq]]), peaks)
        overlaps_end <- findOverlaps(end(insertions@iranges[[seq]]), peaks)
        # First make a csparse mat to add up counts per cell, then save the triplets
        mat <- sparseMatrix(
            i = peak_indices[c(to(overlaps_start), to(overlaps_end))],
            j = as.integer(insertions@cell_ids[[seq]][c(from(overlaps_start), from(overlaps_end))]),
            x = 1L,
            dims = c(length(peakset), length(insertions@cell_levels))
        )
        triplets <- data.frame(
            i = mat@i+1,
            j = rep.int(seq_len(ncol(mat))-1, diff(mat@p))+1,
            x = mat@x
        )
        return(triplets)
    }
    triplets <- mclapply(as.character(unique(seqnames(peakset))), triplets_per_seqlevel, mc.cores=threads) %>%
        do.call(rbind, .)
    ret <- sparseMatrix(
        i = triplets$i, j = triplets$j, x = triplets$x,
        dims = c(length(peakset), length(insertions@cell_levels)))
    colnames(ret) <- insertions@cell_levels
    return(ret)
}

insertionProfile <- function(insertions, regions, cell_groups, cell_weights=NULL, threads=1) {
    if(length(setdiff(seqnames(regions), names(insertions@iranges))) != 0) {
        stop("regions contains seqnames not present in x")
    }
    width <- width(regions)[1]
    if(!all(width(regions) == width)) {
        stop("regions must all be identical width")
    }
    if(length(cell_groups) != length(insertions@cell_levels)) {
        stop("length(cell_groups) must equal length(insertions@cell_levels)")
    }
    if(!is.null(cell_weights) && length(cell_weights) != length(insertions@cell_levels)) {
        stop("lengnth(cell_weights) must equal length(insertions@cell_levels)")
    }
    cell_groups <- as.factor(cell_groups)
    triplets_per_seqlevel <- function(seq) {
        region_indices <- which(as.vector(seqnames(regions) == seq))
        region_iranges <- ranges(regions[region_indices])

        insertion_coords <- c(start(insertions@iranges[[seq]]), end(insertions@iranges[[seq]]))
        cell_ids <- as.vector(as.integer(c(insertions@cell_ids[[seq]], insertions@cell_ids[[seq]])))

        overlaps <- findOverlaps(insertion_coords, region_iranges)
        is_minus_strand <- as.vector(strand(regions[region_indices])[to(overlaps)] == "-")
        dist <- (insertion_coords[from(overlaps)] - start(region_iranges)[to(overlaps)]) *
            ifelse(is_minus_strand, -1, 1) +
            (width-1) * is_minus_strand
        
        if (is.null(cell_weights)) {
            x <- 1L
        } else {
            x <- cell_weights[cell_ids[from(overlaps)]]
        }
        # First make a csparse mat to add up counts, then save the triplets
        mat <- sparseMatrix(
            i = dist + 1,
            j = as.integer(cell_groups[cell_ids[from(overlaps)]]),
            x = x,
            dims = c(width, length(levels(cell_groups)))
        )
        
        triplets <- data.frame(
            i = mat@i+1,
            j = rep.int(seq_len(ncol(mat))-1, diff(mat@p))+1,
            x = mat@x
        )
        return(triplets)
    }
    triplets <- mclapply(as.character(unique(seqnames(regions))), triplets_per_seqlevel, mc.cores=threads) %>%
        do.call(rbind, .)
    ret <- sparseMatrix(
        i = triplets$i, j = triplets$j, x = triplets$x,
        dims = c(width, length(levels(cell_groups))))
    colnames(ret) <- levels(cell_groups)
    return(ret)
}

#' Read fragments from an arrow file. Runs a little slow still, but
#' an initial working draft. (One issue is that it has to use method="slow"
#' on .getFragsFromArrow since that fails if you only have read-access to files)
#' @param read_method "fast" if you have write access too, "slow" otherwise
ArchR_to_insertion <- function(proj, cellNames=NULL, read_method="fast") {
  if (is.null(cellNames)) {
    cellNames <- proj$cellNames
  }
  chrs <- as.character(seqnames(getGenomeAnnotation(proj)$chromSizes))
  iranges <- list()
  cell_ids <- list()
  for (chr in chrs) {
    iranges_chr <- list()
    for (ArrowFile in proj@sampleColData$ArrowFiles) {
      iranges_chr[[ArrowFile]] <- ArchR:::.getFragsFromArrow(
        ArrowFile, chr=chr, out="IRanges", cellNames=cellNames, method=read_method
      )
    }
    names(iranges_chr) <-  NULL
    sorted_iranges <- sort(do.call(c, iranges_chr))
    iranges[[chr]] <- sorted_iranges
    cell_ids[[chr]] <- as(factor(mcols(sorted_iranges)$RG, cellNames), "Rle")
  }
  new("Insertions", cell_levels=cellNames, cell_ids=cell_ids, iranges=iranges)
}