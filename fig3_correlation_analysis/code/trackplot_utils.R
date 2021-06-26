
# Design: The aim is to make genome-browser style plots,
# where stacked plots are aligned to share a horizontal axis representing
# genome location. 

# Each component plot function returns either a plain ggplot or a list of ggplot objects.
# These are then combined using draw_trackplot_grid, which handles alignment of 
# all the plots with labels

#' Combine ggplot track plots into an aligned grid.
#' Uses patchwork to perform the alignment
#' @param ... Plots in order from top to bottom, generally plain ggplots. 
#'    To better accomodate many bulk tracks,
#'    it is also possible to pass a list of plots as an argument. For plots in a list, the
#'    labels argument will be ignored and labels will be taken from the list names. 
#'    All list elements will be given the same relative height from the rel_heights list. 
#'    Additionally, only a single shared legend will be displayed in the legends section.
#' @param labels Text labels to display for each track
#' @param title Text for overarching title of the plot
#' @param heights Relative heights for each component plot. It is suggested to use 1 as standard height of a 
#'    pseudobulk track.
#' @param label_width Fraction of width that should be used for labels relative to the main track area
#' @param label_style Arguments to pass to geom_text to adjust label text style
#' @return A plot object with aligned genome plots. Each aligned row has
#'    the text label, y-axis, and plot body. The relative height of each row is given
#'    by heights. A shared title and x-axis are put at the top.
draw_trackplot_grid <- function(..., labels, title=NULL, 
        heights=rep(1, length(plots)), 
        label_width=0.2,
        label_style=list(fontface="bold", size=4)
        ) {
    if (!requireNamespace("patchwork", quietly=TRUE))
        stop("patchwork must be installed to se the grid functionality")
    plots <- list(...)
    if (length(plots) != length(labels) || 
        length(plots) != length(heights)) {
            stop("plots, labels, and heights must all have same length")
    }

    # Handle lists of plots by flattening appropriately
    # This is a bit awkward, but will hopefully save our function callers from inconvenience
    list_length <- vapply(plots, function(p) {if(is(p, "list")) length(p) else 0}, 0)
    if (any(list_length != 0)) {
        total_elements <- sum(list_length) + sum(list_length == 0)
        next_index <- 1

        new_plots <- list()
        new_labels <- rep("", total_elements)
        new_heights <- rep(0, total_elements)

        for(i in seq_along(plots)) {
            if(list_length[i] == 0) {
                new_plots[[next_index]] <- plots[[i]]
                new_labels[next_index] <- labels[i]
                new_heights[next_index] <- heights[i]
                next_index <- next_index + 1
            } else {
                for (j in seq_along(plots[[i]])) {
                    new_plots[[next_index]] <- plots[[i]][[j]]
                    new_labels[next_index] <- names(plots[[i]])[j]
                    new_heights[next_index] <- heights[i]
                    next_index <- next_index + 1
                }
            }
        }
        plots <- new_plots
        labels <- new_labels
        heights <- new_heights
    }


    x_axis_height <- 0.5

    # Make column of plot labels
    labels_plots <- lapply(labels, function(label_text) {
        ggplot(NULL, aes(0,0,label=label_text)) + do.call(geom_text, label_style) +
            theme_void()
    })
    
    data_plots <- lapply(plots, function(p) {
        p + theme(
            axis.title.x = element_blank(),
            axis.ticks.length.x = unit(0, "pt"),
            axis.text.x = element_blank()
        )
    })
    data_plots[[1]] <- plots[[1]]

    patch <- Reduce(`+`, c(labels_plots, data_plots)) + 
        patchwork::plot_layout(ncol=2, byrow=FALSE, widths=c(label_width, 1), heights=heights, guides="collect")
    if(!is.null(title))
        patch <- patch + patchwork::plot_annotation(title=title, theme=theme(plot.title=element_text(hjust=0.5)))
    return(patch)
}

#' Construct pseudobulk trackplot 
#' @param region GRanges of length 1 with region to plot
#' @param insertions Insertions object
#' @param groups Group assignments for each cell in "insertions"
#' @param group_norms Named numeric vector of normalization factors to divide insertion counts by for each group
#' @param palette Named character vector of color values for each group
#' @param bin_width Bin width in basepairs for counting insertions
#' @param crop_quantile Quantile of values for clipping y-axis limits. Default of 0.999 will crop out
#'    just the most extreme outliers across the region
#' @param legend_label Custom label to put on the legend
#' @param counts_ylabel Boolean whether to display raw counts on y-axis label (TRUE), 
#'    or counts scaled by group_norms (FALSE, default)
#' @return Named list of bulk genome track plots
trackplot_bulk <- function(region, insertions, groups, group_norms, palette, bin_width=200, crop_quantile=0.999,
                           legend_label = "group", counts_ylabel=FALSE) {
    # TODO: remove top+bottom margins & equalize height between the tracks
    if (!is.factor(groups))
        groups <- as.factor(groups)
    group_levels <- levels(groups)

    if (length(setdiff(group_levels, names(group_norms))) > 0) {
        stop("group_norms must have named entries for each level in groups")
    }

    tiles <- GenomicRanges::slidingWindows(region, bin_width, step=bin_width)[[1]]
    # Chop off any trailing bin that is not a full bin_width
    if (width(tiles[length(tiles)]) < bin_width)
        tiles <- tiles[-length(tiles)]
    cell_counts <- peakMatrix(insertions, tiles)
    group_counts <- lapply(group_levels, function(group) {
        rowSums(cell_counts[,groups == group])
    })
    group_counts <- do.call(rbind, group_counts)
    rownames(group_counts) <- group_levels
    rm(cell_counts)
    norm_group_counts <- group_counts / as.vector(group_norms[group_levels])

    ymax <- quantile(as.vector(norm_group_counts), crop_quantile)

    track_data <- tibble(
        x = rep(start(tiles) + bin_width/2, each=length(group_levels)),
        y = pmin(as.vector(norm_group_counts), ymax),
        group = factor(rep(group_levels, length(tiles)), levels=group_levels),
        group_norm = rep(group_norms[group_levels], length(tiles))
    )

    trackplots <- list()
    
    for (i in seq_along(group_levels)) {
        group <- group_levels[i]
        if (counts_ylabel) {
            ylim <- c(0, ymax * group_norms[group])
            names(ylim) <- c("", round(ylim[2]))
            mapping <- aes(x, y * group_norm, fill=group)
        } else {
            ylim <- c(0, ymax)
            names(ylim) <- c("", sprintf("%0.2e", ylim[2]))
            mapping <- aes(x, y, fill=group)
        }
        
        plot <- ggplot(
                filter(track_data, group == !!group), 
                mapping
            ) +
            geom_area() +
            scale_fill_manual(values=palette, drop=FALSE) +
            scale_x_continuous(limits=c(start(region), end(region)), expand=c(0,0), position="top") +
            scale_y_continuous(breaks=ylim, limits=ylim, expand=c(0,0)) +
            labs(x=NULL, y=NULL, fill=legend_label) +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                plot.margin = unit(c(0,0,0,0), "pt")
            )

        trackplots[[group]] <- plot
    }
    return(trackplots)
}


#' Plot single cell binarized accessibility track
#' @param region GRanges of length 1 with region to plot
#' @param insertions Insertions object
#' @param cell_order Character vector of cell_ids in order to display from top to bottom
#' @param bin_width width of bins to use
#' @return Single cell insertions plot
trackplot_singlecell <- function(region, insertions, cell_order, bin_width=250) {
    tiles <- GenomicRanges::slidingWindows(region, bin_width, step=bin_width)[[1]]

    counts <- peakMatrix(insertions, tiles)[,cell_order]
    #binarize counts
    counts@x <- pmin(counts@x, 1)
    counts <- as(counts, "dgTMatrix")

    data <- tibble(
        x = start(tiles)[counts@i] + bin_width/2,
        y = length(cell_order) + 1 - counts@j
    )
    plot <- ggplot(data, aes(x=x, y=y)) +
        geom_tile(width=bin_width, height=1) +
        labs(x=NULL, y=NULL) +
        scale_x_continuous(limits=c(start(region), end(region)), expand=c(0,0)) +
        scale_y_continuous(labels=NULL, breaks=NULL, expand=c(0,0)) +
        theme_bw() + theme(panel.grid=element_blank())

    return(plot)
}

#' @param values Character vector of discrete labels to present
#' @param palette If discrete is TRUE, named character vector of colors for each label
#'                If discrete is FALSE, character vector of colors to interpolate with
#' @param discrete Boolean for whether to use continuous or discrete color scale
sideplot_singlecell_colorbar <- function(values, palette, discrete=TRUE, label_text=NULL) {
    data <- tibble(
        y = rev(seq_along(values)),
        x = 1,
        color = values
    )
    if (discrete) fill_scale <- scale_fill_manual(values=palette)
    else          fill_scale <- scale_fill_gradientn(colors=palette)
    plot <- ggplot(data, aes(x=x, y=y, fill=color)) +
        geom_raster() +
        fill_scale +
        labs(x=NULL, y=NULL) +
        scale_x_continuous(breaks=NULL) +
        scale_y_continuous(labels=NULL, breaks=NULL, expand=c(0,0))

    legend <- cowplot::get_legend(plot)

    plot <- plot + guides(fill=FALSE)

    return(list(
        label = cowplot::ggdraw() + cowplot::draw_label(label_text, angle=90),
        plot = plot,
        legend = legend
    ))
}



#' Plot transcript models for a trackplot
#' @param region GRanges of length 1 with region to plot
#' @param transcripts GRanges of transcript features to plot. This will generally include
#'   several items per transcript in order to list exons, introns, coding sequences, etc.
#' @param transcript_ids Character list of transcript IDs for each item in transcripts
#' @param labels Character vector with labels for each item in transcripts. NA for items that should not be labeled
#' @param widths Numeric vector with relative line width of each item in transcripts (generally exons are drawn wider)
#' @param label_size size for transcript labels
#' @param width_range Length 2 vector with minimum and maximum width to use
#' @return Plot of gene locations
trackplot_gene <- function(region, transcripts, transcript_ids, labels, widths, label_size=2, width_range=c(1,3)) {
    in_region <- GenomicRanges::countOverlaps(transcripts, region) > 0
    transcripts <- transcripts[in_region] 

    
    data <- tibble(
        strand = as.vector(strand(transcripts)), 
        start = ifelse(strand == "-", end(transcripts), start(transcripts)),
        end = ifelse(strand == "-", start(transcripts), end(transcripts)),
        transcript_id = transcript_ids[in_region],
        label = labels[in_region],
        size = widths[in_region]
    )
    # Adjust the endpoints of any partially overlapping elements to fit within
    # the plot boundaries
    data <- mutate(
        data,
        start = pmax(start(region), pmin(end(region), start)),
        end =   pmax(start(region), pmin(end(region), end)),
    )
    plot <- ggplot(
        data,
        aes(x=start, xend=end, 
            y=transcript_id, yend=transcript_id, 
            size=size)
    ) +
        geom_segment(aes(color=strand)) +
        ggrepel::geom_text_repel(
            data=filter(data, !is.na(label)), 
            aes(label=label), 
            size=label_size,
            position=ggrepel::position_nudge_repel(y=0.25)
        ) +
        scale_size(range=c(1,3)) +
        scale_color_manual(values=c("+"="black", "-"="darkgrey")) +
        scale_x_continuous(limits=c(start(region), end(region)), expand=c(0,0)) +
        scale_y_discrete(labels=NULL, breaks=NULL) +
        labs(x=NULL, y=NULL) +
        guides(size=FALSE) +
        theme_bw() + theme(panel.grid=element_blank())

    return(plot)
}


#' Plot peak locations for a trackplot
#' @param region GRanges of length 1 with region to plot
#' @param peaks GRanges of peaks to plot.
#' @param peak_groups Optional character vector of group names for each peak
#' @param palette Optional named character vector with a color for each item in peak_groups
#' @param width Number with width for peak line segments
#' @return List with 3 entries: 
#'   - "label", Label for the track (will be the text "Peaks")
#'   - "plot", Peak location track plot object
#'   - "legend", color legend
trackplot_peak <- function(region, peaks, peak_groups, palette, width=2) {
    in_region <- GenomicRanges::countOverlaps(peaks, region) > 0
    data <- tibble(
        start = start(peaks[in_region]),
        end = end(peaks[in_region]),
        group = peak_groups[in_region]
    )
    plot <- ggplot(data, aes(x=start, xend=end, 
                             y=group, yend=group, color=group)) +
        geom_segment(size=width) +
        scale_color_manual(values=palette) +
        labs(x=NULL, y=NULL) +
        scale_x_continuous(limits=c(start(region), end(region)), expand=c(0,0)) +
        theme_bw() + theme(panel.grid=element_blank())

    return(plot)
}


# Loop track
# - pass in lists of start, end, and score
# - One function for loop track, other for color legend

#' Plot loops for a trackplot
#' @param region GRanges of length 1 with region to plot
#' @param loops GRanges of loops to plot. Each item will turn into a loop from start coordinate to end coordinate
#' @param scores Optional numeric vector of scores to color each loop
#' @param palette Optional character vector of color levels to interpolate
#' @param width Number with width for loops
#' @return List with 3 entries: 
#'   - "label", Label for the track (will be the text "Loops")
#'   - "plot", Loop track plot object
#'   - "legend", color legend
trackplot_loop <- function(region, loops, scores, palette, width=2) {
    in_region <- GenomicRanges::countOverlaps(loops, region, type="within") > 0
    
    # It turns out to be difficult to make geom_curve display well for a loop track,
    # since it cuts off the bottoms of loops and will only place the endpoints vertically in the middle.
    # As a work-around, we'll follow the ArchR approach of calculating points in an arc from a circle,
    # then manually drawing lines in the shape of an arc for each loop
    degrees <- 90 * pi/180
    
    curve <- tibble(
        angle = seq(-pi/2 - degrees/2, -pi/2 + degrees/2, length.out=100),
        x_raw = cos(angle),
        y_raw = sin(angle),
        # Adjust x and y coordinates so we have an arc that goes from 0 to 1 on the x axis
        x = (x_raw - min(x_raw)) / (max(x_raw) - min(x_raw)),
        y = (y_raw - max(y_raw)) / (max(x_raw) - min(x_raw))
    )

    loops <- tibble(
        start = start(loops[in_region]),
        end = end(loops[in_region]),
        score = scores[in_region],
        id = row_number(score) # Order loops for display
    )

    # Duplicate our curve coordinates for each loop, then recalculate the points
    data <- expand_grid(curve, loops)
    data <- transmute(data, 
        x = x * (end-start) + start,
        y = y * (end-start),
        id = id,
        score=score
    )
    data <- arrange(data, id)

    # if(nrow(data) == 0) {
    #     return(list(
    #         label = cowplot::ggdraw() + cowplot::draw_label(label_text),
    #         plot = NULL,
    #         legend = NULL
    #     ))
    # }
    plot <- ggplot(data, aes(x=x, y=y, group=id, color=score)) +
        geom_line(size=width) +
        scale_color_gradientn(colors=palette) +
        labs(x=NULL, y=NULL) +
        scale_x_continuous(limits=c(start(region), end(region)), expand=c(0,0)) +
        scale_y_continuous(limits=c(min(data$y * 1.05), 0), labels=NULL, breaks=NULL, expand=c(0,0)) +
        theme_bw() + theme(panel.grid=element_blank())

    return(plot)
}
