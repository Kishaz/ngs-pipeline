#!/usr/bin/env Rscript

# =============================================================================
# Ancestry Pipeline — Companion Visualization Script
#
# Reads the output files from the ancestry pipeline and generates:
#   1. PCA projection plots (PC1-PC5 pairwise panels)
#   2. Combined ancestry barplot figure:
#        [metadata strips]        ← optional (Race, Ethnicity, etc.)
#        [super-population bars]  ← K=5 continental ancestry
#        [sub-population bars]    ← K=23 sub-population ancestry (if available)
#      All panels share the same x-axis, aligned and ordered by descending
#      African (AFR) ancestry proportion.
#
# All plots saved as SVG and PNG (300 dpi).
#
# Usage:
#   Rscript plot_ancestry.R <results_directory> [output_prefix] [metadata_file]
#
# Examples:
#   # Without metadata
#   Rscript plot_ancestry.R ./results
#   Rscript plot_ancestry.R ./results my_cohort
#
#   # With metadata — adds annotation strips above the barplots
#   Rscript plot_ancestry.R ./results my_cohort metadata.tsv
#
# Metadata file format (tab-delimited, header required):
#   sample_id    Race
#   SKMEL1094A   White
#   WM4304       Black
#   YUCRATE      Hispanic
#
#   - First column MUST be "sample_id" (matching pipeline output)
#   - Additional columns become annotation strips (Race, Ethnicity, Sex, etc.)
#
# Required input files (in results_directory):
#   ancestry_pca_coordinates.tsv      — PCA coordinates
#   ancestry_summary_superpops.tsv    — Continental ancestry proportions
#   ancestry_summary_subpops.tsv      — Sub-population proportions (optional)
#
# Author: Stephen Mwambu (@Kishaz)
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(grid)
})

# -- Parse arguments ----------------------------------------------------------
args <- commandArgs(TRUE)
if (length(args) < 1) {
    cat("Usage: Rscript plot_ancestry.R <results_dir> [prefix] [metadata] [-c col1,col2,...]\n")
    cat("\n")
    cat("  results_dir   Path to ancestry pipeline output\n")
    cat("  prefix        Prefix for output filenames (default: ancestry)\n")
    cat("  metadata      File with sample_id + metadata columns (CSV/TSV/XLSX)\n")
    cat("  -c COLUMNS    Comma-separated column names to include as annotation strips\n")
    cat("                (default: all columns except sample_id)\n")
    cat("\n")
    cat("Examples:\n")
    cat("  Rscript plot_ancestry.R ./results cohort metadata.csv -c race\n")
    cat("  Rscript plot_ancestry.R ./results cohort metadata.xlsx -c race,gender,treatment\n")
    cat("  Rscript plot_ancestry.R ./results cohort metadata.csv    # all columns\n")
    quit(status = 1)
}

# Separate positional args from -c flag
positional <- c()
selected_cols <- NULL
i <- 1
while (i <= length(args)) {
    if (args[i] == "-c" && i < length(args)) {
        selected_cols <- trimws(unlist(strsplit(args[i + 1], ",")))
        i <- i + 2
    } else {
        positional <- c(positional, args[i])
        i <- i + 1
    }
}

results_dir   <- positional[1]
prefix        <- if (length(positional) >= 2) positional[2] else "ancestry"
metadata_file <- if (length(positional) >= 3) positional[3] else NULL
outdir        <- results_dir

# -- File paths ----------------------------------------------------------------
pca_file      <- file.path(results_dir, "ancestry_pca_coordinates.tsv")
superpop_file <- file.path(results_dir, "ancestry_summary_superpops.tsv")
subpop_file   <- file.path(results_dir, "ancestry_summary_subpops.tsv")

if (!file.exists(superpop_file)) stop("Required file not found: ", superpop_file)
has_pca     <- file.exists(pca_file)
has_subpops <- file.exists(subpop_file)

# -- Load metadata if provided -------------------------------------------------
has_metadata <- FALSE
metadata <- NULL
meta_cols <- c()

if (!is.null(metadata_file)) {
    if (!file.exists(metadata_file)) {
        stop("Metadata file not found: ", metadata_file)
    }

    # Auto-detect file format: .xlsx, .csv, .tsv, or tab-delimited
    ext <- tolower(tools::file_ext(metadata_file))
    if (ext %in% c("xlsx", "xls")) {
        if (!requireNamespace("readxl", quietly = TRUE))
            stop("readxl package required for Excel files. Install with: install.packages('readxl')")
        metadata <- as.data.frame(readxl::read_excel(metadata_file), stringsAsFactors = FALSE)
        cat("  Read metadata as Excel (.xlsx)\n")
    } else {
        # Sniff delimiter: read first line and check for tabs vs commas
        first_line <- readLines(metadata_file, n = 1, warn = FALSE)
        if (grepl("\t", first_line)) {
            metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
            cat("  Read metadata as TSV (tab-delimited)\n")
        } else if (grepl(",", first_line)) {
            metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
            cat("  Read metadata as CSV (comma-delimited)\n")
        } else {
            # Fallback: try read.delim, then read.csv
            metadata <- tryCatch(
                read.delim(metadata_file, stringsAsFactors = FALSE),
                error = function(e) read.csv(metadata_file, stringsAsFactors = FALSE)
            )
            cat("  Read metadata (auto-detected delimiter)\n")
        }
    }

    # Trim whitespace from column names (common issue with Excel exports)
    colnames(metadata) <- trimws(colnames(metadata))

    if (!"sample_id" %in% colnames(metadata)) {
        stop("Metadata file must have a 'sample_id' column. Found columns: ",
             paste(colnames(metadata), collapse = ", "))
    }
    all_meta_cols <- setdiff(colnames(metadata), "sample_id")

    # Apply -c column filter if specified
    if (!is.null(selected_cols)) {
        # Case-insensitive matching
        col_map <- setNames(all_meta_cols, tolower(all_meta_cols))
        matched <- col_map[tolower(selected_cols)]
        missing <- selected_cols[is.na(matched)]
        if (length(missing) > 0) {
            warning("Columns not found in metadata (ignored): ", paste(missing, collapse = ", "),
                    "\n  Available: ", paste(all_meta_cols, collapse = ", "))
        }
        meta_cols <- unname(matched[!is.na(matched)])
        cat("  Column filter (-c): showing ", paste(meta_cols, collapse = ", "), "\n")
    } else {
        meta_cols <- all_meta_cols
    }

    if (length(meta_cols) == 0) {
        warning("No valid metadata columns to display - ignoring metadata")
    } else {
        has_metadata <- TRUE
    }
}

cat("Results directory:", results_dir, "\n")
cat("Output prefix:    ", prefix, "\n")
cat("PCA coordinates:  ", ifelse(has_pca, "yes", "no (skipping PCA plots)"), "\n")
cat("Sub-populations:  ", ifelse(has_subpops, "yes", "no"), "\n")
cat("Metadata:         ", ifelse(has_metadata,
    paste0("yes (", paste(meta_cols, collapse = ", "), ")"), "no"), "\n\n")

# =============================================================================
# COLOR DEFINITIONS
# =============================================================================

superpop_colors <- c(
    "AFR" = "#E41A1C",
    "AMR" = "#FF7F00",
    "EAS" = "#4DAF4A",
    "EUR" = "#377EB8",
    "SAS" = "#984EA3"
)

subpop_to_superpop <- c(
    "YRI" = "AFR", "LWK" = "AFR", "GWD" = "AFR", "MSL" = "AFR",
    "ESN" = "AFR", "ACB" = "AFR", "ASW" = "AFR",
    "CLM" = "AMR", "PEL" = "AMR", "PUR" = "AMR", "MXL" = "AMR",
    "CHB" = "EAS", "JPT" = "EAS", "CDX" = "EAS", "CHS" = "EAS", "KHV" = "EAS",
    "CEU" = "EUR", "GBR" = "EUR", "TSI" = "EUR", "IBS" = "EUR", "FIN" = "EUR",
    "BEB" = "SAS", "GIH" = "SAS", "STU" = "SAS", "ITU" = "SAS", "PJL" = "SAS"
)

make_gradient <- function(base_hex, n) {
    base_rgb <- col2rgb(base_hex)
    sapply(seq(0.35, 0.95, length.out = max(n, 1)), function(f) {
        r <- as.integer(255 - f * (255 - base_rgb[1, ]))
        g <- as.integer(255 - f * (255 - base_rgb[2, ]))
        b <- as.integer(255 - f * (255 - base_rgb[3, ]))
        rgb(r, g, b, maxColorValue = 255)
    })
}

subpop_colors <- c()
for (sp in names(superpop_colors)) {
    subs <- sort(names(subpop_to_superpop[subpop_to_superpop == sp]))
    cols <- make_gradient(superpop_colors[sp], length(subs))
    names(cols) <- subs
    subpop_colors <- c(subpop_colors, cols)
}

# =============================================================================
# COMMON THEME — base size 14pt, proportional scaling
# =============================================================================

BASE_SIZE <- 14

theme_ancestry <- theme_bw(base_size = BASE_SIZE) +
    theme(
        plot.title       = element_text(size = BASE_SIZE + 4, face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(size = BASE_SIZE, hjust = 0.5, color = "grey40"),
        axis.title       = element_text(size = BASE_SIZE, face = "bold"),
        axis.text        = element_text(size = BASE_SIZE - 2),
        legend.text      = element_text(size = BASE_SIZE - 2),
        legend.title     = element_text(size = BASE_SIZE, face = "bold"),
        legend.key.size  = unit(0.5, "cm"),
        strip.text       = element_text(size = BASE_SIZE, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin      = margin(12, 12, 12, 12)
    )

# -- Save helpers --------------------------------------------------------------
save_plot <- function(plot_obj, filepath, width, height) {
    ggsave(paste0(filepath, ".png"), plot = plot_obj,
           width = width, height = height, units = "in", dpi = 300)
    ggsave(paste0(filepath, ".svg"), plot = plot_obj,
           width = width, height = height, units = "in")
    cat("  Saved:", paste0(filepath, ".png"), "\n")
    cat("  Saved:", paste0(filepath, ".svg"), "\n")
}

save_grob <- function(grob_obj, filepath, width, height) {
    png(paste0(filepath, ".png"), width = width, height = height,
        units = "in", res = 300)
    grid.draw(grob_obj)
    dev.off()

    svg(paste0(filepath, ".svg"), width = width, height = height)
    grid.draw(grob_obj)
    dev.off()

    cat("  Saved:", paste0(filepath, ".png"), "\n")
    cat("  Saved:", paste0(filepath, ".svg"), "\n")
}

# =============================================================================
# SAMPLE ORDERING — always by decreasing African (AFR) ancestry
# =============================================================================

superpops <- read.delim(superpop_file, stringsAsFactors = FALSE)

if ("AFR" %in% colnames(superpops) && "EUR" %in% colnames(superpops)) {
    # Order: high AFR (left) → admixed (center) → high EUR (right)
    # Sort by AFR - EUR descending: AFR-dominant = positive, EUR-dominant = negative
    superpops$afr_eur_diff <- superpops$AFR - superpops$EUR
    sample_order <- superpops$sample_id[order(-superpops$afr_eur_diff)]
    superpops$afr_eur_diff <- NULL
} else if ("AFR" %in% colnames(superpops)) {
    sample_order <- superpops$sample_id[order(-superpops$AFR)]
} else {
    sample_order <- superpops$sample_id
}

cat("Sample order (AFR-high → admixed → EUR-high):\n")
for (s in sample_order) {
    afr_val <- if ("AFR" %in% colnames(superpops)) superpops$AFR[superpops$sample_id == s] else NA
    eur_val <- if ("EUR" %in% colnames(superpops)) superpops$EUR[superpops$sample_id == s] else NA
    cat(sprintf("  %-20s  AFR = %5.1f%%  EUR = %5.1f%%\n", s,
                afr_val * 100, eur_val * 100))
}
cat("\n")

# -- Prepare matched metadata --------------------------------------------------
meta_matched <- NULL
if (has_metadata) {
    meta_matched <- metadata[metadata$sample_id %in% sample_order, ]
    n_matched <- nrow(meta_matched)
    n_total <- length(sample_order)
    if (n_matched < n_total) {
        missing <- setdiff(sample_order, meta_matched$sample_id)
        warning(sprintf("Metadata missing for %d sample(s): %s",
                        length(missing), paste(missing, collapse = ", ")))
        missing_df <- data.frame(sample_id = missing, stringsAsFactors = FALSE)
        for (mc in setdiff(colnames(meta_matched), "sample_id"))
            missing_df[[mc]] <- "N/A"
        meta_matched <- rbind(meta_matched,
                              missing_df[, colnames(meta_matched), drop = FALSE])
    }
}

# =============================================================================
# METADATA STRIP BUILDER
# =============================================================================

build_meta_strips <- function(meta_df, meta_columns, sample_levels) {
    strips <- list()
    for (col in meta_columns) {
        strip_df <- data.frame(
            sample_id = factor(meta_df$sample_id, levels = sample_levels),
            value     = as.character(meta_df[[col]]),
            stringsAsFactors = FALSE
        )

        p_strip <- ggplot(strip_df, aes(x = sample_id, y = 1, fill = value)) +
            geom_tile(colour = "white", linewidth = 0.3) +
            scale_y_continuous(expand = c(0, 0)) +
            scale_x_discrete(drop = FALSE) +
            coord_fixed(ratio = 1) +
            labs(x = NULL, y = col, fill = col) +
            theme_minimal(base_size = BASE_SIZE) +
            theme(
                axis.text.x      = element_blank(),
                axis.ticks.x     = element_blank(),
                axis.text.y      = element_blank(),
                axis.ticks.y     = element_blank(),
                axis.title.y     = element_text(size = BASE_SIZE, face = "bold",
                                                angle = 0, vjust = 0.5),
                panel.grid       = element_blank(),
                legend.position  = "right",
                legend.key.size  = unit(0.5, "cm"),
                legend.text      = element_text(size = BASE_SIZE - 2),
                legend.title     = element_text(size = BASE_SIZE, face = "bold"),
                plot.margin      = margin(4, 5.5, 4, 5.5)
            )

        strips[[col]] <- p_strip
    }
    return(strips)
}

# =============================================================================
# 1. PCA PLOTS (skipped if ancestry_pca_coordinates.tsv not found)
# =============================================================================

if (has_pca) {
    cat("Generating PCA plots...\n")
    pca <- read.delim(pca_file, stringsAsFactors = FALSE)

    # Handle is_study_sample as logical, string ("True"/"False"), or integer (0/1)
    pca$is_study <- tolower(as.character(pca$is_study_sample)) %in% c("true", "1", "yes")
    ref_data   <- pca[!pca$is_study, ]
    study_data <- pca[pca$is_study, ]

    cat("  Reference samples:", nrow(ref_data), "\n")
    cat("  Study samples:    ", nrow(study_data), "\n")

    if (nrow(study_data) == 0) {
        cat("  WARNING: No study samples found in PCA coordinates.\n")
        cat("  Column 'is_study_sample' values:", head(unique(pca$is_study_sample), 10), "\n")
        cat("  Skipping PCA plots.\n\n")
        has_pca <- FALSE
    }
}

if (has_pca) {
    pop_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS", "Study Sample")
    pca_colors <- c(superpop_colors, "Study Sample" = "#000000")
    pca_shapes <- c("AFR" = 16, "AMR" = 16, "EAS" = 16, "EUR" = 16,
                    "SAS" = 16, "Study Sample" = 18)
    pca_sizes  <- c("AFR" = 1.8, "AMR" = 1.8, "EAS" = 1.8, "EUR" = 1.8,
                    "SAS" = 1.8, "Study Sample" = 5.5)
    pca_alphas <- c("AFR" = 0.4, "AMR" = 0.4, "EAS" = 0.4, "EUR" = 0.4,
                    "SAS" = 0.4, "Study Sample" = 1.0)

    ref_data$pop   <- ref_data$population
    study_data$pop <- "Study Sample"
    pca_all <- rbind(ref_data, study_data)
    pca_all$pop <- factor(pca_all$pop, levels = pop_levels)

    make_pca_panel <- function(data, study_pts, xcol, ycol, xlab, ylab, add_labels) {
        p <- ggplot(data, aes(x = .data[[xcol]], y = .data[[ycol]],
                               colour = pop, shape = pop, size = pop)) +
            geom_point(aes(alpha = pop)) +
            scale_colour_manual(name = "Population", values = pca_colors) +
            scale_shape_manual(name = "Population", values = pca_shapes) +
            scale_size_manual(name = "Population", values = pca_sizes) +
            scale_alpha_manual(name = "Population", values = pca_alphas, guide = "none") +
            xlab(xlab) + ylab(ylab) +
            theme_ancestry +
            guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

        if (add_labels && nrow(study_pts) > 0) {
            p <- p + geom_text(data = study_pts,
                               aes(x = .data[[xcol]], y = .data[[ycol]],
                                   label = sample_id),
                               colour = "black", size = BASE_SIZE / 3, fontface = "bold",
                               vjust = -1.2, hjust = 0.5,
                               inherit.aes = FALSE)
        }
        return(p)
    }

    study_pts <- pca_all[pca_all$pop == "Study Sample", ]
    p1 <- make_pca_panel(pca_all, study_pts, "PC1", "PC2", "PC1", "PC2", TRUE)
    p2 <- make_pca_panel(pca_all, study_pts, "PC2", "PC3", "PC2", "PC3", FALSE)
    p3 <- make_pca_panel(pca_all, study_pts, "PC3", "PC4", "PC3", "PC4", FALSE)
    p4 <- make_pca_panel(pca_all, study_pts, "PC4", "PC5", "PC4", "PC5", FALSE)

    extract_legend <- function(p) {
        g <- ggplotGrob(p + theme(legend.position = "bottom"))$grobs
        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
        return(legend)
    }

    pca_legend <- extract_legend(p1)
    pca_lheight <- sum(pca_legend$height)

    pca_grob <- arrangeGrob(
        arrangeGrob(
            p1 + theme(legend.position = "none"),
            p2 + theme(legend.position = "none"),
            p3 + theme(legend.position = "none"),
            p4 + theme(legend.position = "none"),
            ncol = 2
        ),
        pca_legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - pca_lheight, pca_lheight),
        top = textGrob("Ancestry PCA \u2014 Study Samples on 1000 Genomes Reference",
                       gp = gpar(fontsize = BASE_SIZE + 4, fontface = "bold"),
                       vjust = 1)
    )

    save_grob(pca_grob, file.path(outdir, paste0(prefix, "_pca")), 14, 13)
} else {
    cat("PCA coordinates file not found — skipping PCA plots.\n")
    cat("  (This is expected if PLINK PCA failed due to high missingness.)\n\n")
}

# =============================================================================
# 2. COMBINED BARPLOT FIGURE
#    Layout (top to bottom):
#      [metadata strip(s)]         — optional annotation rows
#      [super-population barplot]  — K=5, no x-axis labels
#      [sub-population barplot]    — K=23, x-axis labels at bottom
#    If no sub-populations, super-pop barplot gets the x-axis labels.
# =============================================================================

cat("\nGenerating combined ancestry barplot...\n")

n_samples <- length(sample_order)

# -- Smart width scaling based on sample count --
# Small cohorts (<15):   wide bars, generous spacing
# Medium cohorts (15-50): moderate bars
# Large cohorts (50+):   compact bars, capped at 48 inches (under ggplot's 50in limit)
if (n_samples <= 15) {
    bar_width <- max(10, n_samples * 0.7 + 4)
} else if (n_samples <= 50) {
    bar_width <- n_samples * 0.45 + 4
} else {
    bar_width <- min(48, n_samples * 0.3 + 6)
}

# Scale font and bar parameters for large sample counts
if (n_samples > 60) {
    x_label_size <- max(5, BASE_SIZE - 6)
    bar_line_width <- 0.05
} else if (n_samples > 30) {
    x_label_size <- BASE_SIZE - 3
    bar_line_width <- 0.1
} else {
    x_label_size <- BASE_SIZE - 1
    bar_line_width <- 0.2
}

cat(sprintf("  Samples: %d | Figure width: %.1f in | Label size: %d pt\n",
            n_samples, bar_width, x_label_size))

pop_order <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# -- Super-population data (long format) --------------------------------------
sp_long <- reshape(superpops,
    varying   = names(superpops)[names(superpops) != "sample_id"],
    v.names   = "proportion",
    timevar   = "population",
    times     = names(superpops)[names(superpops) != "sample_id"],
    direction = "long",
    idvar     = "sample_id"
)
rownames(sp_long) <- NULL
sp_long$sample_id <- factor(sp_long$sample_id, levels = sample_order)
sp_long$population <- factor(sp_long$population, levels = rev(pop_order))

# -- Super-population barplot --------------------------------------------------
# If sub-pops exist, hide x-axis (shared axis at bottom); otherwise show it
p_superpop <- ggplot(sp_long, aes(x = sample_id, y = proportion, fill = population)) +
    geom_bar(stat = "identity", width = 0.92, colour = "white",
             linewidth = bar_line_width) +
    scale_fill_manual(name = "Super-population", values = superpop_colors,
                      breaks = pop_order) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    labs(x = NULL, y = "K=5") +
    theme_ancestry +
    theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.title.y = element_text(size = BASE_SIZE, face = "bold"),
        panel.grid.major.x = element_blank(),
        legend.position = "right"
    )

if (has_subpops) {
    # Hide x-axis labels on superpop — they'll appear on subpop below
    p_superpop <- p_superpop +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
} else {
    p_superpop <- p_superpop +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                         size = x_label_size, face = "bold"))
}

# -- Sub-population data and barplot (if available) ----------------------------
p_subpop <- NULL
subpop_order <- c()

if (has_subpops) {
    subpops <- read.delim(subpop_file, stringsAsFactors = FALSE)
    sub_cols <- names(subpops)[names(subpops) != "sample_id"]

    sub_long <- reshape(subpops,
        varying   = sub_cols,
        v.names   = "proportion",
        timevar   = "population",
        times     = sub_cols,
        direction = "long",
        idvar     = "sample_id"
    )
    rownames(sub_long) <- NULL
    sub_long$sample_id <- factor(sub_long$sample_id, levels = sample_order)

    for (sp in pop_order) {
        subs <- sort(names(subpop_to_superpop[subpop_to_superpop == sp]))
        subs <- subs[subs %in% sub_cols]
        subpop_order <- c(subpop_order, subs)
    }
    sub_long$population <- factor(sub_long$population, levels = rev(subpop_order))

    sub_colors_used <- subpop_colors[names(subpop_colors) %in% sub_cols]

    p_subpop <- ggplot(sub_long, aes(x = sample_id, y = proportion, fill = population)) +
        geom_bar(stat = "identity", width = 0.92, colour = "white",
                 linewidth = bar_line_width) +
        scale_fill_manual(name = "Sub-population", values = sub_colors_used,
                          breaks = subpop_order) +
        scale_x_discrete(drop = FALSE) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           expand = c(0, 0)) +
        labs(x = NULL, y = "K=23") +
        theme_ancestry +
        theme(
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            axis.title.y = element_text(size = BASE_SIZE, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1,
                                       size = x_label_size, face = "bold"),
            panel.grid.major.x = element_blank(),
            legend.position = "right",
            legend.key.size = unit(0.4, "cm"),
            legend.text = element_text(size = BASE_SIZE - 4)
        ) +
        guides(fill = guide_legend(ncol = 1))
}

# =============================================================================
# ASSEMBLE COMBINED FIGURE
# =============================================================================

cat("Assembling combined figure...\n")

# Collect all panels top-to-bottom
panels <- list()
panel_heights <- c()

# (a) Metadata strips
if (has_metadata) {
    strips <- build_meta_strips(meta_matched, meta_cols, sample_order)
    for (col in meta_cols) {
        panels <- c(panels, list(strips[[col]]))
        panel_heights <- c(panel_heights, 0.8)  # thin strip
    }
}

# (b) Super-population barplot
panels <- c(panels, list(p_superpop))
panel_heights <- c(panel_heights, 4)

# (c) Sub-population barplot
if (has_subpops && !is.null(p_subpop)) {
    panels <- c(panels, list(p_subpop))
    panel_heights <- c(panel_heights, 4.5)  # slightly taller for x-axis labels
}

# Use gtable to align all panels so x-axes line up perfectly
panel_grobs <- lapply(panels, ggplotGrob)

# Find the maximum width for each column across all grobs
max_widths <- panel_grobs[[1]]$widths
for (i in seq_along(panel_grobs)) {
    max_widths <- grid::unit.pmax(max_widths, panel_grobs[[i]]$widths)
}

# Apply uniform widths to all grobs
aligned_grobs <- lapply(panel_grobs, function(g) {
    g$widths <- max_widths
    return(g)
})

# Stack vertically
total_height <- sum(panel_heights)

combined <- arrangeGrob(
    grobs = aligned_grobs,
    ncol = 1,
    heights = unit(panel_heights, "in"),
    top = textGrob("Genetic Ancestry Estimation",
                   gp = gpar(fontsize = BASE_SIZE + 6, fontface = "bold"),
                   vjust = 0.5),
    padding = unit(0.5, "in")
)

# Add 1 inch for the title
total_height <- total_height + 1

save_grob(combined,
          file.path(outdir, paste0(prefix, "_barplot_combined")),
          bar_width, total_height)

# =============================================================================
# ALSO SAVE INDIVIDUAL BARPLOTS (for flexibility)
# =============================================================================

cat("\nSaving individual barplots...\n")

# Superpop standalone (with x-axis labels)
p_superpop_standalone <- p_superpop +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                     size = x_label_size, face = "bold")) +
    labs(title = "Continental Ancestry Proportions",
         subtitle = "ADMIXTURE K=5 supervised estimation")

save_plot(p_superpop_standalone,
          file.path(outdir, paste0(prefix, "_barplot_superpops")),
          width = bar_width, height = 6)

if (has_subpops && !is.null(p_subpop)) {
    p_subpop_standalone <- p_subpop +
        labs(title = "Sub-population Ancestry Proportions",
             subtitle = "ADMIXTURE K=23 supervised estimation")

    save_plot(p_subpop_standalone,
              file.path(outdir, paste0(prefix, "_barplot_subpops")),
              width = bar_width, height = 7)
}

# =============================================================================
# DONE
# =============================================================================

cat("\n============================================\n")
cat("Visualization complete!\n")
cat("============================================\n")
cat("Output files:\n")
if (has_pca) {
    cat(sprintf("  %s_pca.png / .svg                PCA projection (4 panels)\n", prefix))
}
cat(sprintf("  %s_barplot_combined.png / .svg    Combined ancestry figure\n", prefix))
cat(sprintf("  %s_barplot_superpops.png / .svg   Continental barplot (standalone)\n", prefix))
if (has_subpops) {
    cat(sprintf("  %s_barplot_subpops.png / .svg     Sub-population barplot (standalone)\n", prefix))
}
if (!has_pca) {
    cat("\n  Note: PCA plots were skipped (ancestry_pca_coordinates.tsv not found).\n")
    cat("  Re-run the pipeline with fixed PCA settings to generate PCA output.\n")
}
if (has_metadata) {
    cat(sprintf("\n  Metadata annotations: %s\n", paste(meta_cols, collapse = ", ")))
}
cat("============================================\n")
