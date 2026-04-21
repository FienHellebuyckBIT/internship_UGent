# Define the personal library path
#personal_lib <- "./R_libs"

# Add personal library to .libPaths
#.libPaths(c(personal_lib, .libPaths()))

# Install BiocManager if needed and load it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
library(BiocManager)

bioc_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

# List of required packages
required_pkgs <- c(
  "remotes", "Cairo", "GenomicRanges", "QDNAseq", "ggplot2", "dplyr", "tidyr", "Biobase", "plotly", "htmlwidgets"
)

# Install or load required packages
lapply(required_pkgs, bioc_package)

# Install and load QDNAseq.hg38 from GitHub if not present
if (!requireNamespace("QDNAseq.hg38", quietly = TRUE)) {
  remotes::install_github("asntech/QDNAseq.hg38@main")
}
#Load the required libraries
library(remotes)
library(GenomicRanges)
library(QDNAseq)ls 
library(Biobase)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)
if (!requireNamespace("QDNAseq.hg38", quietly = TRUE)) {
  remotes::install_github("asntech/QDNAseq.hg38@main")
}
library(QDNAseq.hg38)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

#Base name
base_name <- args[1]

# Bam file path
bam_file <- args[2]

# Bin size
bin_size <-  as.numeric(args[3])

# Bin annotation
bins <- getBinAnnotations(binSize = bin_size, genome = "hg38")

#Output directory
output <- dirname(bam_file)
#output <- file.path(output, paste0(base_name, "_results"))

print(paste("Processing:", base_name))

readCounts <- binReadCounts(bins, bam_file)

readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)

readCountsFiltered <- estimateCorrection(readCountsFiltered)

#Include Sex Chromosomes
readCountsFiltered <- applyFilters(readCountsFiltered, residual = TRUE, blacklist = TRUE, chromosomes = NA)

copyNumbers <- correctBins(readCountsFiltered)

copyNumbersNormalized <- normalizeBins(copyNumbers)

copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)


# Call bins for CNV calling
copyNumbersCalled <- callBins(copyNumbersSegmented)



# ---------------------------
# New visualization
# ---------------------------

# Extract feature coordinates
features <- fData(copyNumbersCalled)

# Extract CNV values
copynumber_values <- copyNumbersCalled@assayData$copynumber[,1]   # raw CNV
segmented_values  <- copyNumbersCalled@assayData$segmented[,1]    # segmented CNV

# Combine into a single data frame
cnv_df <- data.frame(
    chr        = features$chromosome,
    start      = features$start,
    end        = features$end,
    position   = ((features$start + features$end) / 2),
    copynumber = log2(copynumber_values),
    segmented  = log2(segmented_values)
    )


# Clean NAs
clean_cnv_df <- cnv_df %>% drop_na()

# Chromosome order
present_chromosomes <- unique(clean_cnv_df$chr)
chromosome_order <- c(as.character(1:23), "X", "Y")
chromosome_order <- chromosome_order[chromosome_order %in% present_chromosomes]
clean_cnv_df$chr <- factor(clean_cnv_df$chr, levels = chromosome_order)

# Chromosome boundaries for X-axis
chr_boundaries <- aggregate(position ~ chr, data = clean_cnv_df, max)
chr_boundaries <- chr_boundaries[order(match(chr_boundaries$chr, chromosome_order)), ]
chr_offset <- c(0, cumsum(chr_boundaries$position)[-nrow(chr_boundaries)])
names(chr_offset) <- chromosome_order
clean_cnv_df$cumulative_position <- clean_cnv_df$position + chr_offset[clean_cnv_df$chr]
chr_boundaries$cumulative_position <- chr_boundaries$position + chr_offset[chr_boundaries$chr]
chr_boundaries$midpoint <- (c(0, head(chr_boundaries$cumulative_position, -1)) + chr_boundaries$cumulative_position) / 2

# Threshold for gain/loss
threshold <- 0.35
clean_cnv_df$color_group <- cut(clean_cnv_df$segmented, breaks = c(-Inf, -threshold, threshold, Inf),
                                  labels = c("Loss", "Neutral", "Gain"))
clean_cnv_df$cbs_color_group <- cut(clean_cnv_df$segmented, breaks = c(-Inf, -threshold, threshold, Inf),
                                      labels = c("CBS Loss", "CBS Neutral", "CBS Gain"))

# Number of reads
total_reads <- sum(assayData(readCountsFiltered)$counts, na.rm = TRUE)


# plotly

pal <- c("red","black","blue")
pal <- setNames(pal, c("Loss","Neutral","Gain"))

plotly_plot <- plot_ly(
  data = clean_cnv_df,
  x = ~cumulative_position,
  y = ~copynumber,
  type = "scatter",
  mode = "markers",
  color = ~color_group,
  colors = pal,
  marker = list(
    size = 4,       # size of the CNV dots
    opacity = 0.5   # decrease opacity
  )
) %>%
  # Add CBS points as a separate trace
  add_markers(
    y = ~segmented,
    x = ~cumulative_position,
    marker = list(
      color = "orange",
      size = 4,
      opacity = 0.8
    ),
    showlegend = FALSE
  ) %>%
  layout(
    title = list(text = paste("CNV Profile for Sample", base_name, "(", format(total_reads, big.mark = ","), "reads )"), x = 0.5),
    xaxis = list(
      title = list(text = "Chromosome", font = list(size = 14, family = "Arial Black")),
      tickvals = chr_boundaries$midpoint,
      ticktext = chr_boundaries$chr,
      showgrid = FALSE,
      gridcolor = "#d3d3d3",
      gridwidth = 0.5,
      tickfont = list(size = 12, family = "Arial Black")
    ),
    yaxis = list(
      title = list(text = "Log2 Ratio", font = list(size = 14, family = "Arial Black")),
      range = c(-2, 2),
      showgrid = TRUE,
      gridcolor = "#d3d3d3",
      gridwidth = 0.5
    ),
    shapes = lapply(chr_boundaries$cumulative_position, function(x) {
      list(
        type = "line",
        x0 = x, x1 = x,
        y0 = -2, y1 = 2,
        line = list(color = "black", width = 0.5, dash = "solid")
      )
    }),
    showlegend = FALSE
  )


saveWidget(
  plotly_plot,
  file= file.path(output,  paste0(base_name, "_CNV.html")),
  selfcontained = TRUE,
  libdir = NULL
)




print(paste("Process of", base_name, "completed."))