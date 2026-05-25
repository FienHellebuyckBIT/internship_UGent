# CNVnator
################################################################################
# Define the personal library path
personal_lib <- "/data/gent/510/vsc51018/R_libs"

# Add personal library to .libPaths
.libPaths(c(personal_lib, .libPaths()))

#installing required packages if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = personal_lib)
}

bioc_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, lib = personal_lib, ask = FALSE, update = FALSE, type="source")
  }
  library(pkg, character.only = TRUE)
}

# List of required packages
required_pkgs <- c(
  "plotly","dplyr","tidyr","htmlwidgets","GenomicRanges","stringr"
)

# Install or load required packages
lapply(required_pkgs, bioc_package)

#load libraries

library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(GenomicRanges)
library(stringr)

################################################################################
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
#paths
base_name <- args[1]
outdir <- args[2]
BIN_SIZE <- args[3]
tsv_file <- args[4]
calls_file <- args[5]
################################################################################
### pre-processing ###
# load per bin data from rootfile
bins <- read.table(
  tsv_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# load cnv calls
calls <- read.table(
  calls_file,
  fill = TRUE,
  stringsAsFactors = FALSE
)



#Parse CNVnator call regions
calls$type  <- calls$V1
calls$chr   <- str_extract(calls$V2, "^chr[^:]+")
calls$start <- as.numeric(str_extract(calls$V2, "(?<=:)[0-9]+"))
calls$end   <- as.numeric(str_extract(calls$V2, "(?<=-)[0-9]+"))



#Add segmentation column to bins
# Default: normal copy number
bins$segmented <- 0


#Assign segmentation values per bin
for (i in 1:nrow(calls)) {
  
  idx <- bins$chr == calls$chr[i] &
    bins$start >= calls$start[i] &
    bins$end   <= calls$end[i]
  
  if (calls$type[i] == "deletion") {
    bins$segmented[idx] <- -1
  }
  
  if (calls$type[i] == "duplication") {
    bins$segmented[idx] <- 1
  }
}



#Create the dataframe

cnv_df <- data.frame(
  chr        = bins$chr,
  start      = bins$start,
  end        = bins$end,
  position   = bins$position,
  copynumber = bins$log2ratio,
  segmented  = bins$segmented
)

### Preview
head(cnv_df)


##############################################
### visualisations ###
# Remove rows with NA values
clean_cnv_df <- cnv_df %>% drop_na()



# Chromosome order
present_chromosomes <- unique(clean_cnv_df$chr)
chromosome_order <- paste0("chr", c(1:22, "X", "Y"))
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
total_reads <- sum(bins$rd, na.rm = TRUE)


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
  file= file.path(outdir ,paste0(base_name,"_", BIN_SIZE, "_CNV_FREEC.html")),
  selfcontained = TRUE,
  libdir = NULL
)


