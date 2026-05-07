###############################################################################
# HMMcopy for multiple mouse samples, adjusted for use on VSC
###############################################################################
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
    BiocManager::install(pkg, lib = personal_lib)
  }
  library(pkg, character.only = TRUE)
}

# List of required packages
required_pkgs <- c(
  "HMMcopy","plotly","dplyr","tidyr","htmlwidgets","GenomicRanges"
)

# Install or load required packages
lapply(required_pkgs, bioc_package)

#load libraries
library(HMMcopy)
library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(GenomicRanges)


################################################################################
# Bam file path
bam_file <- args[1]
# sample name
base_name <- args[2]
# directory of sample
sample_dir <- dirname(bam_file)
# main directory
mouse_data_dir <- dirname(sample_dir)

################################################################################

#define files
rfile <- file.path(sample_dir, paste0("read_", base_name, ".wig"))
gfile <- file.path(mouse_data_dir, "gc_mm10.wig")
mfile <- file.path(mouse_data_dir, "map_mm10.wig")

################################################################################
#hmmcopy tools:
hmmcopy_tools <- Sys.getenv("HMMCOPY_TOOLS")

################################################################################
  #start
  print(paste("Processing:", base_name))
  
  #make read file
   if (!file.exists(rfile)){
     readCounter <- file.path(hmmcopy_tools, "readCounter")
   system(
     paste0(
       readCounter,
       " -w 1000000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,",
       "chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY ",
       bam_file," > ", file.path(sample_dir, paste0("read_", base_name, ".wig"))
     )
   )}

  ##############################################################################
  
  reads <- wigsToRangedData(rfile, gfile, mfile)
  # Correct reads into copy number
  corrected_readcount <- correctReadcount(reads)
  # segment the read counts
  segments <- HMMsegment(corrected_readcount)
  

  #----------------------
  # New visualization    
  #----------------------
  
  # Extract segment table
  df_segs <- segments$segs
  
  # Convert segments to GRanges
  gr_segs <- GRanges(
    seqnames = df_segs$chr,
    ranges   = IRanges(start = df_segs$start, end = df_segs$end),
    median   = df_segs$median
  )
  
  # Convert bins (corrected_readcount) to GRanges
  gr_bins <- GRanges(
    seqnames = corrected_readcount$chr,
    ranges   = IRanges(start = corrected_readcount$start,
                       end   = corrected_readcount$end)
  )
  
  # Find overlaps between bin and segment
  hits <- findOverlaps(gr_bins, gr_segs)
  
  # Create a per-bin segmented vector
  segmented_vec <- rep(NA, length(gr_bins))
  segmented_vec[queryHits(hits)] <- mcols(gr_segs)$median[subjectHits(hits)]
  
  # Make df
  cnv_df <- data.frame(
    chr        = corrected_readcount$chr,
    start      = corrected_readcount$start,
    end        = corrected_readcount$end,
    position   = ((corrected_readcount$start + corrected_readcount$end) / 2),
    copynumber = corrected_readcount$copy,   # already log2
    segmented  = segmented_vec               # per-bin segmented values
  )
  
  
  # Clean NAs
  clean_cnv_df <- cnv_df %>% drop_na()
  
  # Chromosome order
  present_chromosomes <- unique(clean_cnv_df$chr)
  chromosome_order <- c(paste0("chr", 1:19), "chrX", "chrY")
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
  total_reads <- sum(corrected_readcount$reads, na.rm = TRUE)
  
  
  
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
    file= file.path(sample_dir, paste0(base_name, "_CNV_HMMcopy_mouse.html")), 
    selfcontained = TRUE,
    libdir = NULL
  )


  print(paste("Process of", base_name, "completed."))



print("Done")
