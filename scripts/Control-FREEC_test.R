################################################################################
# Control-FREEC
# https://boevalab.inf.ethz.ch/FREEC/tutorial.html#CONFIG
################################################################################
#installing required packages if not installed
library(BiocManager)

bioc_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}
# List of required packages
required_pkgs <- c(
  "glue","plotly","dplyr","tidyr","htmlwidgets","GenomicRanges"
)

# Install or load required packages
lapply(required_pkgs, bioc_package)

#load libraries

library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(GenomicRanges)
library(glue)
################################################################################
#paths
base_name <- "FD2500483"
data_path <- "/home/guest/internship/data/"
tool_path <- "/home/guest/internship/"
bam_file <- paste0(data_path, "FD2500483_sorted.bam")
################################################################################
windowsize <- 50000
# Config File
config <- glue(
  "[general]
  chrLenFile = {data_path}hg38.fa.fai 
  ploidy = 2 
  gemMappabilityFile = {data_path}out100m2_hg38.gem
  maxThreads = 2 
  outputDir = {data_path}
  uniqueMatch = TRUE 
  window = {windowsize} 
  chrFiles = {data_path}chrFiles_freec
  [sample]
  mateFile = {bam_file}
  inputFormat = BAM
  mateOrientation = FR 
  [control]
")

configFile <- paste0(data_path, "freec_config.txt")
writeLines(config, configFile)

# control-FREEC command
#PathToFREEC/freec -conf myConfig.txt -sample sample.bam -control control.bam
freec <- paste0(tool_path,"FREEC-11.6b/src/freec")
system2(
  command = freec,
  args = c(
    "-conf", configFile,
    "-sample", bam_file
  )
)

################################################################################
# ---------------------------
# New visualization
# ---------------------------

# read ratio file
ratio_data <- read.table(paste0(data_path,"FD2500483_sorted.bam_ratio.txt"), header = TRUE)

# Replace -1 with NA before log2 transformation
ratio_data$Ratio[ratio_data$Ratio == -1] <- NA
ratio_data$MedianRatio[ratio_data$MedianRatio == -1] <- NA

# Create the dataframe
cnv_df <- data.frame(
  chr        = ratio_data$Chromosome,
  start      = ratio_data$Start,
  end        = ratio_data$Start + windowsize-1,
  position   = ratio_data$Start + ((windowsize-1)/2) ,
  copynumber = log2(ratio_data$Ratio),
  segmented  = log2(ratio_data$MedianRatio)
)

# Remove rows with NA values
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

# read .cnp file for read counts
cnp_data <- read.table(paste0(data_path,"FD2500483_sorted.bam_sample.cpn"), header = FALSE)
# sum all read starts (column 3)
total_reads <- sum(cnp_data[, 3], na.rm = TRUE)


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
  file= paste0(base_name, "_CNV_C-FREEC.html"), 
  selfcontained = TRUE,
  libdir = NULL
)




print(paste("Process of", base_name, "completed."))
