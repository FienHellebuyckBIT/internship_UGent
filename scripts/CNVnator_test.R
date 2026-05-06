# CNVnator
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
#paths
base_name <- "FD2500483"
data_path <- "/home/guest/internship/data/"
#tool_path <- "/home/guest/internship/"
bam_file <- paste0(data_path, "FD2500483_sorted.bam")
setwd(data_path)
################################################################################
#commands
# conda activate cnvnator
# 1. Extract read mapping
# system("cnvnator -root FD2500483_cnvnator.root -tree FD2500483_sorted.bam -chrom $(seq -f 'chr%g' 1 22) chrX chrY -lite")
# # 2. Generate histogram
# system("cnvnator -root FD2500483_cnvnator.root -his 1000 -d chrFiles_freec -chrom $(seq -f 'chr%g' 1 22) chrX chrY")
# # 3. Calculate statistics
# system("cnvnator -root FD2500483_cnvnator.root -stat 1000 -chrom $(seq -f 'chr%g' 1 22) chrX chrY")
# # 4 RD SIGNAL PARTITIONING
# system("cnvnator -root FD2500483_cnvnator.root -chrom $(seq -f 'chr%g' 1 22) chrX chrY -partition 1000 ")
# # 5 CNV calling
# system("cnvnator -root FD2500483_cnvnator.root -chrom $(seq -f 'chr%g' 1 22) chrX chrY -call 1000 > FD2500483_cnvnator_calls.txt")

################################################################################
#visualisation
#extract info from calls.txt file
cnv_calls <- read.table("FD2500483_cnvnator_calls.txt",sep = "\t", header = FALSE)
#extract coordinates
chr <- str_extract(cnv_calls$V2, "^chr[^:]+")
start_end <- str_extract(cnv_calls$V2, "[0-9]+-[0-9]+$")
start <- as.numeric(str_extract(start_end, "^[0-9]+"))
end <- as.numeric(str_extract(start_end, "[0-9]+$"))
#extract normalized readdepth
rd <- as.numeric(cnv_calls$V4)
copynumber_log2 <- log2(pmax(rd, 1e-6))


#Create the dataframe
cnv_df <- data.frame(
  chr        = chr,
  start      = start,
  end        = end,
  position   = ((start+end)/2),
  copynumber = copynumber_log2,
  segmented  = copynumber_log2
  
)



##############################################



###############################################

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
  file= paste0(base_name, "_cnvnator.html"), 
  selfcontained = TRUE,
  libdir = NULL
)


