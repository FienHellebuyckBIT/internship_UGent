################################################################################
# CNV Benchmarking Plots
# script made for visualisation of output file from cnv_benchmarking.py
################################################################################
# path to input file
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]


# Define the personal library path
personal_lib <- "/data/gent/510/vsc51018/R_libs"
dir <- dirname(input_file)


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
  "ggplot2","dplyr","tidyr","patchwork")

# Install or load required packages
lapply(required_pkgs, bioc_package)


# Load libraries

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Parse the results file 
parse_cnv_results <- function(filepath) {
  lines <- readLines(filepath)
  
  records <- list()
  i <- 1
  while (i <= length(lines)) {
    if (grepl("^sample:", lines[i])) {
      parts  <- strsplit(lines[i], "\t")[[1]]
      sample <- trimws(sub("sample:", "", parts[1]))
      tool   <- trimws(sub("tool:",   "", parts[2]))
      
      tp <- fp <- fn <- tn <- NA
      f1 <- fdr <- NA
      
      j <- i + 1
      while (j <= length(lines) && !grepl("^#{5}", lines[j])) {
        l <- trimws(lines[j])
        
        if (grepl("^P \\|", l) && is.na(tp)) {
          nums <- as.integer(regmatches(l, gregexpr("[0-9]+", l))[[1]])
          if (length(nums) == 2) { tp <- nums[1]; fp <- nums[2] }
        }
        if (grepl("^N \\|", l) && is.na(fn)) {
          nums <- as.integer(regmatches(l, gregexpr("[0-9]+", l))[[1]])
          if (length(nums) == 2) { fn <- nums[1]; tn <- nums[2] }
        }
        if (grepl("^F1 score:", l) && is.na(f1)) {
          f1 <- as.numeric(sub("F1 score:", "", l))
        }
        if (grepl("^FDR:", l) && is.na(fdr)) {
          fdr <- as.numeric(sub("FDR:", "", l))
        }
        j <- j + 1
      }
      
      tpr <- if (!is.na(tp) && !is.na(fn) && (tp + fn) > 0) tp / (tp + fn) else NA
      
      records[[length(records) + 1]] <- data.frame(
        sample = sample,
        tool   = tool,
        TP = tp, FP = fp, FN = fn, TN = tn,
        TPR = tpr,
        FDR = fdr,
        F1  = f1,
        stringsAsFactors = FALSE
      )
      i <- j
    } else {
      i <- i + 1
    }
  }
  bind_rows(records)
}

df <- parse_cnv_results(input_file)

# rename and order low to high coverage
df$coverage_label <- factor(
  df$sample,
  levels = c("sim_0.1x", "sim_0.5x", "sim_sorted"),
  labels = c("0.1x",     "0.5x",     "1x")
)

#print(df[, c("sample", "tool", "coverage_label", "TPR", "FDR", "F1")])

# colours for tools
# colorpalet from colorbrewer2.org BrBG (colorblind safe & print friendly)
tool_colors <- c(
  "HMMcopy"       = "#a6611a",   
  "Control-FREEC" = "#dfc27d",   
  "QDNAseq"       = "#80cdc1",  
  "CNVnator"      = "#018571"  
)


# Build plots
make_plot <- function(data, y_var, y_label) {
  ggplot(data, aes(x = coverage_label, y = .data[[y_var]],
                   colour = tool, group = tool)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5, shape = 16) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = tool_colors, name = "Tools") +
    labs(x = "Coverage (x)", y = y_label) +
    theme_gray()
}

pA <- make_plot(df, "TPR", "TPR")
pB <- make_plot(df, "FDR", "FDR")
pC <- make_plot(df, "F1",  "F1")

# Combine the plots
combined <- (pA + pB + pC) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") 


# Save the plots
#ggsave(file.path(dir,"cnv_benchmarking_plots.pdf"), combined, width = 10, height = 4)
ggsave(file.path(dir,"cnv_benchmarking_plots.png"), combined, width = 10, height = 4)

