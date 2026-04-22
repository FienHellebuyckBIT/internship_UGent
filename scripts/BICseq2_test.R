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
  
)

# Install or load required packages
lapply(required_pkgs, bioc_package)

#load libraries

################################################################################
#hardcoded information
base_name <- "FD2500483"
data_path <- "/home/guest/internship/data/"
bam_file <- paste0(data_path, "FD2500483_sorted.bam")
################################################################################
print(paste("Processing:", base_name))

# 1. get uniquelly mapped reads - samtools
input  <- file.path(bam_file)
output <- file.path(data_path, "FD2500483_unique.bam")

cmd_samtools <- paste(
  "samtools view -b -q 0 -F 4 -F 256",
  input,
  "-o",
  output
)

system(cmd_samtools)

----------------------------------------------------------
# 2. BICseq2-norm to remove the biases in the data
# install BICseq2-norm
#---------------------------------
# if BICseq2-norm is not installed:
#$ conda create -n bicseq2
#$ conda activate bicseq2
#$ conda install bioconda::bicseq2-norm
#--------------------------------- 
system("conda init")
system("conda activate bicseq2") #problems try maybe reinstal the tool?
system("bicseq2-norm --version")
# command BICseq2-norm
cmd_norm <- 
system()


  
# 3. BICseq2-seg to detect CNVs based on the normalized data
cmd_seq <- 

print("Processing done")
################################################################################
print("Making visualizations")
print("Done")