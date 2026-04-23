################################################################################
# BICseq2
# https://math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/BICseq2.htmlc
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
  
)

# Install or load required packages
lapply(required_pkgs, bioc_package)

#load libraries

################################################################################
#hardcoded information
base_name <- "FD2500483"
data_path <- "/home/guest/internship/data/"
tool_path <- "/home/guest/internship/"
bam_file <- file.path(paste0(data_path, "FD2500483_sorted.bam"))
################################################################################
print(paste("Processing:", base_name))

## 1. get uniquelly mapped reads - samtools
#samfile <- file.path(paste0(data_path,"FD2500483_unique.sam"))
samfile <- file.path(paste0(data_path,"test.sam")) #smaller test file

system2(file.path(paste0(tool_path,"samtools-0.1.7a_getUnique-0.1.3/samtools")),    ##test !!!!!
        args = c("view", bam_file),
        stdout = samfile)

#file.exists(samfile)

#-------------------------------------------------
# make readPosFile
chr <- c()    #empty vectors
pos <- c()

for (line in readLines(samfile)){    #loop over reads
  parts <- strsplit(line, "\t")[[1]]
  chr <- c(chr, parts[3])
  pos <- c(pos, parts[4])
}

readpos_df <- data.frame(chr = chr, pos=pos)    #make df
readPos <- file.path(paste0(data_path,base_name, "_readPos.txt"))
write.table(    #write to file
  readpos_df,
  file=readPos,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)
#-------------------------------------------------
# MapFile
wigfile <- file.path(paste0(data_path, "hg38.fa.mappability_50bp.wig"))

## make mappability files from one wig file, only needed once !
chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY")
dir.create(file.path(paste0(data_path,"mapfiles")))



mapfiles <- c(paste0("chr", 1:22, ".map"), "chrX.map", "chrY.map")

#-------------------------------------------------
# make configFile
chromnames <- sort(unique(chr))
fafile <- file.path(paste0(data_path, "hg38.fa"))


binfilenorm <- file.path(paste0(chromnames,"_norm.bin"))

config <- data.frame(
  chromName = chromnames,
  faFile   = fafile,
  mapFile  = mapfiles,
  readPosFile  = readPos,
  binFileNorm  = binfilenorm
)

config_file <- file.path(paste0(data_path,base_name,"_config.txt"))
write.table(
  config,
  file = config_file,
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

# remove SAM file because very large file
#file.remove(samfile)

#-------------------------------------------------------------------------------
## 2. BICseq2-norm to remove the biases in the data
system2("chmod", args = c("777", "/home/guest/internship/BICseq2-norm_v0.2.6/tmp")) #add rights to tmp directory
system2(
  command = file.path(paste0(tool_path,"BICseq2-norm_v0.2.6/BICseq2-norm.pl")),
  args = c("-b", "50",config_file, "norm_output")
)


#-------------------------------------------------------------------------------
## 3. BICseq2-seg to detect CNVs based on the normalized data
cmd_seq <- 

print("Processing done")
################################################################################
print("Making visualizations")
print("Done")