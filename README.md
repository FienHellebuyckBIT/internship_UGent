# internship_UGent
During my internship I evaluated different tools for detecting Copy Number Variations (CNVs):
* HMMcopy
* QDNAseq
* Control‑FREEC
* BICseq2
* CNVnator
  
All scripts used in this project are available in the scripts directory.

1. *_test.R
   
These scripts were used to explore each tool and run it on a single sample.
Most paths and the sample file are hardcoded, so you’ll need to adjust them to match your own folder structure and data.
All scripts in this group use the human genome (hg38).
The scripts use the same visualisation approach, which makes it easier to compare the tools.
Note that the BICseq2 and CNVnator scripts are not finished yet.

2. *_mouse.R
   
These scripts are a bit more automated. Only a few paths need to be updated at the top of the file.
They are based on the mouse genome (mm10).
Some commands that previously had to be run manually (to create required input files before running the tools) are now integrated directly into the scripts.

Documentation webpages for these tools:
* HMMcopy: https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
* QDNAseq: https://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html
* Control-FREEC: https://boevalab.inf.ethz.ch/FREEC/
* BICseq2: https://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/BICseq2.html


The project data can be requested by contacting me at: fien.hellebuyck@student.howest.be
