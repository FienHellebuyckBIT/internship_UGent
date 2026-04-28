# internship_UGent
During my internship I evaluated four tools for detecting Copy Number Variations (CNVs):
* HMMcopy
* QDNAseq
* Control‑FREEC
* BICseq2
All scripts used in this project are available in the scripts directory.

1. *_test.R
These scripts were used to explore each tool and run it on a single sample.
Most paths and the sample file are hardcoded, so you’ll need to adjust them to match your own folder structure and data.
All scripts in this group use the human genome (hg38).
The BICseq2 script is not finished yet.
2. *_mouse.R
These scripts are a bit more automated. Only a few paths need to be updated at the top of the file.
They are based on the mouse genome (mm10).

