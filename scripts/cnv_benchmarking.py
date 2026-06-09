#!/usr/bin/python3
##################################################################################################################
# True Positive Rate = The number of bins that overlap 60% or more with an actual CNV and is predicted correctly as gain or loss.
# False Positive Rate = The number of bins that overlap less then 60% or not with an actual CNV or that is predicted wrong as gain or loss.
# False Negative Rate = The number of bins that are not detected as a CNV but overlap 60% or more with an actual CNV.
# True Negative Rate = The number of bins that are not detected as a CNV and overlap less then 60% or not with an actual CNV.
##################################################################################################################

import csv
import sys
import os

# fetch command line arguments
input_path = sys.argv[1]
actual_data = sys.argv[2]
outfile = sys.argv[3]


# define statistics:
TPR = 0
FPR = 0
FNR = 0
TNR = 0


# Open file with actual simulated CNVs data.
with open(actual_data, newline="", encoding="utf-8") as actual_cnvdata:
    actual_cnv_rows = list(csv.DictReader(actual_cnvdata, delimiter="\t"))

    # Open the input file with predicted CNVs
    with open(input_path, newline="", encoding="utf-8") as infile:
         # Use DictReader to read each row as a dict.
        predicted_cnv_rows = list(csv.DictReader(infile, delimiter="\t"))
      

        #loop over rows in both files and compare the coordinates and copy number to calculate TPR and TNR.
        for predicted_cnv_row in predicted_cnv_rows:
            #reset 
            overlaps_actual = False
            for actual_cnv_row in actual_cnv_rows:

                # Skip rows with missing coordinates
                if not actual_cnv_row.get("start") or not actual_cnv_row.get("end") or not predicted_cnv_row.get("start") or not predicted_cnv_row.get("end"):
                    continue
                
                # Convert coordinates to integers
                start_actual_cnv = int(float(actual_cnv_row["start"]))
                end_actual_cnv = int(float(actual_cnv_row["end"]))
                start_predicted_cnv = int(float(predicted_cnv_row["start"]))
                end_predicted_cnv = int(float(predicted_cnv_row["end"]))

                # compute overlap length between bin and CNV 
                overlap_start = max(start_actual_cnv, start_predicted_cnv)
                overlap_end = min(end_actual_cnv, end_predicted_cnv)
                overlap_len = overlap_end - overlap_start

                # treshold = 60% of the bin
                threshold = 0.5 * (end_predicted_cnv - start_predicted_cnv)

                # does bin overlap with cnv?
                if ((overlap_len >= threshold) and
                    # check chromosome matches
                    (actual_cnv_row.get("chr") == predicted_cnv_row.get("chr").replace("chr", ""))):
                    overlaps_actual = True
                    break
                
            # assign per bin if it is neutral gain or loss based on log2ratio
            copynumber=float(predicted_cnv_row.get("copynumber"))
            if copynumber >= 0.58:
                CNV="Gain"
            elif copynumber <= -1:
                CNV="Loss"
            else: 
                CNV="Neutral"
                    
            ########## check for FNR and TNR ##########
            if CNV == "Neutral": 
                if overlaps_actual==True:
                    FNR+=1
                else:
                    TNR+=1
                        

            ########## check for TPR and FPR ##########
            # check predicted_cnv is detected
            if CNV != "Neutral":
                if ((float(actual_cnv_row.get("copynumber")) > 2 and CNV == "Gain" ) or 
                (float(actual_cnv_row.get("copynumber")) < 2 and CNV == "Loss" ) and
                overlaps_actual==True
                ):
                    TPR+=1
                else:
                    FPR+=1
                        
              
########## calculate metrics and write to file ########## 
# calculate recall and precision
precision = TPR / (TPR + FPR) 
recall = TPR / (TPR + FNR)

# calculate F1 score
F1= 2 * (precision * recall) / (precision + recall)

#false discovery rate
FDR=(FPR/(FPR+TPR))

# calculate results in percentage
TPRp = (recall*100)
FPRp = ((FPR/(FPR+TNR))*100)
FNRp = ((FNR/(FNR+TPR))*100)
TNRp = ((TNR/(TNR+FPR))*100)
F1p = F1*100
FDRp = FDR*100

# open or create an output file in append mode
with open(outfile, "a", encoding="utf-8") as f:
    #append tool and sample name
    base = os.path.splitext(os.path.basename(input_path))[0]
    parts = base.split("_")
    sample = "_".join(parts[0:2]) 
    tool = parts[-4]
    f.write(f"\nsample: {sample}\ttool: {tool}\n")

    #append metrics
    f.write("\nResults:\n")
    f.write("-" * 32 + "\n")
    f.write("Confusion matrix (counts):\n")
    f.write("-" * 32 + "\n")
    f.write("Predicted (vertical) vs actual (horizontal) CNVs\n")
    f.write("-" * 32 + "\n")
    f.write("   |    P    |    N    |\n")
    f.write("-" * 32 + "\n")
    f.write(f" P | {TPR:7d} | {FPR:7d} |\n")
    f.write(f" N | {FNR:7d} | {TNR:7d} |\n")
    f.write("-" * 32 + "\n")
    f.write(f"F1 score: {F1:.4f}\nFDR: {FDR:.4f}\n")
    #append percentages metrics
    f.write("\nResults in percentage:\n")
    f.write("-" * 32 + "\n")
    f.write("Confusion matrix (percentages):\n")
    f.write("-" * 32 + "\n")
    f.write("Predicted (vertical) vs actual (horizontal) CNVs\n")
    f.write("-" * 32 + "\n")
    f.write("   |    P    |    N    |\n")
    f.write("-" * 32 + "\n")
    f.write(f" P | {TPRp:7.2f} | {FPRp:7.2f} |\n")
    f.write(f" N | {FNRp:7.2f} | {TNRp:7.2f} |\n")
    f.write("-" * 32 + "\n")
    f.write(f"F1 score: {F1p:.2f}\nFDR: {FDRp:.2f}\n")
    f.write("#"*50)


