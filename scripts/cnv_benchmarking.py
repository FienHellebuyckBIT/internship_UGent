#!/usr/bin/python3
##################################################################################################################
# True Positive Rate = The number bins that are correctly classified as a CNV gain or loss by the tool.
# False Positive Rate = The number of bins that were predicted as a CNV by the tool but were not actually a CNV in the data.
# False Negative Rate = The number of bins where one actual CNV overlaps that bin at least half (50kb for bins of 100kb)
# True Negative Rate = The number of bins that were correctly classified as neutral by the tool.
##################################################################################################################

import csv
import sys
import os

# fetch command line arguments
input_path = sys.argv[1]
actual_data = sys.argv[2]
outfile = sys.argv[3]

# determine treshold
bin_size = 100000
treshold = bin_size/2

# define statistics:
TPR = 0
FPR = 0
FNR = 0
TNR = 0
count_predicted_bins =0

#define matched_predicted, set to true if a match is found so it is not counted multiple times.
matched_predicted = False

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
            matched_predicted = False
            overlaps_actual = False
            #count nr of rows
            count_predicted_bins +=1
            for actual_cnv_row in actual_cnv_rows:

                # Skip rows with missing coordinates
                if not actual_cnv_row.get("start") or not actual_cnv_row.get("end") or not predicted_cnv_row.get("start") or not predicted_cnv_row.get("end"):
                    continue
                
                # Convert coordinates to integers
                start_actual_cnv = int(float(actual_cnv_row["start"]))
                end_actual_cnv = int(float(actual_cnv_row["end"]))
                start_predicted_cnv = int(float(predicted_cnv_row["start"]))
                end_predicted_cnv = int(float(predicted_cnv_row["end"]))


                ########## check for TNR and FNR ##########
                if predicted_cnv_row.get("color_group") == "Neutral":

                    # compute overlap length 
                    overlap_start = max(start_actual_cnv, start_predicted_cnv)
                    overlap_end = min(end_actual_cnv, end_predicted_cnv)
                    overlap_len = overlap_end - overlap_start

                    # only count as overlap if >= half bin
                    if ((overlap_len >= treshold) and
                        (actual_cnv_row.get("chr") == predicted_cnv_row.get("chr").replace("chr", ""))):
                        overlaps_actual = True

                    continue


                ########## check for TPR and FPR ##########
                # check predicted_cnv is detected
                if (predicted_cnv_row.get("color_group") != "Neutral"):
                    # check chromosome matches
                    if (actual_cnv_row.get("chr") == predicted_cnv_row.get("chr").replace("chr", "") and
                        # check copynumber gain or loss is right
                        ((float(actual_cnv_row.get("copynumber")) > 2 and predicted_cnv_row.get("color_group") == "Gain" ) or 
                        (float(actual_cnv_row.get("copynumber")) < 2 and predicted_cnv_row.get("color_group") == "Loss" )) and
                        (
                            # the actual_cnv is inbetween or equal to the coordinates of the bin that is correctly classified as cnv
                            (start_predicted_cnv >= start_actual_cnv <= end_predicted_cnv and start_predicted_cnv >= end_actual_cnv <= end_predicted_cnv) or

                            #the predicted_cnv is inbetween the actual_cnv coordinates
                            (start_actual_cnv > start_predicted_cnv < end_actual_cnv and start_actual_cnv > start_predicted_cnv < end_actual_cnv ) or
                            
                            # the actual start position is inbetween the predicted_cnv coordinates and the actual end position is equal or higher.
                            (start_predicted_cnv >= start_actual_cnv <= end_predicted_cnv and end_actual_cnv > start_predicted_cnv) or

                            # the actual start position is equal to or smaller then the predicted start and the end position is inbetween the predicted_cnv coordinates
                            (start_actual_cnv < start_predicted_cnv and start_predicted_cnv >= end_actual_cnv <= end_predicted_cnv)
                        )
                        ): 
                        # only match predicted CNV once
                        if not matched_predicted:
                            TPR += 1
                            matched_predicted = True
                            #print("CNV: {} and SIM:{} ".format(predicted_cnv_row, actual_cnv_row))
                            break
            # false positives        
            if not matched_predicted and predicted_cnv_row.get("color_group") != "Neutral":
                FPR += 1 
            # true and false negatives
            if predicted_cnv_row.get("color_group") == "Neutral":
                if overlaps_actual:
                    FNR += 1   
                else:
                    TNR += 1  
               
                    
# calculate recall and precision
precision = TPR / (TPR + FPR) 
recall = TPR / (TPR + FNR)

# calculate F1 score
F1= 2 * (precision * recall) / (precision + recall)

#false discovery rate
FDR=(FPR/(FPR+TPR))

# calculate results in percentage
TPRp = int(recall*100)
FPRp = int((FPR/(FPR+TNR))*100)
FNRp = int((FNR/(FNR+TPR))*100)
TNRp = int((TNR/(TNR+FPR))*100)
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
    f.write(f" P | {TPRp:7d} | {FPRp:7d} |\n")
    f.write(f" N | {FNRp:7d} | {TNRp:7d} |\n")
    f.write("-" * 32 + "\n")
    f.write(f"F1 score: {F1p:.2f}\nFDR: {FDRp:.2f}\n")
    f.write("#"*50)

