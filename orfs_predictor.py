# ORFS_predictor.py
#---
# Written by Rodrigo Bacigalupe. MSc in Bioinformatics student.
# University of Edinburgh.
#---
# Calculate the number of complete and fragmented genes in assemblies. It takes as input the predict
# file produced by GlimmerMG using phymm.

#-------------------------------------------------------------------------------------------------
# Modules to deal with the file system, parse arguments and use command line arguments
import os, argparse, sys

#----------------------------------- PARSE ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='ORFs predictor.')
parser.add_argument("-p", "--predict", help="Introduce the name of the predict file")
args = parser.parse_args()
# Variable that stores the predict file
predict_file = args.predict

#----------------------------------- CALCULATE ORFS ----------------------------------------
# Define a function that takes the predict file as input and returns a summary table

def quality_stats(predict_file):
    total_orfs = 0
    full_genes = 0
    orf_gene = 0
    fragmented_genes = 0
    
    # Break the file
    for line in open(predict_file):
        if line.startswith('>'):
            break

    for line in open(predict_file):
        if line.startswith('>'):
            contig_length = 0
            contig,contig_length = line.rstrip().split('length=_')

        # Calculate total ORFS, full genes and fragmented genes
        if line.startswith('orf'):
            total_orfs = total_orfs + 1
            orfdata = line.rsplit()
            if int(orfdata[1]) in (range(-3,4)) or int(orfdata[1]) in range(int(contig_length)-3,int(contig_length)+4):
                if int(orfdata[2]) in (range(-3,4)) or int(orfdata[2]) in range(int(contig_length)-3,int(contig_length)+4):
                    orf_gene = orf_gene + 1

                if int(orfdata[2]) not in (range(-3,4)) and int(orfdata[2]) not in range(int(contig_length)-3,int(contig_length)+4):
                    fragmented_genes = fragmented_genes + 1

            if int(orfdata[2]) in (range(-3,4)) or int(orfdata[2]) in range(int(contig_length)-3,int(contig_length)+4):
                if int(orfdata[1]) not in (range(-3,4)) and int(orfdata[1]) not in range(int(contig_length)-3,int(contig_length)+4):
                    fragmented_genes = fragmented_genes + 1

            if int(orfdata[1]) not in (range(-3,4)) and int(orfdata[1]) not in range(int(contig_length)-3,int(contig_length)+4):
                if int(orfdata[2]) not in (range(-3,4)) and int(orfdata[2]) not in range(int(contig_length)-3,int(contig_length)+4):
                    full_genes = full_genes + 1

    if int(total_orfs) == int(full_genes)+int(orf_gene)+int(fragmented_genes):
        return total_orfs, full_genes, orf_gene, fragmented_genes

# Compute the numbers of ORFs/genes
total_orfs, full_genes, orf_gene, fragmented_genes = quality_stats(predict_file)

#--------------------------------------- PRINT THE RESULTS -----------------------------------------
results = ("Total ORFs	Full genes	ORFs=gene	fragmented genes"+'\n'+str(total_orfs)+"	"+str(full_genes)+"	"+str(orf_gene)+"	"+str(fragmented_genes))
print(results)

#----------------------------------- SAVE THE OUTPUT IN A FILE --------------------------------------
open('./'+predict_file+"_"+"stats", 'wt').write(results)
