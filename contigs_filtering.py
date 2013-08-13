# contigs_filtering.py
#---
# Written by Rodrigo Bacigalupe. MSc in Bioinformatics student.
# University of Edinburgh.
#---
# Filter contigs shorter than 300 bp for prediction of genes using Glimmer MG.

#-------------------------------------------------------------------------------------------------
# Modules to use regular expressions and command line arguments.
import re
import sys
# Modules to deal with the file system.
import os
# Module to parse arguments
import argparse
# Biopython modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Filter contigs shorter than a threshold value.')
parser.add_argument("-f", "--fasta", help="Introduce the name of the fasta file containing the assemblies")
parser.add_argument("-t", "--threshold", type=int, default=0, help="Introduce the threshold value")
parser.add_argument("-o", "--output", type=str, default="contigs_filtered", help="Introduce the name of the output file")
args = parser.parse_args()
# Variable that stores fasta sequences
fasta_file = args.fasta
# Variable to store threshold value
threshold = args.threshold
# Variable to store the name of the output file
output = str(args.output)


#----------------------------------- CREATE A FILTERED FILE --------------------------------------
# Define a function that takes the fasta file and the threshold value as input and write the output

def filtercontigs(fasta_file, threshold):
   handle = open(str(fasta_file), "rU")
   contigs_filtered = dict()
   # Iterate through every contig in the fasta file
   for record in SeqIO.parse(handle, "fasta"):
      # Determine size of contigs to analyse
      if len(record.seq) >= threshold:
         contigs_filtered[record.id+"_length=_"+str(len(record.seq))] = str(record.seq)
   return contigs_filtered
   handle.close()

contigsfiltered = filtercontigs(fasta_file, threshold)

# Create a file to save the output with the name given by the user, open it and indicate it will 
# be a text file (wt).
file = open(output+"_"+str(threshold)+".fasta", 'wt')
# Write all the contigs of the filtered dict to a text file
for idname, seq in contigsfiltered.items():
   file.write(idname+"\n"+seq+"\n")
# close the file
file.close
