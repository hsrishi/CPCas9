#! /usr/bin/python

"""
Created March 26, 2017
@author = Harneet Rishi

Use:
Takes in an input folder containing untrimmed NGS files (fastq) and a csv file (parameters) containing 
the length of the variable region (column 1) and the constant primer binding region (column 2)
and returns an output folder containing trimmed reads for BWA. Note that the trimmed reads will
have the constant region removed

Example:
python trim_reads_forBWA.py -i NGS_files -p parameters_withmore.csv -o output -r 14

"""


###Import relevant modules
from __future__ import division
from Bio import SeqIO
import re
import pandas as pd
import numpy as np
import subprocess
import sys
import optparse
from os import listdir, makedirs
from os.path import isfile, join
import os

###Functions for use in main
def get_seq_files(input_folder):
	seq_files = [f for f in listdir(input_folder) if isfile(join(input_folder,f))]
	seq_files = [f for f in seq_files if '.fastq' in f]
	return seq_files

def trim_file(input_folder, seq_file, length_var_region, primer_binding_seq, length_ROI):
	start_pos = length_var_region
	end_pos = length_var_region + len(primer_binding_seq)
	trimmed_reads = (rec[end_pos+3:end_pos+3+length_ROI] for rec in \
		SeqIO.parse('%s/%s' %(input_folder, seq_file), "fastq") \
		if ( (str(rec.seq[start_pos:end_pos])) == primer_binding_seq) and np.mean(rec.letter_annotations["phred_quality"] >= 20)
		)
	return trimmed_reads

if __name__ == "__main__":
	#if no arguments are passed
	if len(sys.argv) < 2:
		raise ValueError('Please provide all required inputs')
	else:
		#set up optional argument parsing
		myparser = optparse.OptionParser()
		myparser.add_option("-i", "--input", action = "store",
			help = "Input folder name containing fastq files to be trimmed")
		myparser.add_option("-p", "--parameters", action = "store",
			help = "Parameters file (csv) containing length of var region and primer binding seq")
		myparser.add_option("-r", "--ROI", action = "store",
			help = "Length of region of interest to keep after constant region")
		myparser.add_option("-o", "--output", action = "store",
			help = "Output folder name to add trimmed files")

	#parsing optional arguments
	(options, args) = myparser.parse_args()
	input_folder = options.input
	parameters_file = options.parameters
	output_folder = options.output
	length_ROI = int(options.ROI)

	if not(os.path.exists(output_folder)):
		os.makedirs(output_folder)

	#get seq files
	seq_files = get_seq_files(input_folder)

	#parse parameters file
	df_parameters = pd.read_csv(parameters_file)
	length_var_region = int(df_parameters['length_var_region'][0])
	primer_binding_seq = df_parameters['primer_binding_seq'][0]

	#trim files
	i = 0
	for seq_file in seq_files:
		print("***** Now analyzing: " + seq_file)
		num_lines = subprocess.check_output(['wc -l %s/%s' %(input_folder, seq_file)], shell = True)
		num_lines_split = num_lines.split(' ')
		num_reads = int(num_lines_split[1])/4
		#num_reads = int(num_lines_split[3])/4 #changed to run on server
		print("%i total reads" %num_reads)

		trimmed_reads = trim_file(input_folder, seq_file, length_var_region, primer_binding_seq, length_ROI)

		count = SeqIO.write(trimmed_reads, "%s/trimmed_forBWA_%smer_%s.fastq" % (output_folder, length_ROI, seq_file.split('.')[0]), "fastq")
		print("Saved %i reads" % count)

		i = i + 1
		print("*Done trimming: "+seq_file)




