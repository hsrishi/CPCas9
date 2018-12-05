#! /usr/bin/python

"""
Created April 27, 2017
@author = Harneet Rishi

Use:
Takes in an input folder containing trimmed NGS files ready for mapping and index and yields:
	1. sai alignment file
	2. sam alignment file
	3. converts sam to bam
	4. sorts the bam file
	5. indexes the bam file
	6. converts sorted bam to sorted sam
for each trimmed read file

Example:
python bwa_samtools.py -i Output -d dCas9_CP
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
from subprocess import call
import os

###Functions for use in main
def get_seq_files(input_folder):
	seq_files = [f for f in listdir(input_folder) if isfile(join(input_folder,f))]
	seq_files = [f for f in seq_files if '.fastq' in f]
	return seq_files

if __name__ == "__main__":
	#if no arguments are passed
	if len(sys.argv) < 2:
		raise ValueError('Please provide all required inputs')
	else:
		#set up optional argument parsing
		myparser = optparse.OptionParser()
		myparser.add_option("-i", "--input", action = "store",
			help = "Input folder name containing fastq files to be mapped")
		myparser.add_option("-d", "--db_index", action = "store",
			help = "File containing bwa index for aligment")

	#parsing optional arguments
	(options, args) = myparser.parse_args()
	input_folder = options.input
	db_index = options.db_index

	#get seq files
	seq_files = get_seq_files(input_folder)

	#change directory to input folder
	os.chdir(os.getcwd()+'/'+input_folder)

	#align using BWA and create alignment files
	for seq_file in seq_files:
		base_name = seq_file.split('.')[0]
		sai_R1 = "mapped_%s_%s.sai" %(db_index, base_name)
		sam_R1 = "mapped_%s_%s.sam" %(db_index, base_name)
		bam_R1 = "mapped_%s_%s.bam" %(db_index, base_name)
		bam_R1_sort_name = "mapped_%s_%s_sort" %(db_index, base_name)
		bam_R1_sort_index = "%s.bam" %(bam_R1_sort_name)
		sam_R1_sort_index = "%s.sam" %(bam_R1_sort_name)

		print("***** Now analyzing: " + seq_file)
		#1. sai alignment file
		print("* aligning with bwa")
		call(['./bwa aln "%s" "%s" > "%s"' % (db_index, seq_file, sai_R1)], shell = True)

		#2. sam alignment file
		print("* generating sam file")
		call(['./bwa samse "%s" "%s" "%s" > "%s"' % (db_index, sai_R1, seq_file, sam_R1)], shell = True)

		#3. convert sam to bam
		print("* converting sam to bam")
		call(['./samtools view -ut "%s" "%s" > "%s"' % (db_index, sam_R1, bam_R1)], shell = True)

		#4. sort bam file
		print("* sorting bam file")
		call(['./samtools sort "%s" "%s"' % (bam_R1, bam_R1_sort_name)], shell = True)

		#5. index bam file
		print("* indexing bam file")
		call(['./samtools index "%s"' % (bam_R1_sort_index)], shell = True)

		#6. convert sorted bam to sorted sam
		print('* converting sorted bam to sorted sam')
		call(['./samtools view -h "%s" > "%s"' % (bam_R1_sort_index, sam_R1_sort_index)], shell = True)
		

