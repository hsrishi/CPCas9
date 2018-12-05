#! /usr/bin/python

"""
Created August 9, 2017
@author = Harneet Rishi

Use: Takes in a sorted same file (e.g. output from bwa_samtools.py; filename hardcoded in line 33) and
provides an nucleotide abundance mapping of functional CP variants (those with in-frame insertions mapping
to the forward strand).



Example:
python process_SAM.py


"""

from __future__ import division

from Bio import SeqIO
from Bio.Seq import Seq
import re

import numpy as np
import pandas as pd
from pandas.io.parsers import *
import os

################################################################
#Import
################################################################
print '*** Importing'
#samples = [f for f in os.listdir('Output') if f.endswith('sort.sam')]
samples = ['mapped_dCas9_CP_trimmed_forBWA_14mer_DFSBLO002D_S4_L002_R1_001_sort.sam'] #hardcode sample b/c doing all at same time is too memory intensive
sample_root = ['%s_%s' %(sample.split('_')[6], sample.split('_')[8]) for sample in samples]
dict_sample_names = {'%s_%s' %(sample.split('_')[6], sample.split('_')[8]):sample for sample in samples}

folder_sam = 'Output'
names = ['QNAME', 'FLAG', 'RNAME', 'POS','MAPQ', 'CIGAR','MRNM','MPOS','ISIZE','SEQ','QUAL',
        'OPT:XT','OPT:NM','OPT:X0','OPT:X1','OPT:XM','OPT:XO','OPT:XG','OPT:MD']

dict_df_sam = {}

for sample in samples:
    name = '%s_%s' %(sample.split('_')[6], sample.split('_')[8])
    dict_df_sam[name] = pd.read_table('%s/%s' %(folder_sam, dict_sample_names[name]),
                       names = names, index_col = False, skiprows = 3)


################################################################
#Filter
################################################################
print '*** Filtering'

dict_map_CP_filter = {}
dict_map_CP_for = {}
dict_map_CP_rev = {}

for key in dict_df_sam.keys():
    df = dict_df_sam[key].copy()
    df = df[df['RNAME'] == 'dCas9_full_CPlinker_NoStop']
    dict_map_CP_filter[key] = df[(df['MAPQ'] == 37) & (df['CIGAR'] == '14M') & (df['OPT:XO'] == 'XO:i:0') & (df['OPT:XG'] == 'XG:i:0')]
    dict_map_CP_for[key] = dict_map_CP_filter[key][dict_map_CP_filter[key]['FLAG'] == 0]
    dict_map_CP_rev[key] = dict_map_CP_filter[key][dict_map_CP_filter[key]['FLAG'] == 16]
    dict_map_CP_for[key]['Out of Frame'] = dict_map_CP_for[key]['POS'].apply(lambda x: (x-1) % 3)

################################################################
#Export dataframes
################################################################
print '*** Exporting Dataframes'

output_folder = 'Output'

for key in dict_map_CP_filter.keys():
	dict_map_CP_filter[key].to_csv('%s/df_map_CP_filter_%s.csv' %(output_folder, key))

for key in dict_map_CP_for.keys():
	dict_map_CP_for[key].to_csv('%s/df_map_CP_for_%s.csv' %(output_folder, key))

for key in dict_map_CP_rev.keys():
	dict_map_CP_rev[key].to_csv('%s/df_map_CP_rev_%s.csv' %(output_folder, key))

################################################################
#Alignment statistics
################################################################
print '*** Calculating Alignment Statistics'

output_folder = 'Output'
for key in dict_map_CP_filter.keys():
	cols = dict_map_CP_filter[key].columns
	cols = [col for col in cols if not((col == 'QNAME') | (col == 'SEQ') | (col == 'QUAL') | (col == 'POS'))]
	Stats_Alignment_file = open('%s/Alignment_Stats_%s.txt'%(output_folder, key), 'w')
	for col in cols:
		Stats_Alignment_file.write('**********'+'\n')
		Stats_Alignment_file.write(col+'\n')
		Stats_Alignment_file.write('   ' + str(set(dict_map_CP_filter[key][col]))+'\n')
		Stats_Alignment_file.write('   ' + str(dict_map_CP_filter[key][col].value_counts())+'\n')
	Stats_Alignment_file.close()

################################################################
#Mapping statistics
################################################################
print '*** Calculating Mapping Statistics'

output_folder = 'Output'
for key in dict_map_CP_filter.keys():
	df = dict_map_CP_filter[key]
	Mapping_Alignment_file = open('%s/Mapping_Stats_%s.txt'%(output_folder, key), 'w')
	Mapping_Alignment_file.write(key+'\n')
	Mapping_Alignment_file.write('Reads filtered have been mapped to dCas9, CIGAR 14M, MAPQ 37, OPT:XO 0 gap open, OPT:XG 0 gap extension'+'\n')
	Mapping_Alignment_file.write('Filtered reads mapping forward, reverse'+'\n')
	Mapping_Alignment_file.write('    total reads after filtering: %d' %(len(df))+'\n')
	Mapping_Alignment_file.write('    forward mapped: %.2f percent' %(len(df[df['FLAG'] == 0])/len(df)*100)+'\n')
	Mapping_Alignment_file.write('    reverse mapped: %.2f percent' %(len(df[df['FLAG'] == 16])/len(df)*100)+'\n')
	Mapping_Alignment_file.write('In-frame'+'\n')
	df = dict_map_CP_for[key]
	Mapping_Alignment_file.write('    %.2f percent of forward reads are in-frame' %(len(df[df['Out of Frame'] == 0])/len(df)*100)+'\n')
	Mapping_Alignment_file.close()

################################################################
#Aggregate data on nucleotide position basis
################################################################
print '*** Aggregating Mapping and Exporting'

output_folder = 'Output'

#Keeping only in-frame positions
df_for_inFrame = {}
df_for_inFrame_gpby = {}
for key in dict_map_CP_for.keys():
    df = dict_map_CP_for[key].copy()
    df_for_inFrame[key] = df[df['Out of Frame'] == 0]
    df_for_inFrame_gpby[key] = df_for_inFrame[key].groupby('POS').count()
    df_for_inFrame_gpby[key]['Pos'] = df_for_inFrame_gpby[key].index
    df_for_inFrame_gpby[key] = df_for_inFrame_gpby[key]['POS']
    df_for_inFrame[key].to_csv('%s/df_for_inFrame_%s.csv' %(output_folder, key))
    df_for_inFrame_gpby[key].to_csv('%s/df_for_inFrame_gpby_%s.csv' %(output_folder, key))















