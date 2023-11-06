from __future__ import division
import pandas as pd
import sys, argparse, os
import fileinput
import re

# input arguments
parser = argparse.ArgumentParser(description="this script is to extract mendelian concordance information")

parser.add_argument('-LCL5', '--LCL5', type=str, help='LCL5 family info',  required=True)
parser.add_argument('-LCL6', '--LCL6', type=str, help='LCL6 family info',  required=True)
parser.add_argument('-family', '--family', type=str, help='family name',  required=True)


args = parser.parse_args()
lcl5 = args.LCL5
lcl6 = args.LCL6
family = args.family


# output file
family_name = family + '.txt'

family_file = open(family_name,'w')

# input files
lcl5_dat = pd.read_table(lcl5)
lcl6_dat = pd.read_table(lcl6)

merged_df = pd.merge(lcl5_dat, lcl6_dat,  how='outer', left_on=['#CHROM','POS'], right_on = ['#CHROM','POS'])

def alt_seq(alt, genotype):
	if genotype == './.':
		seq = './.'
	elif genotype == '0/0':
		seq = '0/0'
	else:
		alt = alt.split(',')
		genotype = genotype.split('/')
		if genotype[0] == '0':
			allele2 = alt[int(genotype[1]) - 1]
			seq = '0/' + allele2
		else:
			allele1 = alt[int(genotype[0]) - 1]
			allele2 = alt[int(genotype[1]) - 1]
			seq = allele1 + '/' +  allele2
	return seq

for row in merged_df.itertuples():
	# correction of multiallele
	if pd.isnull(row.INFO_x) == True or pd.isnull(row.INFO_y) == True:
		mendelian = '.'
	else:
		lcl5_seq = alt_seq(row.ALT_x, row.CHILD_x)
		lcl6_seq = alt_seq(row.ALT_y, row.CHILD_y)
		if lcl5_seq == lcl6_seq:
			mendelian = '1'
		else:
			mendelian = '0'
	if pd.isnull(row.INFO_x) == True:
		mendelian = mendelian + ':.'
	else:
		mendelian = mendelian + ':' + row.INFO_x.split('=')[1]
	if pd.isnull(row.INFO_y) == True:
		mendelian = mendelian + ':.'
	else:
		mendelian = mendelian + ':' + row.INFO_y.split('=')[1]


	outline = row._1 + '\t' + str(row.POS) + '\t' + mendelian + '\n'
	family_file.write(outline)
