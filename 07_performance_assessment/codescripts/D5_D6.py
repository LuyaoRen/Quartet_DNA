from __future__ import division
import pandas as pd
import sys, argparse, os


# input arguments
parser = argparse.ArgumentParser(description="this script is to calculate reproducibility between Quartet_D5 and Quartet_D6s")

parser.add_argument('-sister', '--sister', type=str, help='sister.txt',  required=True)
parser.add_argument('-project', '--project', type=str, help='project name',  required=True)


args = parser.parse_args()
sister_file = args.sister
project_name = args.project

# output file
output_name = project_name + '.sister.reproducibility.txt'

output_file = open(output_name,'w')

# input files
sister_dat = pd.read_table(sister_file)

indel_sister_same = 0
indel_sister_diff = 0
snv_sister_same = 0
snv_sister_diff = 0

for row in sister_dat.itertuples():
	# snv indel
	if ',' in row[4]:
		alt = row[4].split(',')
		alt_len = [len(i) for i in alt]
		alt_max = max(alt_len)
	else:
		alt_max = len(row[4])
	alt = alt_max
	ref = row[3]
	if len(ref) == 1 and alt == 1:
		cate = 'SNV'
	elif len(ref) > alt:
		cate = 'INDEL'
	elif alt > len(ref):
		cate = 'INDEL'
	elif len(ref) == alt:
		cate = 'INDEL'
	# sister
	if row[5] == row[6]:
		if row[5] == './.':
			mendelian = 'noInfo'
			sister_count = "no"
		elif row[5] == '0/0':
			mendelian = 'Ref'
			sister_count = "no"
		else:
			mendelian = '1'
			sister_count = "yes_same"
	else:
		mendelian = '0'
		if (row[5] == './.' or row[5] == '0/0') and (row[6] == './.' or row[6] == '0/0'):
			sister_count = "no"
		else:
			sister_count = "yes_diff"
	if cate == 'SNV':
		if sister_count == 'yes_same':
			snv_sister_same += 1
		elif sister_count == 'yes_diff':
			snv_sister_diff += 1
		else:
			pass
	elif cate == 'INDEL':
		if sister_count == 'yes_same':
			indel_sister_same += 1
		elif sister_count == 'yes_diff':
			indel_sister_diff += 1
		else:
			pass	

indel_sister = indel_sister_same/(indel_sister_same + indel_sister_diff)
snv_sister = snv_sister_same/(snv_sister_same + snv_sister_diff)
outcolumn =  'Project\tReproducibility_D5_D6\n'
indel_outResult = project_name + '.INDEL' + '\t' + str(indel_sister)  + '\n'
snv_outResult = project_name + '.SNV' + '\t' + str(snv_sister)  + '\n'
output_file.write(outcolumn)
output_file.write(indel_outResult)
output_file.write(snv_outResult)



