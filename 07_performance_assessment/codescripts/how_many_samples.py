import pandas as pd
import sys, argparse, os
from operator import itemgetter

parser = argparse.ArgumentParser(description="This script is to get how many samples")

parser.add_argument('-sample', '--sample', type=str, help='quartet_sample',  required=True)
parser.add_argument('-rep', '--rep', type=str, help='quartet_rep',  required=True)
args = parser.parse_args()

# Rename input:
sample = args.sample
rep = args.rep

quartet_sample = pd.read_table(sample,header=None)
quartet_sample = list(quartet_sample[0])
quartet_rep = pd.read_table(rep.header=None)
quartet_rep = quartet_rep[0]

#tags
sister_tag = 'false'
quartet_tag = 'false'

quartet_rep_unique = list(set(quartet_rep))

single_rep = [i for i in range(len(quartet_rep)) if quartet_rep[i] == quartet_rep_unique[0]]
single_batch_sample = itemgetter(*single_rep)(quartet_sample)

num = len(single_batch_sample)
if num == 1:
	sister_tag = 'false'
	quartet_tag = 'false'
elif num == 2:
	if set(single_batch_sample) == set(['LCL5','LCL6']):
		sister_tag = 'true'
		quartet_tag = 'false'
elif num == 3:
	if ('LCL5' in single_batch_sample) and ('LCL6' in single_batch_sample):
		sister_tag = 'true'
		quartet_tag = 'false'
elif num == 4:
	if set(single_batch_sample) == set(['LCL5','LCL6','LCL7','LCL8']):
		sister_tag = 'false'
		quartet_tag = 'true'

sister_outfile = open('sister_tag','w')
quartet_outfile = open('quartet_tag','w')

sister_outfile.write(sister_tag)
quartet_outfile.write(quartet_tag)
