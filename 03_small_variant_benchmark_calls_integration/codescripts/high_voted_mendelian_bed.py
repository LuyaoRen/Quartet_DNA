import pandas as pd
import sys, argparse, os
mut = pd.read_table(sys.argv[1])
outFile = open(sys.argv[2],'w')
for row in mut.itertuples():
#d5
	if ',' in row.V4:
		alt = row.V4.split(',')
		alt_len = [len(i) for i in alt]
		alt_max = max(alt_len)
	else:
		alt_max = len(row.V4)
#d6
	alt = alt_max
	ref = row.V3
	pos = int(row.V2)
	if len(ref) == 1 and alt == 1:
		StartPos = int(pos) -1
		EndPos = int(pos)
		cate = 'SNV'
	elif len(ref) > alt:
		StartPos = int(pos) - 1
		EndPos = int(pos) + (len(ref) - 1)
		cate = 'INDEL'
	elif alt > len(ref):
		StartPos = int(pos) - 1
		EndPos = int(pos) + (alt - 1)
		cate = 'INDEL'
	elif len(ref) == alt:
		StartPos = int(pos) - 1
		EndPos = int(pos) + (alt - 1)
		cate = 'INDEL'
	outline = row.V1 + '\t' + str(StartPos) + '\t' + str(EndPos) + '\t' + str(row.V2) + '\t' + cate + '\n'
	outFile.write(outline)








