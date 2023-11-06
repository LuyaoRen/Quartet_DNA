import pandas as pd
import sys, argparse, os
mut = pd.read_table('/mnt/pgx_src_data_pool_4/home/renluyao/manuscript/MIE/vcf/mutation_type',header=None)
outIndel = open(sys.argv[1],'w')
for row in mut.itertuples():
	if ',' in row._4:
		alt_seq = row._4.split(',')
		alt_len = [len(i) for i in alt_seq]
		alt = max(alt_len)
	else:
		alt = len(row._4)
	ref = row._3
	pos = row._2
	if len(ref) == 1 and alt == 1:
		pass
	elif len(ref) > alt:
		StartPos = int(pos) - 1
		EndPos = int(pos) + (len(ref) - 1)
		outline_indel = row._1 + '\t' + str(StartPos) + '\t' + str(EndPos) + '\n'
		outIndel.write(outline_indel)
	elif alt > len(ref):
		StartPos = int(pos) - 1
		EndPos = int(pos) + (alt - 1)
		outline_indel = row._1 + '\t' + str(StartPos) + '\t' + str(EndPos) + '\n'
		outIndel.write(outline_indel)
	elif len(ref) == alt:
		StartPos = int(pos) - 1
		EndPos = int(pos) + (alt - 1)
		outline_indel = row._1 + '\t' + str(StartPos) + '\t' + str(EndPos) + '\n'
		outIndel.write(outline_indel)







