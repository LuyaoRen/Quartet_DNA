import pandas as pd
import sys, argparse, os
mut = mut = pd.read_table('/mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/vcf/mutation_type',header=None)
vote = pd.read_table('/mnt/pgx_src_data_pool_4/home/renluyao/manuscript/benchmark_calls/all_info/benchmark.vote.mendelian.txt',header=None)
merged_df = pd.merge(vote, mut,  how='inner', left_on=[0,1], right_on = [0,1])
outFile = open(sys.argv[1],'w')
outIndel = open(sys.argv[2],'w')
for row in merged_df.itertuples():
#d5
	if ',' in row._7:
		d5 = row._7.split(',')
		d5_len = [len(i) for i in d5]
		d5_alt = max(d5_len)
	else:
		d5_alt = len(row._7)
#d6
	if ',' in row._15:
		d6 = row._15.split(',')
		d6_len = [len(i) for i in d6]
		d6_alt = max(d6_len)
	else:
		d6_alt = len(row._15)
#f7
	if ',' in row._23:
		f7 = row._23.split(',')
		f7_len = [len(i) for i in f7]
		f7_alt = max(f7_len)
	else:
		f7_alt = len(row._23)
#m8
	if ',' in row._31:
		m8 = row._31.split(',')
		m8_len = [len(i) for i in m8]
		m8_alt = max(m8_len)
	else:
		m8_alt = len(row._31)
	all_length = [d5_alt,d6_alt,f7_alt,m8_alt]
	alt = max(all_length)
	ref = row._35
	pos = int(row._2)
	if len(ref) == 1 and alt == 1:
		StartPos = int(pos) -1
		EndPos = int(pos)
		cate = 'SNV'
	elif len(ref) > alt:
		StartPos = int(pos) - 1
		EndPos = int(pos) + (len(ref) - 1)
		cate = 'INDEL'
		outline_indel = row._1 + '\t' + str(StartPos) + '\t' + str(EndPos) + '\n'
		outIndel.write(outline_indel)
	elif alt > len(ref):
		StartPos = int(pos) - 1
		EndPos = int(pos) + (alt - 1)
		cate = 'INDEL'
		outline_indel = row._1 + '\t' + str(StartPos) + '\t' + str(EndPos) + '\n'
		outIndel.write(outline_indel)
	elif len(ref) == alt:
		StartPos = int(pos) - 1
		EndPos = int(pos) + (alt - 1)
		cate = 'INDEL'
		outline_indel = row._1 + '\t' + str(StartPos) + '\t' + str(EndPos) + '\n'
		outIndel.write(outline_indel)
	outline = row._1 + '\t' + str(StartPos) + '\t' + str(EndPos) + '\t' + str(row._2) + '\t' + cate + '\n'
	outFile.write(outline)








