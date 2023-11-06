import sys,getopt
from itertools import islice

over_50_outfile = open("indel_lenth_over_50.txt",'w')
less_50_outfile = open("benchmark.men.vote.diffbed.lengthlessthan50.txt","w")

def process(line):
	strings = line.strip().split('\t')
#d5
	if ',' in strings[6]:
		d5 = strings[6].split(',')
		d5_len = [len(i) for i in d5]
		d5_alt = max(d5_len)
	else:
		d5_alt = len(strings[6])
#d6
	if ',' in strings[14]:
		d6 = strings[14].split(',')
		d6_len = [len(i) for i in d6]
		d6_alt = max(d6_len)
	else:
		d6_alt = len(strings[14])
#f7
	if ',' in strings[22]:
		f7 = strings[22].split(',')
		f7_len = [len(i) for i in f7]
		f7_alt = max(f7_len)
	else:
		f7_alt = len(strings[22])
#m8
	if ',' in strings[30]:
		m8 = strings[30].split(',')
		m8_len = [len(i) for i in m8]
		m8_alt = max(m8_len)
	else:
		m8_alt = len(strings[30])
#ref
	ref_len = len(strings[34])
	if (d5_alt > 50) or (d6_alt > 50) or (f7_alt > 50) or (m8_alt > 50) or (ref_len > 50):
		over_50_outfile.write(line)
	else:
		less_50_outfile.write(line)


input_file = open(sys.argv[1])  
for line in islice(input_file, 1, None):  
	process(line)