import sys,getopt
import fileinput

def process(line):
	strings = line.strip().split('\t')
	pos2 = int(strings[2])
	pos1 = int(strings[1])
	c = pos2 - pos1
	return c
	

result = 0

for line in fileinput.input(sys.argv[1]):
	C = process(line)
	result = result + C

print(result)