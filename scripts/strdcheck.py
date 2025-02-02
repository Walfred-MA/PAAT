#!/usr/bin/env python3

import os 
import re
import argparse
import sys


def checkstrand(inputfile):
	
	with open(inputfile, mode= 'r') as f:
		align = f.read()
	f.close()
	
	score0=0.0
	try:
		score0=float(re.search(r'(?<=# Score: )\s*-*[0-9]+\.*[0-9]*', align).group().strip())
	except:
		score0=0.0

	score1 = 0	
	if score0 < 0:
		aseq = re.search(r'(?<=\[-asequence\]\s).+', align).group().strip()
		bseq = re.search(r'(?<=\[-bsequence\]\s).+', align).group().strip()
	
		globalfile2 = inputfile + "_reverse.txt"
	
		cm='stretcher {} {} -snucleotide2  -gapopen 16  -gapextend 4 -sreverse1  {}'.format(aseq,bseq,globalfile2)
		os.system(cm)
		
		with open(globalfile2, mode= 'r') as f:
			align = f.read()
		f.close()
		
		score1=0.0
		try:
			score1=float(re.search(r'(?<=# Score: )\s*-*[0-9]+\.*[0-9]*', align).group().strip())
		except:
			score1=0.0
		
		if score1 > score0:
			os.system("mv {} {}".format(inputfile, inputfile+"_reverse"))
			os.system("mv {} {}".format(globalfile2, inputfile))
		else:
			os.system("rm {} ".format(globalfile2))
		
	return int(score1 > score0)


def strand_check(args):
		
	checkstrand(args.input)
			
			
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.set_defaults(func=strand_check)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()

