#!/usr/bin/env python3


import argparse
import os
import subprocess
import datetime
import re
import string
from StretcherReader import findglobalalign
from graphcigar import graphcigar



def getalignrregion(fullcigar, rsize):
	
	
	segments = re.findall(r'[<>][^<^>]+',fullcigar)
	
	cigars = re.findall(r'[<>]*\d+[MIDXH=]',fullcigar)
	
	rregions = []
	
	for segment in segments:
		
		rposi = 0
		strd = -1 if segment[0] == '<' else 1
		
		cigars = re.findall(r'\d+[MIDXH=]',segment)
		if strd  == -1:
			cigars = cigars[::-1]
			
		for cigar in cigars:
			
			thetype = cigar[-1]
			thesize = int(cigar[:-1])
			if thetype in ['M','H','D','X','=']:
				
				if thetype in ['M','='] or (thetype in ['X'] and thesize < 20):
					rregions.append([rposi,rposi + thesize])
					
				rposi += thesize 
				
			lasttype = thetype
			
	new_rregions = []
	lastregion = [-31,-31]
	for region in rregions:
		
		if abs(region[0] - lastregion[1])<30:
			new_rregions[-1][1] = region[1]
			
		else:
			new_rregions.append(region)
			
		lastregion = region
		
		
	return new_rregions


def getseqs(inputfile):
	
	with open(inputfile, mode ='r') as f:
		read = f.read()
		
	if read=='':
		return ["","",""]
	
	alllines=[a for a in read.split('\n')]
	#ref: original genome, query:assemblies
	alllinesindex=[a for a in range(len(alllines)) if len(alllines[a])>0 and alllines[a][0]!='#' and alllines[a][0]!=':' and  re.match(r'\S+\s\S+',alllines[a])!=None]
	
	eachline=alllinesindex[::2]
	qlines=[a for a in eachline]
	rlines=[a+2 for a in eachline]
	alines=[a+1 for a in eachline]
	if qlines==[]:
		return ["","",""]
	
	line0=alllines[qlines[0]]
	lineend = len(line0.strip())
	linestart = len(line0.split()[0]) + 1
	
	query=''.join([alllines[a][linestart:lineend].strip() for a in qlines])
	ref=''.join([alllines[a][linestart:lineend].strip() for a in rlines])
	align=''.join([(alllines[a])[linestart:lineend] for a in alines])[:len(ref)]
	
	query = query.replace("-","")
	ref = ref.replace("-","")
	align = align.replace("-","")
	
	return [query, ref, align]


def aligngraph(graphfile, queryfile,nthreads=4):
	
	tempfolder = graphfile + "_tempfolder/"
	
	graphfilename = graphfile.split("/")[-1]
	
	dbfile = "{}_db".format(graphfile)
	
	alignout = "{}_align.out".format(queryfile)
	
	dbcmd = "makeblastdb -in {}  -dbtype nucl -parse_seqids -out {}".format(graphfile,  dbfile)
	
	os.system(dbcmd)
	
	cmd = "blastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30 -perc_identity 90 -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( queryfile,  dbfile, alignout ,nthreads)
	
	os.system(cmd)
	
	queryname, fullpath, fullcigar, pathranges,qranges = graphcigar(graphfile, queryfile , alignout)
	
	return queryname, fullpath, fullcigar, pathranges,qranges





def main(args):
	
	score, cigartext = findglobalalign(args.input)
	
	query, ref, align = getseqs(args.input)
	
	allcigars = re.findall(r'[<>]*\d+[IDXM][A-Za-z]*',cigartext)
	
	tempfolder = args.input + "_tempfolder/"
	
	try:
		os.mkdir(tempfolder)
	except:
		pass
		
	queryfile = "{}/query.fa".format(tempfolder)
	with open(queryfile, mode = 'w') as f:
		
		f.write(">query\n{}".format(query))
		
	reffile = "{}/ref.fa".format(tempfolder)
	with open(reffile, mode = 'w') as f:
		
		f.write(">ref\n{}".format(ref))
		
	allinsertfile = "{}/allinserts.fa".format(tempfolder)
	alldeletefile = "{}/alldeletes.fa".format(tempfolder)

	rposi = 0
	qposi = 0
	with open(allinsertfile, mode = 'w') as allinsertfile_: 
		with open(alldeletefile, mode = 'w') as alldeletefile_: 
			for i,cigar in enumerate(allcigars):
			
				size = int(re.findall(r'\d+',cigar)[0])
				if "I" in cigar:
				
					if size > 200:
									
						cigarseq = re.findall(r'[a-zA-Z]+',cigar)[0][1:]
						
						allinsertfile_.write(">Insertion_{}\n{}\n".format(i,cigarseq))
					qposi += size
			
				elif "D" in cigar:
					
					if size > 200:
					
						cigarseq = re.findall(r'[a-zA-Z]+',cigar)[0][1:]
						
						alldeletefile_.write(">Deletion_{}\n{}\n".format(i,ref[rposi:(rposi+size)]))
					rposi += size
									
				else:
					qposi += size
					rposi += size	
	
	SVs = []
	rposi = 0
	qposi = 0
	for i,cigar in enumerate(allcigars):
		
		size = int(re.findall(r'\d+',cigar)[0])
		if "I" in cigar:
			
			if size > 200:
				
				seqfile = tempfolder+"Insert_{}.fa".format(i)
				
				cigarseq = re.findall(r'[a-zA-Z]+',cigar)[0][1:]
				
				with open(seqfile, mode = 'w') as f: 
					
					f.write(">Insertion_{}\n{}".format(i,cigarseq))
				
					
				queryname, fullpath, fullcigar, pathranges,qranges = aligngraph(reffile, seqfile)
				
				alignregions = getalignrregion(fullcigar, len(ref))
				rstarts = ";".join([str(x[0]) for x in alignregions])
				rends = ";".join([str(x[1]) for x in alignregions])
				
				if  len(rstarts) == 0:
					rstarts = -1
					rends = -1
					
				matchedsize = sum([int(x[:-1]) for x in re.findall(r'\d+[=M]', fullcigar)] + [0])
				
				typename = "Insertion"
				if matchedsize > 0.5 * size:
					typename = "Duplication" 
				
				else:
					queryname2, fullpath2, fullcigar2, pathranges2,qranges2 = aligngraph(alldeletefile,seqfile)
				
					matchedsize2 = sum([int(x[:-1]) for x in re.findall(r'\d+[=M]', fullcigar2)] + [0])

					if matchedsize2 > 0.5 * size:
						typename = "TransInsertion"
				
					
				SVs.append([i,typename, qposi, qposi+size, rstarts, rends , size, matchedsize, fullcigar])
				
			qposi += size
			
		elif "D" in cigar:
			if size > 200:
				
				seqfile = tempfolder+"Deletion_{}.fa".format(i)
				
				cigarseq = re.findall(r'[a-zA-Z]+',cigar)[0][1:]
				
				with open(seqfile, mode = 'w') as f: 
					
					f.write(">Deletion_{}\n{}".format(i,ref[rposi:(rposi+size)]))
					
				queryname, fullpath, fullcigar, pathranges,qranges = aligngraph(queryfile,seqfile)
				
				alignregions = getalignrregion(fullcigar, len(ref))
				qstarts = ";".join([str(x[0]) for x in alignregions])
				qends = ";".join([str(x[1]) for x in alignregions])
				
				if  len(qstarts) == 0:
					qstarts = -1
					qends = -1
					
				matchedsize = sum([int(x[:-1]) for x in re.findall(r'\d+[=M]', fullcigar)] + [0])
				
				typename = "Deletion"
				if matchedsize > 0.5 * size:
					typename = "Contraction"
					
				else:
					queryname2, fullpath2, fullcigar2, pathranges2,qranges2 = aligngraph(allinsertfile,seqfile)
					
					matchedsize2 = sum([int(x[:-1]) for x in re.findall(r'\d+[=M]', fullcigar2)] + [0])
					
					if matchedsize2 > 0.5 * size:
						typename = "TransDeletion"
				
				fullcigar = fullcigar.translate(str.maketrans('DITCGNtcgn','IDAAAAAAAA')).replace('A','') 		
					
				SVs.append([i, typename, qstarts, qends,rposi, rposi + size, size, matchedsize,fullcigar])
				
			rposi += size 
		else:
			qposi += size
			rposi += size
		
	text_out = "\n".join(["\t".join([str(x) for x in sv]) for sv in SVs]) + '\n'
	
	
	with open(args.output, mode = 'w') as f:
		f.write( text_out)
		
	os.system("rm -rf "+tempfolder)
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	
	parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
	parser.add_argument("-p", "--edgepolish", help="path to output file", dest="ifpolish", type=int,default = 1)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
