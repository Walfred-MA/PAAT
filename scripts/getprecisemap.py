#!/usr/bin/env python3

import os
import argparse
import re
import bisect
	

def cigar_findrrange(cigars, qstart , qend, outstrd = '+'):
	
	chunks = cigars.replace("<","><").split(">")[1:]
	refsizes = [sum([int(x[:-1]) for x in  re.findall(r'\d+[MX=DH]', chunk)]) for chunk in chunks]
	
	cigars = re.findall(r'[><]*\d+[a-zA-Z=][actgACTG]*', cigars)
	
	refposi = 0
	posi = 0
	
	chunkindex = -1
	newcigars = []
	
	refsize = 0
	lastrefsize = 0
	
	strd = 1
	laststrd  = 1
	
	breaks = []
	
	refstart = -1
	refend = -1

	matchsize = 0 
	
	for i,cigar in enumerate(cigars):
		
		lastrposi = refposi
		lastposi = posi
		
		thesize = re.findall(r'\d+', cigar)[0]
		thetype = cigar[len(thesize):len(thesize)+1]
		theseq = cigar[(len(thesize)+1):]
		
		thesize = int(thesize)
		
		if cigar[0] in ["<",'>'] :
			
			
			chunkindex += 1
			lastrefsize  = refsize 
			refsize = refsizes[chunkindex]
			
			newcigars.append([])
			matchsize = 0
			
			laststrd = strd
			if cigar[0]=='>':
				strd = 1
			else:
				strd = -1
				
			refposi = thesize
			
		elif thetype == "=" or thetype == "M":
			
			posi += thesize
			refposi += thesize
			
		elif thetype == "X" or thetype == "S":
			
			posi += thesize
			refposi += thesize
			
		elif thetype == "D" :
			
			refposi += thesize
			
		elif thetype == "I" :
			
			posi += thesize
			
			
		if posi <= qstart:
			
			continue
	
	
		if cigar[0] in ["<",'>'] :
			
			if laststrd == -1:
				breaks.append(refsize -  lastrposi)
			else:
				breaks.append(lastrposi)
				
			if strd == -1:
				breaks.append(refsize - refposi)
			else:
				breaks.append(refposi)
				
		if lastposi <= qstart:
			
			if thetype in ['M','=', 'X','S', 'I','H']:
				
				refstart = refposi - thesize + qstart - lastposi
			else:
				refstart = refposi

				
			if strd == -1:
				breaks.append(refsize - refstart)
			else:
				breaks.append( refstart)
				
			corrsize = min(qend,posi) - qstart

			newcigars[-1].append("{}{}H".format('>' if strd > 0 else '<', refstart))
			newcigars[-1].append("{}{}".format( corrsize,thetype))
			
			
			
		if posi >= qend:
			
			left = 0
			if thetype in ['M','=', 'X','S', 'I','H']:
				
				refend = refposi - posi + qend
			else:
				refend = refposi
				
				
			if strd == -1:
				breaks.append(refsize - refend)
			else:
				breaks.append( refend)
				
				
			if lastposi > qstart:
				corrsize = qend - max(qstart, lastposi)
				
				newcigars[-1].append("{}{}".format(corrsize,thetype))
			
			newcigars[-1].append("{}H".format(refsize - refend))
			break
		
		if lastposi > qstart:
			
			newcigars[-1].append(cigar)
			
	if posi <= qstart:
		
		posi = qstart
		breaks.append(refposi)
		
	if posi < qend:
		
		newcigars[-1].append("{}I".format(qend-posi))
		
		breaks.append(refposi)
	
	newcigars = [x for x in newcigars if x != []]
	
	cigars_exclude = []
	leftover = []
	exclude = []
	for i, chunk in enumerate(newcigars):

		if chunk[0][-1] != 'H':
			chunk = [chunk[0][0] + "0H"] +chunk[1:]
		matchsize = sum([int(x[:-1]) for x in chunk if x[-1] in ['M','=']])
		if matchsize < 200:
			exclude.append(i)
			insert = sum([int(x[:-1]) for x in chunk if x[-1] in ['M','=','I']])
			insert_cigar = ""			
			if insert:
				insert_cigar = "{}I".format(insert)

			if len(cigars_exclude):
				cigars_exclude[-1].insert(-1,insert_cigar)

			else:
				leftover.append(insert_cigar)

		else:
			if len(leftover):
				chunk = [chunk[0]] + leftover + chunk[1:]
				leftover.clear()

			cigars_exclude.append(chunk)
		
	
	if len(leftover) == 0:
		newcigars = cigars_exclude
		breaks = [x for i,x in enumerate(breaks) if i//2 not in exclude]


	
	if outstrd  == '+':
		fullcigar = "".join(["".join(x) for x in newcigars])
	else:

		fullcigar = "".join(["".join(x) for x in newcigars[::-1]][::-1]).translate(str.maketrans('<>', '><'))
	
	
	return  breaks, fullcigar

def polishstartend(cigars):

	cigars = re.findall(r'\d+[MXID]',cigars)

	lindex = 0	
	score = 0
	for lindex, cigar in enumerate(cigars):
		
		size = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == 'M':
			
			score += size
			
			if score > 200 or size > 200:
				
				break
			
		elif thetype in ['X','D']:
			
			score -= 4 * size
			
	ltrim = cigars[:lindex]
	lrstart = sum([int(x[:-1]) for x in cigars[:lindex] if x[-1] in ['D','M','X']])
	lqstart = sum([int(x[:-1]) for x in cigars[:lindex] if x[-1] in ['I','M','X']])
	cigars = (cigars[lindex:])[::-1]
	
	rindex = 0
	score = 0
	for rindex, cigar in enumerate(cigars):
		
		size = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == 'M':
			
			score += size
			
			if score > 200 or size > 200:
				
				break
			
		elif thetype in ['X','D']:
			
			score -= 4 * size
			
	rtrim = cigars[:rindex]
	rrstart = sum([int(x[:-1]) for x in cigars[:rindex] if x[-1] in ['D','M','X']])
	rqstart = sum([int(x[:-1]) for x in cigars[:rindex] if x[-1] in ['I','M','X']])	
	cigars = (cigars[lindex:])[::-1]
	
	cigars = (["{}I".format(lqstart)] if lqstart > 0 else [])+cigars+(["{}I".format(rqstart)] if rqstart >0 else [])
	
	newcigars = []
	lasttype = ""
	lastsize = 0
	for cigar in cigars:
		size = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == lasttype:
			newcigars[-1] = "{}{}".format(size+lastsize, thetype)
			lasttype = thetype
			lastsize += size
		else:
			newcigars.append(cigar)
			lasttype = thetype
			lastsize = size
			
			
	return lrstart,rrstart,"".join(newcigars)


def cigartoscore(cigars):
	
	allcigars = re.findall(r'\d+[A-Z=]',cigars)
	
	score = 0
	
	match = 0
	mismatch = 0
	sv = 0
	for cigar in allcigars:
		
		size = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype in ['=','M']:
			
			score += size
			match += size
		elif thetype in ['X']:
			
			score -= 20*size 
			mismatch -= 20*size
		elif thetype in ['D']:
			
			score -= (size + 20)
			
			sv -=  (size//2+ 20)
		elif thetype in ['I']:
			score -= (size//2 + 20)
			sv -= 20
			
	#print("score",match, mismatch,sv,score)
	return score

def cigar_findqrange(cigars, qstart , qend):
	
	cigars = re.findall(r'\d+[HSMIDX=][A-Za-z]*', cigars)
	
	
	refposi = 0
	posi = 0
	newcigars = []
	refstart = -1
	refend = -1
	for cigar in cigars:
		
		lastposi = posi
		
		thesize = re.findall(r'\d+',cigar)[0]
		thetype = cigar[len(thesize):][0]
		seq = cigar[(len(thesize)+1):]
		thesize = int(thesize)
		
		if thetype == "=" or thetype == "M":
			
			posi += thesize
			refposi += thesize
			
			
		elif thetype == "X" or thetype == "S":
			
			posi += thesize
			refposi += thesize
			
		elif thetype == "D" or thetype == "H":
			
			refposi += thesize
			
		elif thetype == "I":
			
			posi += thesize
			
		if posi <= qstart:
			
			continue
		
		if lastposi <= qstart:
			
			if thetype in ['M','=', 'X','S', 'D','H']:
				
				refstart = refposi - thesize + qstart - lastposi
			else:
				refstart = refposi
				
			corrsize = min(qend,posi) - qstart
			seq = ""
			if len(seq):
				seq = seq[(qstart-lastposi):(min(qend,posi)-lastposi)]
			newcigars.append("{}{}{}".format(corrsize,thetype,seq))
			
			
		if posi >= qend:
			
			if thetype in ['M','=', 'X','S', 'D','H']:
				
				refend = refposi - posi + qend
			else:
				refend = refposi
				
			if lastposi > qstart:
				corrsize = qend - max(qstart, lastposi)
				
				seq = ""
				if len(seq):
					seq = seq[:corrsize]
				
				newcigars.append("{}{}{}".format(corrsize,thetype,seq))
				
			break
		
		if lastposi > qstart:
			
			newcigars.append(cigar)
	
	
			
	return  [refstart, refend, "".join(newcigars)]

def main(args):

	genename = ""
	misssizes = {}
	genesize = {}	
	geneorder = {}
	generange = {}

	queries = args.query.split(",")


	f = open(queries[0], mode = 'r')
	for line in f:
		if len(line) and line[0] == ">":
			line = line.split()
			locrange = line[1].split(":")[1][:-1].split('-')
		
			loc_s = int("-".join(locrange[:-1]))
			loc_e = int(locrange[-1])
			generange[line[0][1:]] = [loc_s,loc_e]
			geneorder[line[0][1:]] = len(geneorder)
			genename = line[0][1:]
			misssizes[genename] = [loc_s, - loc_e + max(0, loc_s)  ] 
		else:
			misssizes[genename][1] += len(line.strip())
	f.close()

	for query in queries[1:]:
		f = open(query, mode = 'r')
		for line in f:
			if len(line) and line[0] == ">":

				line = line.split()
				locrange = line[1].split(":")[1][:-1].split('-')
		
				loc_s = int("-".join(locrange[:-1]))
				loc_e = int(locrange[-1])

				genename = line[0][1:]
				misssizes[genename] =  [ loc_s, - loc_e + max(0, loc_s) ]  
			else:
				misssizes[genename][1] += len(line.strip())
		f.close()

	
	globalalign = {}
	with open(args.align, mode = 'r') as f:
		for line in f:
			if len(line):
				line = line.split()
				globalalign[line[0]] = line[1:]
	
	output = []
	with open(args.raw, mode = 'r') as f:
		for line in f:
			if len(line.strip()) == 0:
				continue

			line = line.split()
			
			name,chrom,chrom_s,chrom_e, chrom_strand,mapchrom, mapchrom_s, mapchrom_e, mapchrom_strand,maplen, mapscore,rawcigar = line[:12]

			chrom_s,chrom_e,mapchrom_s, mapchrom_e = int(chrom_s),int(chrom_e),int(mapchrom_s), int(mapchrom_e)
			
			gene_s, gene_e = generange[name]
			gene_s, gene_e = int(gene_s), int(gene_e)
		
			gene_eold, chrom_eold,mapchrom_eold = gene_e, chrom_e,mapchrom_e

			score, cigar = globalalign[line[0]]

			qsize = sum([int(x[:-1]) for x in re.findall(r'\d+[MXI]',cigar)])
			rsize = sum([int(x[:-1]) for x in re.findall(r'\d+[MXD]',cigar)])


			anchormiss = misssizes[name]


			#gene_e = min(gene_e,gene_s + genesize -1) 

			mapchrom_e = min(mapchrom_e, mapchrom_s + rsize)
	
			if chrom_strand == "+":
				relative_s = gene_s - chrom_s
				relative_e = gene_e - chrom_s
			else:
				relative_s = chrom_e - gene_e   
				relative_e = chrom_e - gene_s

			rstart_rela, rend_rela, cigar_query = cigar_findqrange(cigar, relative_s , relative_e )

			start, end, polishedcigar = polishstartend(cigar_query)
	
			rstart_rela += start
			rend_rela -= end
			
			if mapchrom_strand == "+":
				globalstart = mapchrom_s + rstart_rela
				globalend = mapchrom_s + rend_rela
			else:
				globalstart = mapchrom_e - rend_rela
				globalend =  mapchrom_e - rstart_rela
		
			score = cigartoscore(polishedcigar)

			breaks, graphcigar =cigar_findrrange(rawcigar, relative_s , relative_e )

			min_break = min(breaks)
			max_break = max(breaks)

			strdscore = sum([breaks[2*i +1] - breaks[2*i] for i in range(len(breaks)//2)])

			strd_breaks = '+' if strdscore >= 0 else '-'

			if chrom_strand == '-':
				strd_breaks = '+' if strd_breaks =='-' else '-'

			graphscore = cigartoscore(graphcigar)

			globalscore = cigartoscore(cigar)

			line = [line[0]] + [chrom ,chrom_strand,gene_s, gene_e,  abs(gene_s - gene_e), mapchrom, mapchrom_strand, globalstart, globalend,  abs(globalstart- globalend), score, strd_breaks, min_break, max_break, max_break- min_break, graphscore] + [polishedcigar, graphcigar] + [chrom,chrom_s,chrom_e,mapchrom, mapchrom_s, mapchrom_e, mapchrom_strand, maplen, mapscore, globalscore, int(anchormiss[0] > 0) + int(anchormiss[1] > 0) ] 
			
			output.append(line)
	
	output = sorted(output, key = lambda x:geneorder[x[0]])
	
	with open(args.output, mode = 'w') as f:
		for line in output:
			f.write("\t".join(list(map(str,line))) + '\n')
	
	with open(args.output+"_summary.txt", mode = 'w') as f:
		for line in output:
			name, chrom ,chrom_strand,gene_s, gene_e,  genesize, mapchrom, mapchrom_strand, globalstart, globalend,  galignsize, globalscore, strd_breaks, min_break, max_break, graphsize, graphscore, polishedcigar, graphcigar , chrom,chrom_s,chrom_e,mapchrom, mapchrom_s, mapchrom_e, mapchrom_strand, maplen, mapscore, rawglobalscore, enddis = line 


			alignloc = []
			#if globalscore > 0.9 * galignsize or globalscore >= graphscore:
			alignloc = [mapchrom, mapchrom_strand, globalstart, globalend, globalscore, maplen, mapscore, enddis, polishedcigar]
			#else:
				#alignloc = [mapchrom, strd_breaks, min_break, max_break,graphscore, maplen, mapscore, enddis, graphcigar]

			line = "\t".join(list(map(str, [name, chrom ,chrom_strand,gene_s, gene_e] + alignloc   )))
	
			f.write(line+ "\n")


def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	
	parser.add_argument("-r", "--raw", help="raw location", dest="raw", type=str,required=True)
	parser.add_argument("-q", "--query", help="path to output file", dest="query", type=str,required=True)
	parser.add_argument("-a", "--align", help="path to output file", dest="align", type=str,required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
