#!/usr/bin/env python3

import os
import argparse
import re
import bisect


class kmeraligns:
	
	def __init__(self, cigars, ksize = 30, cutoff = 50):
		
		self.cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
		self.ksize = ksize
		self.cutoff = cutoff
		self.min = -cutoff
		self.max = cutoff*3
		self.consensus_f = []
		self.consensus_r = []
	
	def getalignsegs(self):
		
		lastzero = 0
		self.consensus_f = [[0,1]]
		rposi = 0
		qposi = 0 
		lastvar = 0
		lastscore = 0
		score = 0 
		currseg = 1
		for index, cigar in enumerate(self.cigars):
			
			lastscore = score
			thesize = int(cigar[:-1])
			thetype = cigar[-1]
			if thetype == "=" or thetype == "M":
				
				qposi += thesize
				rposi += thesize
					
				score += thesize - self.ksize
				score = min(self.max,max(self.min, score))
				
			elif thetype == "X":
				
				qposi += thesize
				rposi += thesize
					
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			elif thetype == "H":
				
				rposi += thesize
				continue
			
			elif thetype == "I" :
				
				rposi += thesize
				
				score -= ( min(self.ksize, qposi - lastvar)  )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			elif thetype == "D":
				
				qposi += thesize
					
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			if lastscore * score < 0 or score == 0:
				
				lastzero = index
				
			if currseg <= 0 and score >= self.cutoff:
				
				self.consensus_f.append([lastzero,lastzero])
				currseg = 1
				
			if currseg > 0 and score <= 0:
				
				
				self.consensus_f[-1][1] = index 
				currseg = -1
				
		if currseg > 0:
			if score >= lastscore:
				self.consensus_f[-1][1] = index + 1
			else:
				self.consensus_f[-1][1] = index 
		
		self.consensus_f = [x for x in self.consensus_f if x[1] > x[0]]
		
		lastzero = len(self.cigars)
		self.consensus_r = [[lastzero-1,lastzero]]
		rposi = 0
		qposi = 0 
		lastvar = 0
		lastscore = 0
		score = 0 
		currseg = 1
		for index, cigar in enumerate(self.cigars[::-1]):
			
			index = len(self.cigars) - index
			lastscore = score
			thesize = int(cigar[:-1])
			thetype = cigar[-1]
			
			if thetype == "=" or thetype == "M":
				
				qposi += thesize
				rposi += thesize
					
				score += thesize - self.ksize
				score = min(self.max,max(self.min, score))
				
			elif thetype == "X":
				
				qposi += thesize
				rposi += thesize
				
				if score > 0 and thesize > self.cutoff:
					score = 0
					
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			elif thetype == "H":
				
				rposi += thesize
				continue
			
			elif thetype == "I" :
				
				rposi += thesize
				
				score -= ( min(self.ksize, qposi - lastvar)  )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			elif thetype == "D":
				
				qposi += thesize
				
				if score > 0 and thesize > self.cutoff:
					score = 0
					
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = min(self.max,max(self.min, score))
				
				lastvar = qposi
				
			if lastscore * score < 0 or score == 0:
				
				lastzero = index
				
			if currseg <= 0 and score >= self.cutoff:
				
				self.consensus_r.append([lastzero,lastzero])
				currseg = 1
				
			if currseg > 0 and score <= 0:
				
				self.consensus_r[-1][0] = index 
				currseg = -1
		
		if currseg > 0:
			if score >= lastscore:
				self.consensus_r[-1][0] = 0
			else:
				self.consensus_r[-1][0] = 1
		
		self.consensus_r = [x for x in self.consensus_r if x[1] > x[0]]
				
		return self
		
	def combineseg(self):
		
		consensus = self.consensus_f + self.consensus_r
		consensus_pos = [a for b in consensus for a in b]
		
		sortindexes = sorted(range(len(consensus_pos)), key = lambda x: consensus_pos[x])
		
		num_curr = 0
		self.consensus_comb = []
		for index, sortindex in enumerate(sortindexes):
			
			if sortindex % 2 == 0:
				num_curr += 1
				
				if num_curr == 1:
					
					gapsize = 10000000
					if len(self.consensus_comb):
						lastend = self.consensus_comb[-1][-1]
						gapsize = sum([int(self.cigars[x][:-1]) for x in range(lastend,sortindex) if self.cigars[x][-1]])
					
					if gapsize > self.cutoff:
						self.consensus_comb.append([consensus_pos[sortindex], consensus_pos[sortindex]])
					
			else:
				num_curr -= 1
				
				if num_curr == 0:
					
					self.consensus_comb[-1][-1] = consensus_pos[sortindex]
		
		return self
	
	def trim(self):
		
		self.consensus_trimed = []
		for segment in self.consensus_comb:
			
			lastvar = 0
			score = 0
			lastzero = 0
			rposi = 0
			qposi = 0 
			lastvar = 0
			lindex = segment[1]-segment[0]
			rindex = segment[1]-segment[0]
						
			for index in range(segment[0],segment[1]):
				
				cigar = self.cigars[index] 
				thesize = int(cigar[:-1])
				thetype = cigar[-1]
				
				if thetype == "=" or thetype == "M":
					
					if score < 0 and thesize > self.cutoff:
						score = 0
						
					score += thesize - self.ksize
					score = min(self.max,max(self.min, score))
					
					if score > self.ksize:
						lindex = index
						break
					
				elif thetype == "X":
					
					if score > 0 and thesize > self.cutoff:
						score = 0
						
					score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
					
				elif thetype == "H":
					
					rposi += thesize
					continue
				
				elif thetype == "I" :
					
					rposi += thesize
					
					score -= ( min(self.ksize, qposi - lastvar)  )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
					
				elif thetype == "D":
					
					qposi += thesize
						
					score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
			
			score = 0
			lastzero = 0
			rposi = 0
			qposi = 0 
			lastvar = 0
			
			for index in list(range(segment[0],segment[1]))[::-1]:
				
				cigar = self.cigars[index] 
				thesize = int(cigar[:-1])
				thetype = cigar[-1]
				
				if thetype == "=" or thetype == "M":
					
					if score < 0 and thesize > self.cutoff:
						score = 0
						
					score += thesize - self.ksize
					score = min(self.max,max(self.min, score))
					
					if score >= self.ksize:
						rindex = index
						break
					
				elif thetype == "X":
					
					score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
					
				elif thetype == "H":
					
					rposi += thesize
					continue
				
				elif thetype == "I" :
					
					rposi += thesize
					
					score -= ( min(self.ksize, qposi - lastvar)  )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
					
				elif thetype == "D":
					
					qposi += thesize
						
					score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
					score = min(self.max,max(self.min, score))
					
					lastvar = qposi
						
			if rindex >= lindex:
				self.consensus_trimed.append([lindex, rindex+1])
			
		return self
	def coordinate(self):
		
		self.segment_coordinates = [[] for x in self.consensus_comb]
		
		allstartends = {}
		for index,segment in enumerate(self.consensus_trimed):
			
			allstartends[2*segment[0]] = index
			allstartends[2*segment[1]-1] = index
			
		
		peakscore = 0
		score = 0
		rposi = 0
		qposi = 0 
		lastvar = 0
		for index, cigar in enumerate(self.cigars):
			
			if 2*index in allstartends:
				self.segment_coordinates[allstartends[2*index]].append(qposi)
				score = 0
				peakscore = 0
			
			thesize = int(cigar[:-1])
			thetype = cigar[-1]
			
			if thetype == "=" or thetype == "M":
				
				qposi += thesize
				rposi += thesize
				score += thesize - self.ksize
				score = max(self.min, score)
				peakscore = max(peakscore, score)
				
			elif thetype == "X":
				
				qposi += thesize
				rposi += thesize
				score -= ( min(self.ksize, qposi - lastvar) + thesize - 1 )
				score = max(self.min, score)
				
				lastvar = qposi
				
			elif thetype == "H":
				
				rposi += thesize
				continue
			
			elif thetype == "I" :
				
				rposi += thesize
				
				score -= ( min(self.ksize, qposi - lastvar)  )
				score = max(self.min, score)
				
				lastvar = qposi
				
			elif thetype == "D":
				
				qposi += thesize
				
				score -= ( min(self.ksize, qposi - lastvar)  )
				score = max(self.min, score)
				
				lastvar = qposi
				
			if 2*index + 1 in allstartends:
				
				segindex = allstartends[2*index + 1]
								
				self.segment_coordinates[segindex].append(qposi)
				self.segment_coordinates[segindex].append(peakscore)
				cigar_segs = "".join([self.cigars[x] for x in range(*self.consensus_trimed[segindex])])
				
				self.segment_coordinates[segindex].append(cigar_segs)
		
		return self
	

def findkmeraligns(cigars, ksize = 30, cutoff = 100):
	
	consensus = kmeraligns(cigars, ksize, cutoff).getalignsegs().combineseg().trim().coordinate()
	
	return  consensus.segment_coordinates

def makereverse(seq):
	
	tran=str.maketrans('ATCGatcg', 'TAGCtagc')
	
	return seq[::-1].translate(tran)

def pathtoseq(path, pathstrs):
	
	pathes = re.findall(r'[><]s\d+', path)
	
	fullseq = ""
	for eachpath in pathes:
		
		strd = eachpath[0]
		
		seq = pathstrs[eachpath[1:]]
		
		if strd == "<":
			
			seq = makereverse(seq)
			
		fullseq += seq
		
	return fullseq

def alignregion(allaligns, allow_gap = 0):


	all_coordinates = [x for y in allaligns for x in y[:2]]
	scores = [y[2] for y in allaligns]
	lineindex = [y[3] for y in allaligns]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	current_alignscores = []
	allgroups = [[0,0,0,0]]
		
	number_curr_seq = 0
	last_coordinate = - allow_gap - 1
	last_index = -2
	last_score = -1
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		score = scores[index//2]


		if index %2 == 0 :

			last_score = current_alignscores[-1][0] if len(current_alignscores) else 0
			if last_index%2 and len(current_alignscores) and coordinate - last_coordinate>allow_gap:

				current_alignscores.remove([scores[last_index//2], last_index//2])

			
			if  len(current_alignscores) == 0  or score > last_score :
				
				if number_curr_seq:
					allgroups[-1][-1] = coordinate
					
				allgroups.append([score, index//2, coordinate, coordinate])
				
				
			bisect.insort(current_alignscores, [score, index//2])
		
	
			number_curr_seq += 1
			
		else:
			
			allgroups[-1][-1] = coordinate
			
			number_curr_seq -= 1
		
		last_index =  index
	
		last_coordinate = coordinate

	allgroups = [x[2:]+allaligns[x[1]][2:]+[lineindex[x[1]]] for x in allgroups[1:]]

	return allgroups


def cigaredit(cigars, queryseq, refseq):
	
	cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
	
	refposi = 0
	posi = 0
	newcigars = []
	
	newref = ""
	for cigar in cigars:
		
		thesize = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == "=":
			
			
			matchormismatch = 0
			last_i = 0
			newcigar = []
			for i, (r,q) in enumerate(zip(refseq[refposi:(refposi+thesize)], queryseq[posi:(posi+thesize)])):

				if matchormismatch ==0 and r.upper() != q.upper():

					newcigar.append("={}".format(i-last_i))
					last_i = i
					matchormismatch = 1
					check = 1

				elif matchormismatch ==1 and r.upper() == q.upper():

					newcigar.append("X{}".format(i-last_i))
					last_i = i
					matchormismatch = 0

			i = thesize
			if matchormismatch ==0:
				newcigar.append("={}".format(i-last_i))
			else:
				newcigar.append("X{}".format(i-last_i))

			cigar = "".join(newcigar)
			
			
			newcigars.append(cigar)
			
			newref += refseq[refposi:(refposi+thesize)]
			
			posi += thesize
			
			refposi += thesize
			
			
		elif thetype == "X":
			
			newcigars.append(cigar+queryseq[posi:(posi+thesize)])
			posi += thesize
			
			newref += "X"*thesize
			refposi += thesize
			
		elif thetype == "D":
			
			newcigars.append(cigar)
			
			refposi += thesize
			
		elif thetype == "I":
			
			newcigars.append(cigar+queryseq[posi:(posi+thesize)])
			posi += thesize
			
			newref += "-"*thesize
			
			
	return "".join(newcigars)



def lineartograph(graphfile, queryfile, alignfile, output):
	
	with open(queryfile, mode = 'r') as f:
		reads = [read.splitlines() for read in f.read().split(">")[1:]]
		reads = {read[0].split()[0]:"".join(read[1:]) for read in reads}
		
	with open(graphfile, mode = 'r') as f:
		refs = [read.splitlines() for read in f.read().split(">")[1:]]
		refs = {read[0].split()[0]:"".join(read[1:]) for read in refs}
	
	
	allaligns = []
	with open(alignfile, mode = 'r') as f:
		
		lines = f.read().splitlines()


	for index, line in enumerate(lines):
			
		if len(line) == 0 or line[0]=="@":
			continue
			
		elements = line.strip().split()
			#name, qsize, qstart, qend, path, rsize, rstart, rend = elements[:8]
			
			
		identity =float( [x for i,x in enumerate(elements) if x[:5] == "PI:f:"][0][5:])
		match = int( [x for i,x in enumerate(elements) if x[:5] == "AS:i:"][0][5:] ) 
			
		if identity < 90 or match < 100 or elements[4] != '255':
			continue
			
		strand, qstart, cigar = elements[1], int(elements[3]), elements[5]
		strand = 1 if strand == '0' else -1
			
		aligns = findkmeraligns(cigar, qstart)
			
		allaligns.extend([x+[match, index] for x in aligns])
			
			#newcigar = cigaredit(cigar, queryseq[int(qstart):int(qend)], refseq)
			
			#elements[index] = "{}I{}".format(qstart,queryseq[:int(qstart)])+ newcigar + "{}I{}".format(int(qsize)-int(qend),queryseq[int(qend):])
	
	alignregions = alignregion(allaligns)

	for region in alignregions:

		line = lines [region[0]] 	
		
			
	exit(0)
	
	w= open(output, mode = 'w')
	
	w.write("\t".join(elements)+"\n")
	w.close()
	
def main(args):
	
	lineartograph(args.ref, args.query, args.input, args.output)
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")

	parser.add_argument("-r", "--ref", help="path to output file", dest="ref", type=str,required=True) 
	parser.add_argument("-q", "--query", help="path to output file", dest="query", type=str,required=True) 	
	parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
