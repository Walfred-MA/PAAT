#!/usr/bin/env python3


import collections as cl
import pandas as pd
import re
import argparse
import os

alterchrom = set(['chr10_GL383545v1_alt', 'chr10_GL383546v1_alt', 'chr10_KI270824v1_alt', 'chr10_KI270825v1_alt', 'chr10_KN196480v1_fix', 'chr10_KQ090021v1_fix', 'chr10_ML143354v1_fix', 'chr10_ML143355v1_fix', 'chr11_JH159136v1_alt', 'chr11_JH159137v1_alt', 'chr11_KI270721v1_random', 'chr11_KI270829v1_alt', 'chr11_KI270830v1_alt', 'chr11_KI270831v1_alt', 'chr11_KI270832v1_alt', 'chr11_KI270902v1_alt', 'chr11_KI270903v1_alt', 'chr11_KI270927v1_alt', 'chr11_KN196481v1_fix', 'chr11_KN538368v1_alt', 'chr11_KQ090022v1_fix', 'chr11_KQ759759v1_fix', 'chr11_KV766195v1_fix', 'chr11_KZ559108v1_fix', 'chr11_KZ559109v1_fix', 'chr11_KZ559110v1_alt', 'chr11_KZ559111v1_alt', 'chr11_ML143357v1_fix', 'chr11_ML143358v1_fix', 'chr11_ML143359v1_fix', 'chr11_ML143360v1_fix', 'chr12_GL383550v2_alt', 'chr12_GL383551v1_alt', 'chr12_GL383552v1_alt', 'chr12_GL383553v2_alt', 'chr12_GL877875v1_alt', 'chr12_GL877876v1_alt', 'chr12_KI270833v1_alt', 'chr12_KI270834v1_alt', 'chr12_KI270835v1_alt', 'chr12_KI270837v1_alt', 'chr12_KI270904v1_alt', 'chr12_KN538369v1_fix', 'chr12_KN538370v1_fix', 'chr12_KQ759760v1_fix', 'chr12_KZ208916v1_fix', 'chr12_KZ208917v1_fix', 'chr12_KZ559112v1_alt', 'chr12_ML143361v1_fix', 'chr12_ML143362v1_fix', 'chr13_KI270838v1_alt', 'chr13_KI270840v1_alt', 'chr13_KI270842v1_alt', 'chr13_KN538371v1_fix', 'chr13_KN538372v1_fix', 'chr13_KN538373v1_fix', 'chr13_ML143365v1_fix', 'chr13_ML143366v1_fix', 'chr14_GL000009v2_random', 'chr14_GL000194v1_random', 'chr14_GL000225v1_random', 'chr14_KI270726v1_random', 'chr14_KI270844v1_alt', 'chr14_KI270845v1_alt', 'chr14_KI270846v1_alt', 'chr14_KI270847v1_alt', 'chr14_KZ208919v1_alt', 'chr14_KZ208920v1_fix', 'chr14_ML143367v1_fix', 'chr15_GL383554v1_alt', 'chr15_GL383555v2_alt', 'chr15_KI270727v1_random', 'chr15_KI270848v1_alt', 'chr15_KI270849v1_alt', 'chr15_KI270850v1_alt', 'chr15_KI270851v1_alt', 'chr15_KI270852v1_alt', 'chr15_KI270905v1_alt', 'chr15_KI270906v1_alt', 'chr15_KQ031389v1_alt', 'chr15_ML143369v1_fix', 'chr15_ML143370v1_fix', 'chr15_ML143371v1_fix', 'chr15_ML143372v1_fix', 'chr16_GL383556v1_alt', 'chr16_GL383557v1_alt', 'chr16_KI270728v1_random', 'chr16_KI270853v1_alt', 'chr16_KI270854v1_alt', 'chr16_KI270855v1_alt', 'chr16_KI270856v1_alt', 'chr16_KQ090026v1_alt', 'chr16_KQ090027v1_alt', 'chr16_KV880768v1_fix', 'chr16_KZ208921v1_alt', 'chr16_KZ559113v1_fix', 'chr16_ML143373v1_fix', 'chr17_GL000205v2_random', 'chr17_GL000258v2_alt', 'chr17_GL383563v3_alt', 'chr17_GL383564v2_alt', 'chr17_GL383565v1_alt', 'chr17_GL383566v1_alt', 'chr17_JH159146v1_alt', 'chr17_JH159147v1_alt', 'chr17_JH159148v1_alt', 'chr17_KI270857v1_alt', 'chr17_KI270858v1_alt', 'chr17_KI270859v1_alt', 'chr17_KI270860v1_alt', 'chr17_KI270861v1_alt', 'chr17_KI270862v1_alt', 'chr17_KI270907v1_alt', 'chr17_KI270908v1_alt', 'chr17_KI270909v1_alt', 'chr17_KI270910v1_alt', 'chr17_KV575245v1_fix', 'chr17_KV766196v1_fix', 'chr17_KV766198v1_alt', 'chr17_KZ559114v1_alt', 'chr17_ML143374v1_fix', 'chr17_ML143375v1_fix', 'chr18_GL383567v1_alt', 'chr18_GL383569v1_alt', 'chr18_GL383570v1_alt', 'chr18_GL383571v1_alt', 'chr18_GL383572v1_alt', 'chr18_KI270863v1_alt', 'chr18_KI270911v1_alt', 'chr18_KI270912v1_alt', 'chr18_KQ090028v1_fix', 'chr18_KQ458385v1_alt', 'chr18_KZ208922v1_fix', 'chr18_KZ559115v1_fix', 'chr18_KZ559116v1_alt', 'chr19_GL000209v2_alt', 'chr19_GL383573v1_alt', 'chr19_GL383574v1_alt', 'chr19_GL383575v2_alt', 'chr19_GL383576v1_alt', 'chr19_GL949746v1_alt', 'chr19_GL949747v2_alt', 'chr19_GL949748v2_alt', 'chr19_GL949749v2_alt', 'chr19_GL949750v2_alt', 'chr19_GL949751v2_alt', 'chr19_GL949752v1_alt', 'chr19_GL949753v2_alt', 'chr19_KI270865v1_alt', 'chr19_KI270866v1_alt', 'chr19_KI270867v1_alt', 'chr19_KI270868v1_alt', 'chr19_KI270882v1_alt', 'chr19_KI270883v1_alt', 'chr19_KI270884v1_alt', 'chr19_KI270885v1_alt', 'chr19_KI270886v1_alt', 'chr19_KI270887v1_alt', 'chr19_KI270888v1_alt', 'chr19_KI270889v1_alt', 'chr19_KI270890v1_alt', 'chr19_KI270891v1_alt', 'chr19_KI270914v1_alt', 'chr19_KI270915v1_alt', 'chr19_KI270916v1_alt', 'chr19_KI270917v1_alt', 'chr19_KI270918v1_alt', 'chr19_KI270919v1_alt', 'chr19_KI270920v1_alt', 'chr19_KI270921v1_alt', 'chr19_KI270922v1_alt', 'chr19_KI270923v1_alt', 'chr19_KI270929v1_alt', 'chr19_KI270930v1_alt', 'chr19_KI270931v1_alt', 'chr19_KI270932v1_alt', 'chr19_KI270933v1_alt', 'chr19_KI270938v1_alt', 'chr19_KN196484v1_fix', 'chr19_KQ458386v1_fix', 'chr19_KV575246v1_alt', 'chr19_KV575247v1_alt', 'chr19_KV575248v1_alt', 'chr19_KV575249v1_alt', 'chr19_KV575250v1_alt', 'chr19_KV575251v1_alt', 'chr19_KV575252v1_alt', 'chr19_KV575253v1_alt', 'chr19_KV575254v1_alt', 'chr19_KV575255v1_alt', 'chr19_KV575256v1_alt', 'chr19_KV575257v1_alt', 'chr19_KV575258v1_alt', 'chr19_KV575259v1_alt', 'chr19_KV575260v1_alt', 'chr19_ML143376v1_fix', 'chr1_GL383518v1_alt', 'chr1_GL383519v1_alt', 'chr1_GL383520v2_alt', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270759v1_alt', 'chr1_KI270760v1_alt', 'chr1_KI270761v1_alt', 'chr1_KI270762v1_alt', 'chr1_KI270763v1_alt', 'chr1_KI270765v1_alt', 'chr1_KI270766v1_alt', 'chr1_KI270892v1_alt', 'chr1_KN196472v1_fix', 'chr1_KN196473v1_fix', 'chr1_KN196474v1_fix', 'chr1_KN538360v1_fix', 'chr1_KN538361v1_fix', 'chr1_KQ031383v1_fix', 'chr1_KQ458382v1_alt', 'chr1_KQ458383v1_alt', 'chr1_KQ458384v1_alt', 'chr1_KQ983255v1_alt', 'chr1_KV880763v1_alt', 'chr1_KZ208904v1_alt', 'chr1_KZ208905v1_alt', 'chr1_KZ208906v1_fix', 'chr20_GL383577v2_alt', 'chr20_KI270869v1_alt', 'chr20_KI270870v1_alt', 'chr20_KI270871v1_alt', 'chr21_GL383579v2_alt', 'chr21_GL383580v2_alt', 'chr21_GL383581v2_alt', 'chr21_KI270872v1_alt', 'chr21_KI270873v1_alt', 'chr21_KI270874v1_alt', 'chr21_ML143377v1_fix', 'chr22_GL383582v2_alt', 'chr22_GL383583v2_alt', 'chr22_KB663609v1_alt', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270875v1_alt', 'chr22_KI270876v1_alt', 'chr22_KI270877v1_alt', 'chr22_KI270878v1_alt', 'chr22_KI270879v1_alt', 'chr22_KI270928v1_alt', 'chr22_KN196485v1_alt', 'chr22_KN196486v1_alt', 'chr22_KQ458387v1_alt', 'chr22_KQ458388v1_alt', 'chr22_KQ759761v1_alt', 'chr22_KQ759762v1_fix', 'chr22_ML143378v1_fix', 'chr22_ML143380v1_fix', 'chr2_GL383521v1_alt', 'chr2_GL383522v1_alt', 'chr2_GL582966v2_alt', 'chr2_KI270767v1_alt', 'chr2_KI270768v1_alt', 'chr2_KI270769v1_alt', 'chr2_KI270770v1_alt', 'chr2_KI270772v1_alt', 'chr2_KI270773v1_alt', 'chr2_KI270774v1_alt', 'chr2_KI270776v1_alt', 'chr2_KI270893v1_alt', 'chr2_KI270894v1_alt', 'chr2_KN538362v1_fix', 'chr2_KN538363v1_fix', 'chr2_KQ031384v1_fix', 'chr2_KQ983256v1_alt', 'chr2_KZ208907v1_alt', 'chr2_KZ208908v1_alt', 'chr2_ML143341v1_fix', 'chr2_ML143342v1_fix', 'chr3_GL383526v1_alt', 'chr3_JH636055v2_alt', 'chr3_KI270777v1_alt', 'chr3_KI270779v1_alt', 'chr3_KI270780v1_alt', 'chr3_KI270781v1_alt', 'chr3_KI270782v1_alt', 'chr3_KI270783v1_alt', 'chr3_KI270784v1_alt', 'chr3_KI270895v1_alt', 'chr3_KI270924v1_alt', 'chr3_KI270934v1_alt', 'chr3_KI270935v1_alt', 'chr3_KI270936v1_alt', 'chr3_KI270937v1_alt', 'chr3_KN196475v1_fix', 'chr3_KN196476v1_fix', 'chr3_KN538364v1_fix', 'chr3_KQ031385v1_fix', 'chr3_KV766192v1_fix', 'chr3_KZ208909v1_alt', 'chr3_KZ559102v1_alt', 'chr3_KZ559103v1_alt', 'chr3_KZ559104v1_fix', 'chr3_KZ559105v1_alt', 'chr3_ML143343v1_alt', 'chr4_GL000257v2_alt', 'chr4_GL383527v1_alt', 'chr4_GL383528v1_alt', 'chr4_KI270785v1_alt', 'chr4_KI270786v1_alt', 'chr4_KI270788v1_alt', 'chr4_KI270789v1_alt', 'chr4_KI270790v1_alt', 'chr4_KI270896v1_alt', 'chr4_KI270925v1_alt', 'chr4_KQ090014v1_alt', 'chr4_KQ090015v1_alt', 'chr4_KQ983257v1_fix', 'chr4_KQ983258v1_alt', 'chr4_KV766193v1_alt', 'chr4_ML143344v1_fix', 'chr4_ML143345v1_fix', 'chr4_ML143347v1_fix', 'chr4_ML143349v1_fix', 'chr5_GL339449v2_alt', 'chr5_GL383531v1_alt', 'chr5_GL383532v1_alt', 'chr5_GL949742v1_alt', 'chr5_KI270791v1_alt', 'chr5_KI270792v1_alt', 'chr5_KI270793v1_alt', 'chr5_KI270794v1_alt', 'chr5_KI270795v1_alt', 'chr5_KI270796v1_alt', 'chr5_KI270897v1_alt', 'chr5_KI270898v1_alt', 'chr5_KN196477v1_alt', 'chr5_KV575243v1_alt', 'chr5_KV575244v1_fix', 'chr5_ML143350v1_fix', 'chr6_GL000250v2_alt', 'chr6_GL000251v2_alt', 'chr6_GL000252v2_alt', 'chr6_GL000253v2_alt', 'chr6_GL000254v2_alt', 'chr6_GL000255v2_alt', 'chr6_GL000256v2_alt', 'chr6_GL383533v1_alt', 'chr6_KB021644v2_alt', 'chr6_KI270758v1_alt', 'chr6_KI270797v1_alt', 'chr6_KI270798v1_alt', 'chr6_KI270799v1_alt', 'chr6_KI270800v1_alt', 'chr6_KI270801v1_alt', 'chr6_KI270802v1_alt', 'chr6_KN196478v1_fix', 'chr6_KQ031387v1_fix', 'chr6_KQ090016v1_fix', 'chr6_KV766194v1_fix', 'chr6_KZ208911v1_fix', 'chr7_GL383534v2_alt', 'chr7_KI270803v1_alt', 'chr7_KI270804v1_alt', 'chr7_KI270805v1_alt', 'chr7_KI270806v1_alt', 'chr7_KI270807v1_alt', 'chr7_KI270808v1_alt', 'chr7_KI270809v1_alt', 'chr7_KI270899v1_alt', 'chr7_KQ031388v1_fix', 'chr7_KV880764v1_fix', 'chr7_KV880765v1_fix', 'chr7_KZ208912v1_fix', 'chr7_KZ208913v1_alt', 'chr7_KZ559106v1_alt', 'chr7_ML143352v1_fix', 'chr8_KI270810v1_alt', 'chr8_KI270811v1_alt', 'chr8_KI270812v1_alt', 'chr8_KI270813v1_alt', 'chr8_KI270814v1_alt', 'chr8_KI270815v1_alt', 'chr8_KI270816v1_alt', 'chr8_KI270817v1_alt', 'chr8_KI270818v1_alt', 'chr8_KI270819v1_alt', 'chr8_KI270821v1_alt', 'chr8_KI270822v1_alt', 'chr8_KI270900v1_alt', 'chr8_KI270901v1_alt', 'chr8_KI270926v1_alt', 'chr8_KV880766v1_fix', 'chr8_KZ208914v1_fix', 'chr8_KZ208915v1_fix', 'chr8_KZ559107v1_alt', 'chr9_GL383539v1_alt', 'chr9_GL383540v1_alt', 'chr9_GL383541v1_alt', 'chr9_GL383542v1_alt', 'chr9_KI270823v1_alt', 'chr9_KN196479v1_fix', 'chr9_KQ090018v1_alt', 'chrUn_GL000195v1', 'chrUn_GL000213v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_KI270442v1', 'chrUn_KI270744v1', 'chrUn_KI270750v1', 'chrX_KI270880v1_alt', 'chrX_KI270881v1_alt', 'chrX_KI270913v1_alt', 'chrX_ML143381v1_fix', 'chrY_KN196487v1_fix', 'chrY_KZ208924v1_fix'])

def reversecigar(cigars):

	return  "".join(["".join(['>' if chunk[0] == '<' else '<']+ re.findall(r'\d+[a-zA-Z=][actgACTG]*', chunk)[::-1]) for chunk in cigars.replace("<","&<").replace(">","&>").split("&")[1:]][::-1])

def getlocations(inputfile, queryfile, rawloc, forcemain = 1):

	qsize = cl.defaultdict(int)	
	qname = ""
	qinfo = {}

	if os.path.isfile(queryfile):
		f = open(queryfile, mode = 'r')
		
		for line in f:
			
			if len(line) == 0 or line[0] != '>':

				qsize[qname] += len(line.strip())
				continue
			
			line = line.strip().split()
			
			qname = line[0][1:]
			qinfo [qname] = line[1]
		f.close()
		

	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			line = line.split()
			
			score = int(line[6])
			#score -= abs(int(line[1]) - (int(line[4])-int(line[3])))
			
			if "_" in line[2] and "NC_0609" not in line[2]  and "_h" not in line[0]:
				
				continue
		
			if score <= 0:
				continue

			if forcemain or ( "v" in line[0].split("_")[2] or "NC_0609"  in line[0] ):	
				if "NC_0609" in line[2]:
					score *= 0.000001
				elif "_"  in line[2] or "v" in line[2]:
					score *= 0
			else:
				if "NC_0609" in line[2] :
					score *= 1.0
				elif "_"  in line[2] or "v" in line[2]:
					score *= 0.9			

					
			size = int(line[1])

			qname = line[0]

			qstrd = qinfo[line[0]][-1]

			if qstrd == '-':
				line[-1] = reversecigar(line[-1])

			rawloc [line[0]].append([(score), [line[0],line[2], line[3], line[4], line[5], line[1],line[6], line[-1]], [qinfo[line[0]],qsize[qname] ] ])

	return rawloc

def main(args):
	
	inputs = args.input.split(",")
	queries = args.query.split(",")
	
	rawloc = cl.defaultdict(list)
	for infile, queryfile in zip(inputs, queries):
	
		rawloc = getlocations(infile, queryfile, rawloc, args.main)
	
	with open(args.output, mode = 'w') as f:
		for query, locs in rawloc.items():
			
			locs = sorted(locs, key = lambda x: x[0],reverse  = True)
			qname = locs[0][1][0]
			size = locs[0][1][5]
				
			qinfo = locs[0][2][0]
			qsize = locs[0][2][1]
			qstrd = qinfo[-1]
			info = qinfo [:-1]
			qcontig = qinfo.split(":")[0]
			qrange = qinfo.split(":")[1]
			qstart = 0 if qrange[0] == '-' else  int(qrange.split("-")[0])

			#qend = int(qrange.split("-")[:-1])
			qend = qstart + qsize
	
			rstrd = locs[0][1][4]
			
			if qstrd == '-':
				rstrd = '-' if rstrd == '+' else '+'
				
			locs[0][1][4] = rstrd

			
			locs[0][1][5] = str(qend - qstart + 1)
			out = [qname,qcontig, str(qstart), str(qend), qstrd] + locs[0][1][1:]   
			
			f.write("\t".join(out)+"\n")
			
			
			
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-q", "--query", help="path to input data file",dest="query", type=str, required=True)
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
	parser.add_argument("-m", "--main", help="path to output file", dest="main", type=int, default = 1)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
