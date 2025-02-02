#!/usr/bin/env python3


import argparse
import os
import bisect
import collections as cl
import re
import multiprocessing as mul

manager = mul.Manager()
chm13chrom = ['NC_060925.1', 'NC_060926.1', 'NC_060927.1', 'NC_060928.1', 'NC_060929.1', 'NC_060930.1', 'NC_060931.1', 'NC_060932.1', 'NC_060933.1', 'NC_060934.1', 'NC_060935.1', 'NC_060936.1', 'NC_060937.1', 'NC_060938.1', 'NC_060939.1', 'NC_060940.1', 'NC_060941.1', 'NC_060942.1', 'NC_060943.1', 'NC_060944.1', 'NC_060945.1', 'NC_060946.1', 'NC_060947.1', 'NC_060948.1']
mainchrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrM', 'chrX', 'chrY', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
alterchrom = set(['chr10_GL383545v1_alt', 'chr10_GL383546v1_alt', 'chr10_KI270824v1_alt', 'chr10_KI270825v1_alt', 'chr10_KN196480v1_fix', 'chr10_KQ090021v1_fix', 'chr10_ML143354v1_fix', 'chr10_ML143355v1_fix', 'chr11_JH159136v1_alt', 'chr11_JH159137v1_alt', 'chr11_KI270721v1_random', 'chr11_KI270829v1_alt', 'chr11_KI270830v1_alt', 'chr11_KI270831v1_alt', 'chr11_KI270832v1_alt', 'chr11_KI270902v1_alt', 'chr11_KI270903v1_alt', 'chr11_KI270927v1_alt', 'chr11_KN196481v1_fix', 'chr11_KN538368v1_alt', 'chr11_KQ090022v1_fix', 'chr11_KQ759759v1_fix', 'chr11_KV766195v1_fix', 'chr11_KZ559108v1_fix', 'chr11_KZ559109v1_fix', 'chr11_KZ559110v1_alt', 'chr11_KZ559111v1_alt', 'chr11_ML143357v1_fix', 'chr11_ML143358v1_fix', 'chr11_ML143359v1_fix', 'chr11_ML143360v1_fix', 'chr12_GL383550v2_alt', 'chr12_GL383551v1_alt', 'chr12_GL383552v1_alt', 'chr12_GL383553v2_alt', 'chr12_GL877875v1_alt', 'chr12_GL877876v1_alt', 'chr12_KI270833v1_alt', 'chr12_KI270834v1_alt', 'chr12_KI270835v1_alt', 'chr12_KI270837v1_alt', 'chr12_KI270904v1_alt', 'chr12_KN538369v1_fix', 'chr12_KN538370v1_fix', 'chr12_KQ759760v1_fix', 'chr12_KZ208916v1_fix', 'chr12_KZ208917v1_fix', 'chr12_KZ559112v1_alt', 'chr12_ML143361v1_fix', 'chr12_ML143362v1_fix', 'chr13_KI270838v1_alt', 'chr13_KI270840v1_alt', 'chr13_KI270842v1_alt', 'chr13_KN538371v1_fix', 'chr13_KN538372v1_fix', 'chr13_KN538373v1_fix', 'chr13_ML143365v1_fix', 'chr13_ML143366v1_fix', 'chr14_GL000009v2_random', 'chr14_GL000194v1_random', 'chr14_GL000225v1_random', 'chr14_KI270726v1_random', 'chr14_KI270844v1_alt', 'chr14_KI270845v1_alt', 'chr14_KI270846v1_alt', 'chr14_KI270847v1_alt', 'chr14_KZ208919v1_alt', 'chr14_KZ208920v1_fix', 'chr14_ML143367v1_fix', 'chr15_GL383554v1_alt', 'chr15_GL383555v2_alt', 'chr15_KI270727v1_random', 'chr15_KI270848v1_alt', 'chr15_KI270849v1_alt', 'chr15_KI270850v1_alt', 'chr15_KI270851v1_alt', 'chr15_KI270852v1_alt', 'chr15_KI270905v1_alt', 'chr15_KI270906v1_alt', 'chr15_KQ031389v1_alt', 'chr15_ML143369v1_fix', 'chr15_ML143370v1_fix', 'chr15_ML143371v1_fix', 'chr15_ML143372v1_fix', 'chr16_GL383556v1_alt', 'chr16_GL383557v1_alt', 'chr16_KI270728v1_random', 'chr16_KI270853v1_alt', 'chr16_KI270854v1_alt', 'chr16_KI270855v1_alt', 'chr16_KI270856v1_alt', 'chr16_KQ090026v1_alt', 'chr16_KQ090027v1_alt', 'chr16_KV880768v1_fix', 'chr16_KZ208921v1_alt', 'chr16_KZ559113v1_fix', 'chr16_ML143373v1_fix', 'chr17_GL000205v2_random', 'chr17_GL000258v2_alt', 'chr17_GL383563v3_alt', 'chr17_GL383564v2_alt', 'chr17_GL383565v1_alt', 'chr17_GL383566v1_alt', 'chr17_JH159146v1_alt', 'chr17_JH159147v1_alt', 'chr17_JH159148v1_alt', 'chr17_KI270857v1_alt', 'chr17_KI270858v1_alt', 'chr17_KI270859v1_alt', 'chr17_KI270860v1_alt', 'chr17_KI270861v1_alt', 'chr17_KI270862v1_alt', 'chr17_KI270907v1_alt', 'chr17_KI270908v1_alt', 'chr17_KI270909v1_alt', 'chr17_KI270910v1_alt', 'chr17_KV575245v1_fix', 'chr17_KV766196v1_fix', 'chr17_KV766198v1_alt', 'chr17_KZ559114v1_alt', 'chr17_ML143374v1_fix', 'chr17_ML143375v1_fix', 'chr18_GL383567v1_alt', 'chr18_GL383569v1_alt', 'chr18_GL383570v1_alt', 'chr18_GL383571v1_alt', 'chr18_GL383572v1_alt', 'chr18_KI270863v1_alt', 'chr18_KI270911v1_alt', 'chr18_KI270912v1_alt', 'chr18_KQ090028v1_fix', 'chr18_KQ458385v1_alt', 'chr18_KZ208922v1_fix', 'chr18_KZ559115v1_fix', 'chr18_KZ559116v1_alt', 'chr19_GL000209v2_alt', 'chr19_GL383573v1_alt', 'chr19_GL383574v1_alt', 'chr19_GL383575v2_alt', 'chr19_GL383576v1_alt', 'chr19_GL949746v1_alt', 'chr19_GL949747v2_alt', 'chr19_GL949748v2_alt', 'chr19_GL949749v2_alt', 'chr19_GL949750v2_alt', 'chr19_GL949751v2_alt', 'chr19_GL949752v1_alt', 'chr19_GL949753v2_alt', 'chr19_KI270865v1_alt', 'chr19_KI270866v1_alt', 'chr19_KI270867v1_alt', 'chr19_KI270868v1_alt', 'chr19_KI270882v1_alt', 'chr19_KI270883v1_alt', 'chr19_KI270884v1_alt', 'chr19_KI270885v1_alt', 'chr19_KI270886v1_alt', 'chr19_KI270887v1_alt', 'chr19_KI270888v1_alt', 'chr19_KI270889v1_alt', 'chr19_KI270890v1_alt', 'chr19_KI270891v1_alt', 'chr19_KI270914v1_alt', 'chr19_KI270915v1_alt', 'chr19_KI270916v1_alt', 'chr19_KI270917v1_alt', 'chr19_KI270918v1_alt', 'chr19_KI270919v1_alt', 'chr19_KI270920v1_alt', 'chr19_KI270921v1_alt', 'chr19_KI270922v1_alt', 'chr19_KI270923v1_alt', 'chr19_KI270929v1_alt', 'chr19_KI270930v1_alt', 'chr19_KI270931v1_alt', 'chr19_KI270932v1_alt', 'chr19_KI270933v1_alt', 'chr19_KI270938v1_alt', 'chr19_KN196484v1_fix', 'chr19_KQ458386v1_fix', 'chr19_KV575246v1_alt', 'chr19_KV575247v1_alt', 'chr19_KV575248v1_alt', 'chr19_KV575249v1_alt', 'chr19_KV575250v1_alt', 'chr19_KV575251v1_alt', 'chr19_KV575252v1_alt', 'chr19_KV575253v1_alt', 'chr19_KV575254v1_alt', 'chr19_KV575255v1_alt', 'chr19_KV575256v1_alt', 'chr19_KV575257v1_alt', 'chr19_KV575258v1_alt', 'chr19_KV575259v1_alt', 'chr19_KV575260v1_alt', 'chr19_ML143376v1_fix', 'chr1_GL383518v1_alt', 'chr1_GL383519v1_alt', 'chr1_GL383520v2_alt', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270759v1_alt', 'chr1_KI270760v1_alt', 'chr1_KI270761v1_alt', 'chr1_KI270762v1_alt', 'chr1_KI270763v1_alt', 'chr1_KI270765v1_alt', 'chr1_KI270766v1_alt', 'chr1_KI270892v1_alt', 'chr1_KN196472v1_fix', 'chr1_KN196473v1_fix', 'chr1_KN196474v1_fix', 'chr1_KN538360v1_fix', 'chr1_KN538361v1_fix', 'chr1_KQ031383v1_fix', 'chr1_KQ458382v1_alt', 'chr1_KQ458383v1_alt', 'chr1_KQ458384v1_alt', 'chr1_KQ983255v1_alt', 'chr1_KV880763v1_alt', 'chr1_KZ208904v1_alt', 'chr1_KZ208905v1_alt', 'chr1_KZ208906v1_fix', 'chr20_GL383577v2_alt', 'chr20_KI270869v1_alt', 'chr20_KI270870v1_alt', 'chr20_KI270871v1_alt', 'chr21_GL383579v2_alt', 'chr21_GL383580v2_alt', 'chr21_GL383581v2_alt', 'chr21_KI270872v1_alt', 'chr21_KI270873v1_alt', 'chr21_KI270874v1_alt', 'chr21_ML143377v1_fix', 'chr22_GL383582v2_alt', 'chr22_GL383583v2_alt', 'chr22_KB663609v1_alt', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270875v1_alt', 'chr22_KI270876v1_alt', 'chr22_KI270877v1_alt', 'chr22_KI270878v1_alt', 'chr22_KI270879v1_alt', 'chr22_KI270928v1_alt', 'chr22_KN196485v1_alt', 'chr22_KN196486v1_alt', 'chr22_KQ458387v1_alt', 'chr22_KQ458388v1_alt', 'chr22_KQ759761v1_alt', 'chr22_KQ759762v1_fix', 'chr22_ML143378v1_fix', 'chr22_ML143380v1_fix', 'chr2_GL383521v1_alt', 'chr2_GL383522v1_alt', 'chr2_GL582966v2_alt', 'chr2_KI270767v1_alt', 'chr2_KI270768v1_alt', 'chr2_KI270769v1_alt', 'chr2_KI270770v1_alt', 'chr2_KI270772v1_alt', 'chr2_KI270773v1_alt', 'chr2_KI270774v1_alt', 'chr2_KI270776v1_alt', 'chr2_KI270893v1_alt', 'chr2_KI270894v1_alt', 'chr2_KN538362v1_fix', 'chr2_KN538363v1_fix', 'chr2_KQ031384v1_fix', 'chr2_KQ983256v1_alt', 'chr2_KZ208907v1_alt', 'chr2_KZ208908v1_alt', 'chr2_ML143341v1_fix', 'chr2_ML143342v1_fix', 'chr3_GL383526v1_alt', 'chr3_JH636055v2_alt', 'chr3_KI270777v1_alt', 'chr3_KI270779v1_alt', 'chr3_KI270780v1_alt', 'chr3_KI270781v1_alt', 'chr3_KI270782v1_alt', 'chr3_KI270783v1_alt', 'chr3_KI270784v1_alt', 'chr3_KI270895v1_alt', 'chr3_KI270924v1_alt', 'chr3_KI270934v1_alt', 'chr3_KI270935v1_alt', 'chr3_KI270936v1_alt', 'chr3_KI270937v1_alt', 'chr3_KN196475v1_fix', 'chr3_KN196476v1_fix', 'chr3_KN538364v1_fix', 'chr3_KQ031385v1_fix', 'chr3_KV766192v1_fix', 'chr3_KZ208909v1_alt', 'chr3_KZ559102v1_alt', 'chr3_KZ559103v1_alt', 'chr3_KZ559104v1_fix', 'chr3_KZ559105v1_alt', 'chr3_ML143343v1_alt', 'chr4_GL000257v2_alt', 'chr4_GL383527v1_alt', 'chr4_GL383528v1_alt', 'chr4_KI270785v1_alt', 'chr4_KI270786v1_alt', 'chr4_KI270788v1_alt', 'chr4_KI270789v1_alt', 'chr4_KI270790v1_alt', 'chr4_KI270896v1_alt', 'chr4_KI270925v1_alt', 'chr4_KQ090014v1_alt', 'chr4_KQ090015v1_alt', 'chr4_KQ983257v1_fix', 'chr4_KQ983258v1_alt', 'chr4_KV766193v1_alt', 'chr4_ML143344v1_fix', 'chr4_ML143345v1_fix', 'chr4_ML143347v1_fix', 'chr4_ML143349v1_fix', 'chr5_GL339449v2_alt', 'chr5_GL383531v1_alt', 'chr5_GL383532v1_alt', 'chr5_GL949742v1_alt', 'chr5_KI270791v1_alt', 'chr5_KI270792v1_alt', 'chr5_KI270793v1_alt', 'chr5_KI270794v1_alt', 'chr5_KI270795v1_alt', 'chr5_KI270796v1_alt', 'chr5_KI270897v1_alt', 'chr5_KI270898v1_alt', 'chr5_KN196477v1_alt', 'chr5_KV575243v1_alt', 'chr5_KV575244v1_fix', 'chr5_ML143350v1_fix', 'chr6_GL000250v2_alt', 'chr6_GL000251v2_alt', 'chr6_GL000252v2_alt', 'chr6_GL000253v2_alt', 'chr6_GL000254v2_alt', 'chr6_GL000255v2_alt', 'chr6_GL000256v2_alt', 'chr6_GL383533v1_alt', 'chr6_KB021644v2_alt', 'chr6_KI270758v1_alt', 'chr6_KI270797v1_alt', 'chr6_KI270798v1_alt', 'chr6_KI270799v1_alt', 'chr6_KI270800v1_alt', 'chr6_KI270801v1_alt', 'chr6_KI270802v1_alt', 'chr6_KN196478v1_fix', 'chr6_KQ031387v1_fix', 'chr6_KQ090016v1_fix', 'chr6_KV766194v1_fix', 'chr6_KZ208911v1_fix', 'chr7_GL383534v2_alt', 'chr7_KI270803v1_alt', 'chr7_KI270804v1_alt', 'chr7_KI270805v1_alt', 'chr7_KI270806v1_alt', 'chr7_KI270807v1_alt', 'chr7_KI270808v1_alt', 'chr7_KI270809v1_alt', 'chr7_KI270899v1_alt', 'chr7_KQ031388v1_fix', 'chr7_KV880764v1_fix', 'chr7_KV880765v1_fix', 'chr7_KZ208912v1_fix', 'chr7_KZ208913v1_alt', 'chr7_KZ559106v1_alt', 'chr7_ML143352v1_fix', 'chr8_KI270810v1_alt', 'chr8_KI270811v1_alt', 'chr8_KI270812v1_alt', 'chr8_KI270813v1_alt', 'chr8_KI270814v1_alt', 'chr8_KI270815v1_alt', 'chr8_KI270816v1_alt', 'chr8_KI270817v1_alt', 'chr8_KI270818v1_alt', 'chr8_KI270819v1_alt', 'chr8_KI270821v1_alt', 'chr8_KI270822v1_alt', 'chr8_KI270900v1_alt', 'chr8_KI270901v1_alt', 'chr8_KI270926v1_alt', 'chr8_KV880766v1_fix', 'chr8_KZ208914v1_fix', 'chr8_KZ208915v1_fix', 'chr8_KZ559107v1_alt', 'chr9_GL383539v1_alt', 'chr9_GL383540v1_alt', 'chr9_GL383541v1_alt', 'chr9_GL383542v1_alt', 'chr9_KI270823v1_alt', 'chr9_KN196479v1_fix', 'chr9_KQ090018v1_alt', 'chrUn_GL000195v1', 'chrUn_GL000213v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_KI270442v1', 'chrUn_KI270744v1', 'chrUn_KI270750v1', 'chrX_KI270880v1_alt', 'chrX_KI270881v1_alt', 'chrX_KI270913v1_alt', 'chrX_ML143381v1_fix', 'chrY_KN196487v1_fix', 'chrY_KZ208924v1_fix'])
validchrom = set(mainchrom).union(chm13chrom).union(alterchrom)


def runblastn(queryfile, reffile, alignfile):
	
	cmd = "blastn -query {} -task megablast -db {} -gapopen 10 -gapextend 2  -outfmt 17 -out {}  -num_threads 1 -max_target_seqs 100".format( queryfile, reffile , alignfile)
	
	if os.path.isfile(alignfile) and os.stat(alignfile).st_size > 10 :
		
		pass
		return 1
	else:
		pass
		#os.system(cmd)
		return 0
	
def makereverse(seq):
	
	tran=str.maketrans('ATCGatcg', 'TAGCtagc')
	
	return seq[::-1].translate(tran)

def cigarextend(cigars, refseq, queryseq):
	
	cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
	
	rposi = 0
	qposi = 0
	newcigars = []
	
	newref = ""
	for cigar in cigars:
		
		
		thesize = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == "=" or thetype == "M":
			
			
			matchormismatch = 0
			last_i = 0
			newcigar = []
			for i, (r,q) in enumerate(zip(refseq[rposi:(rposi+thesize)], queryseq[qposi:(qposi+thesize)])):
				
				if matchormismatch ==0 and r.upper() != q.upper():
					
					newcigar.append("{}=".format(i-last_i))
					last_i = i
					matchormismatch = 1
					check = 1
					
				elif matchormismatch ==1 and r.upper() == q.upper():
					
					#newcigar.append("{}X{}".format(i-last_i, queryseq[(qposi+last_i):(qposi+i)]))
					newcigar.append("{}X".format(i-last_i))
					last_i = i
					matchormismatch = 0
					
					
			i = thesize
			if matchormismatch ==0:
				newcigar.append("{}=".format(i-last_i))
			else:
				#newcigar.append("{}X{}".format(i-last_i,queryseq[(qposi+last_i):(qposi+i)]))
				newcigar.append("{}X".format(i-last_i))
				
			cigar = "".join(newcigar)
			
			newcigars.append(cigar)
			
			qposi += thesize
			
			#newref += refseq[refposi:(refposi+thesize)]
			rposi += thesize
			
			
		elif thetype == "X" or thetype == "S":
			
			#newcigars.append(cigar+queryseq[qposi:(qposi+thesize)])
			newcigars.append(cigar)
			qposi += thesize
			
			#newref += "X"*thesize
			rposi += thesize
			
		elif thetype == "I" or thetype == "H":
			
			newcigars.append(cigar)
			
			rposi += thesize
			
		elif thetype == "D":
			
			#newcigars.append(cigar+queryseq[qposi:(qposi+thesize)])
			newcigars.append(cigar)
			qposi += thesize
			#newref += "-"*thesize
			
			
	return "".join(newcigars)


def cigar_findqrange(cigars, qstart , qend):
	
	cigars = re.findall(r'\d+[a-zA-Z=]', cigars)
	
	
	rposi = 0
	qposi = 0
	newcigars = []
	refstart = -1
	refend = -1
	newref = ""
	for cigar in cigars:
		
		lastposi = qposi
		
		thesize = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype == "=" or thetype == "M":
			
			qposi += thesize
			rposi += thesize
			
			
		elif thetype == "X" or thetype == "S":
			
			qposi += thesize
			rposi += thesize
			
		elif thetype == "I" or thetype == "H":
			
			rposi += thesize
			
		elif thetype == "D":
			
			qposi += thesize
			newref += "-"*thesize
			
		if qposi <= qstart:
			
			continue
		
		if lastposi <= qstart:
			
			if thetype in ['M','=', 'X','S', 'I','H']:
				
				refstart = rposi - thesize + qstart - lastposi
			else:
				refstart = rposi
				
			corrsize = min(qend,qposi) - qstart
			newcigars.append("{}{}".format(corrsize,thetype))
			
			
		if qposi >= qend:
			
			if thetype in ['M','=', 'X','S', 'I','H']:
				
				refend = rposi - qposi + qend
			else:
				
				refend = rposi
				
			if lastposi > qstart:
				corrsize = qend - max(qstart, lastposi)
				
				newcigars.append("{}{}".format(corrsize,thetype))
				
			break
		
		if lastposi > qstart:
			
			newcigars.append(cigar)
			
			
	return  [refstart, refend, "".join(newcigars)]

def cigartoscore(cigars):
	
	allcigars = re.findall(r'\d+[A-Z=]',cigars)
	
	score = 0

	rposi = 0
	qposi = 0
	
	match = 0
	mismatch = 0
	sv = 0
	for cigar in allcigars:
		
		size = int(cigar[:-1])
		thetype = cigar[-1]
		
		if thetype in ['=','M']:
			
			#usize = sum(1 for c in rseq[rposi:(rposi+size)] if c.isupper())

			score += size
			rposi += size
			qposi += size


		elif thetype in ['X']:
			#usize = sum(1 for c in rseq[rposi:(rposi+size)] if c.isupper())

			score -= 20*size
			rposi += size
			qposi += size
 
		elif thetype in ['D']:
			#usize = sum(1 for c in rseq[rposi:(rposi+size)] if c.isupper())
			score -= (size + 20)
			rposi += size
			
			sv -=  (size+ 20)
		elif thetype in ['I']:
			#usize = sum(1 for c in qseq[qposi:(qposi+size)] if c.isupper())
			score -= (size//10 + 20)
			sv -= 20
			qposi += size
			
	#print("score",match, mismatch,sv,score)
	return score



def getcigar(alignpath, valid_aligns, queryseq ,refseq_f, refseq_r):
	
	pathes = []
	for path, aligns in zip(alignpath, valid_aligns):
		
		index, qstart, qend = path
		
		elements = aligns[-1]
		
		rname, strand, qname, alignment_qstart, cigar = elements[0], elements[1], elements[2], int(elements[3]) -1, elements[5]
		
		newcigar = cigar_findqrange(cigar, qstart - alignment_qstart , qend  - alignment_qstart ) 
		
		if strand == '+':
			rseq = refseq_f[newcigar[0]:]
		else:
			rseq = refseq_r[newcigar[0]:]
			
		qseq = queryseq[qstart:]
		
		newcigar[-1] = cigarextend(newcigar[-1], rseq, qseq )
		
		pathes.append(path + [rname, strand] + newcigar )
		
	fullpath = ""
	fullcigar = []
	reftemplate = []
	last_rname = ""
	last_qend = 0
	last_rend = 0
	last_strand = ""
	last_rclip = ""
	
	qstarts = []
	pathrange = []
	for path in pathes:
		
		index, qstart, qend, rname, strand, rstart, rend, cigar = path
		
		qstarts.append("{}_{}".format(qstart,qend))
		cigar = cigar.replace("H","I")
		
		if last_rname == rname and last_strand == strand and rstart >= last_rend:
			
			fullcigar[-1] = ""
			if rstart - last_rend:
				fullcigar.append(str(rstart - last_rend) + "I")
				
				
			if qstart - last_qend:
				#fullcigar.append(str(qstart - last_qend) + "D" + queryseq[last_qend:qstart])
				fullcigar.append(str(qstart - last_qend) + "D")
					
			fullcigar.extend(re.findall(r'\d+[a-zA-Z=][actgACTG]*', cigar))
			
			pathrange[-1][1] = rend
			
			
		else:
			if strand == "+":
				fullcigar.append(">")
			else:
				fullcigar.append("<")
				
			if rstart:
				fullcigar.append(str(rstart) + "H")
				pass
				
			if qstart - last_qend:
				#fullcigar.append(str(qstart - last_qend) + "D" + queryseq[last_qend:qstart])
				fullcigar.append(str(qstart - last_qend) + "D")
				
			fullcigar.extend(re.findall(r'\d+[a-zA-Z=][actgACTG]*', cigar))
			
			if strand == "+":
				fullpath += ">{}".format(rname)
			else:
				fullpath += "<{}".format(rname)
				
			pathrange.append([rstart, rend])
			
		if len(refseq_f) - rend:
			fullcigar.append(str(len(refseq_f) - rend) + "H")
			
		last_qend = qend
		last_rend = rend
		last_rname = rname
		last_strand = strand
		
	qstarts.append(len(queryseq))
	fullcigar = "".join(fullcigar).replace("D",'@').replace('I','D').replace('@','I')
	
	return fullcigar,pathes




def localization(self, outfile):
	
	self.localized = {}
	alllocations = []
	
	for qchr, index in self.allpairs.items():
		
		locations = []
		for rchr, subindex in index.items():
			
			locations = getmaplocation(self.data.iloc[subindex].reset_index(drop=True))
			
			
			
		alllocations.extend(locations)
		self.localized[qchr] = max(locations, key=lambda x:(x[0], x[1].iloc[0]["ref"] if len(x[1]) else float('-inf') ))[0]
		
		
	mapped = []
	unmapped = []
	location_count=cl.defaultdict(int)
	
	mapped_list = set()
	
	mapped_score = cl.defaultdict(int)
	for score,data in alllocations :
		
		if len(data) == 0:
			continue
		
		allqstart = min(data["qstart"])
		allqend = max(data["qend"])
		
		allrstart = min(data["rstart"])
		allrend = max(data["rend"])
		
		qlength,rlength = data.iloc[0]["qlength"],data.iloc[0]["rlength"]
		
		query,ref = data.iloc[0]["query"],data.iloc[0]["ref"]
		
		location_count[query]+=1
		
		totalscores = sum(data["match"])
		quality = max(data["quality"])
		
		strand_size= {"+":0,'-':0}
		
		for strand, start, end in data[["strand","rstart","rend"]].values.tolist():
			
			strand_size[strand] += abs(start - end)
			
		if strand_size["+"] >= strand_size["-"]:
			strand = "+"
		else:
			strand = "-"
			
			
			
		if totalscores < self.coverage * (qlength - self.Ngaps[query]):
			
			unmapped.append([query, qlength, allqstart , allqend, strand, ref, rlength,allrstart, allrend, totalscores, quality])
			
		else:
			if mapped_score[query] < totalscores:
				mapped.append([query, qlength, allqstart , allqend, strand, ref, rlength,allrstart, allrend, totalscores, quality])
				mapped_score[query] = totalscores
				mapped_list.add(query)  
				
				
				
	unmapped = [x for x in unmapped if x[0] not in mapped_list]
	
	pd.DataFrame.from_records(mapped).to_csv(outfile,mode='a',header=None,index=False, sep = "\t")
	
	pd.DataFrame.from_records(unmapped).to_csv(outfile+"_unmapped.txt",mode='a',header=None,index=False, sep = "\t")
	
	return self


def resolveoverlap(allaligns):
	
	all_coordinates = [x for y in allaligns for x in y[:2]]
	scores = [y[4] for y in allaligns]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	current_alignscores = []
	allgroups = [[0,0,0,0]]
	
	ended_indexes = []
	last_score = -1
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		score = scores[index//2]
		
		if index %2 == 0 :
			
			last_score = current_alignscores[-1][0] if len(current_alignscores) else 0
			if len(current_alignscores) ==0 or score > last_score :
				
				if len(current_alignscores):
					allgroups[-1][-1] = coordinate
					
				allgroups.append([score, index//2, coordinate, coordinate])
				
			bisect.insort(current_alignscores, [score, -(index//2)])
			
		else:
			
			allgroups[-1][-1] = coordinate
			
			highest_index = -current_alignscores[-1][1]
			
			current_alignscores.remove([scores[index//2], -(index//2)])
			
			if  index//2 == highest_index  and len(current_alignscores):
				
				allgroups.append([current_alignscores[-1][0], -current_alignscores[-1][1]]+[coordinate, coordinate])
				
	allgroups = [[x[1]]+x[2:] for x in allgroups[1:]]
	
	last_index = -1
	
	new_allgroups = []
	for segment in allgroups:
		
		index,start, end = segment[0], segment[1], segment[2]
		
		if index == last_index :
			new_allgroups[-1][-1] = end
		else:
			new_allgroups.append(segment)
			
		last_index = index 
		last_end = end
		
	allgroups = [x for x in new_allgroups if x[2] - x[1] > 50]
	
	return allgroups

def getmaplocation(allaligns, allowgap = 50000):
	
	starts = [x[2] for x in allaligns]
	ends = [x[3] for x in allaligns]
	
	l = len(starts)
	allposi = starts+ends
	
	allposi_sortindex = sorted(range(len(starts+ends)), key = lambda i: allposi[i])
	
	allgroups = []
	current_group = [-allowgap - 1,-allowgap - 1,0,[]]
	last_endpoint = -allowgap - 1
	last_startpoint = -allowgap - 1
	overlap_count = 0
	
	for rank, index in enumerate( allposi_sortindex):
		
		if index >= l :
			
			overlap_count -= 1
			
			if overlap_count == 0:
				
				last_endpoint = allposi[index]
				current_group[2] += last_endpoint - last_startpoint
				
				
		else:
			
			overlap_count += 1
			
			if overlap_count == 1:
				last_startpoint = allposi[index]
				
			if overlap_count == 1 and allposi[index] - last_endpoint > allowgap:
				
				current_group[1] = last_endpoint
				allgroups.append(current_group)
				current_group = [allposi[index] , allposi[index] , 0 , []]
				
			current_group[-1].append(index)
			
	current_group[1] = last_endpoint
	allgroups.append(current_group)
	
	return allgroups[1:]


def readfasta(reffiles, userefnames = set([])):
	
	refseqs = manager.dict()
	refnames = []
	for reffile in reffiles.split(","):
		
		ifname = 0
		refname = ""
		refseq = ""
		with open(reffile, mode = 'r') as f:
			
			for line in f:
				
				if len(line) == 0:
					continue
				
				if line[0] == '>':
					
					refseqs[refname] = refseq
					
					refseq = ""
					refname = line.split()[0][1:]
					refnames.append(refname)
					
					if len(userefnames) ==0 or refname in userefnames:
						ifname = 1
					else:
						ifname = 0
						
				elif ifname:
					refseq += line.strip()
					
			refseqs[refname] = refseq
			
	for refname in list(refseqs.keys()):
		refseqs[refname + "_r"] = makereverse(refseqs[refname])
		
		
	return refseqs, refnames

def approlocation(alignfile, queryseqs):
	
	lindex = 0
	allrefnames = set([])
	allaligns = cl.defaultdict(list)
	refalignsizes = cl.defaultdict(int)
	with open(alignfile, mode = 'r') as rfile:
		
		for line in rfile:
			
			if len(line.strip()) == 0 or line[0]=="@":
				continue
			
			
			lindex += 1
			
			elements = line.strip().split()
			#name, qsize, qstart, qend, path, rsize, rstart, rend = elements[:8]
			identity =float( [x for i,x in enumerate(elements) if x[:5] == "PI:f:"][0][5:])
			match = int( [x for i,x in enumerate(elements) if x[:5] == "AS:i:"][0][5:] ) 
			if identity < 95 or elements[4] != '255':
				
				continue
			
			rname, strand, qname, qstart, cigar = elements[0], elements[1], elements[2], int(elements[3])-1, elements[5]
			
			qname = qname.split("_group")
			if "NC_0609" in qname[1]:
				qname[1] = qname[1].replace("v",'.')
			qname = qname[0].replace("_",'-').replace("v",".") + "_group" +qname[1]
			
			
			if qname not in queryseqs:
				qname = [x for x in queryseqs.keys() if x.replace("_",'-').replace("v",".") == qname.replace("_",'-').replace("v",".")][0]
				
			strand = "+" if strand == '0' else '-'
			elements[1] = strand
			
			# or qname not in [ "SMN_group1_HG02723_h1_198","SMN_group1_chr5_KI270897v1_alt_254","SMN_group1_HG02148_h2_154"] 
			if (rname not in validchrom ):
				continue
			
			allrefnames.add(rname)
			#qname = qnames[ int(qname.split('_')[1]) - 1 ]
			#elements[2] = qname
			
			
			rsize = sum([int(x[:-1]) for x in re.findall(r'\d+[MXI]',cigar)])
			qsize = sum([int(x[:-1]) for x in re.findall(r'\d+[SMXD]',cigar)])
			
			
			if strand == '+':
				leftclip = re.findall(r'^\d+H',cigar)
			else:
				leftclip = re.findall(r'\d+H$',cigar)
				
			if len(leftclip) == 0:
				leftclip = 0
			else:
				leftclip = int(leftclip[0][:-1])
				
			alignment = [qstart, qstart + qsize, leftclip, leftclip + rsize]
			
			refalignsizes[(qname,rname)] += rsize
			
			allaligns[(qname,rname)].append( alignment + [ match, lindex] )
			
			
	Ngaps = {}
	
	allalignedlocations = cl.defaultdict(list)
	
	loadindex = set([])
	loadrefnames = set([])
	for (qname, rname) , aligns in allaligns.items():
		
		if qname not in Ngaps:
			nsize = queryseqs[qname].count('N') + queryseqs[qname].count('n')
			Ngaps[qname] = nsize
			
		nsize = Ngaps[qname]
		
		qsize = len(queryseqs[qname]) - nsize
		
		if refalignsizes[(qname,rname)] < qsize//5:
			continue
		
		#locations = getmaplocation(aligns, len(queryseqs[qname]))
		locations = getmaplocation(aligns, len(queryseqs[qname])/2)
		
		
		locations = [x[3] for x in locations if x[2] > qsize//5]
		
		if len(locations):
			
			loadrefnames.add(rname)
			
		for location in locations:
			
			local_aligns = [aligns[i] for i in location]
			
			alignpath = resolveoverlap(local_aligns)
			
			valid_aligns = [local_aligns[x[0]] for x in alignpath]
			
			loadindex.update([x[-1] for x in valid_aligns])
			
			allalignedlocations[(qname, rname)].append([alignpath, valid_aligns])
			
	return allalignedlocations, loadindex, loadrefnames,Ngaps

def getaligndetails(qname, rname, nsize, alignpath, valid_aligns, queryseqs, refseqs):
	
	
	queryseq, refseq, refseq_r = queryseqs[qname], refseqs[rname], refseqs[rname+"_r"]
	
	
	loca_s = min([x[2] for x in valid_aligns])
	loca_e = max([x[3] for x in valid_aligns])
	
	fullcigar, pathes = getcigar(alignpath, valid_aligns, queryseq, refseq, refseq_r)
	
	score = cigartoscore(fullcigar) + 11*nsize//10
	
	#index, qstart, qend, rname, strand, rstart, rend, cigar
	#strd_score = sum([ ( path[2]-path[1] ) * (1 if path[4] =='+' else -1)  for path in pathes])

	cigar2 = fullcigar.replace("<","><")
	cigar2  = "".join([x for x in  cigar2.split(">")[1:] if x[0] != "<"])
	size_f  = sum([int(x[:-1]) for x in re.findall(r'\d+[=M]',cigar2)]) 
	size = sum([int(x[:-1]) for x in re.findall(r'\d+[=M]',fullcigar)])
	strd_score = 2*size_f - size

	
	return [qname, len( queryseq) - nsize, rname, loca_s, loca_e, [ '-' , '+'][int(strd_score >= 0)], score, fullcigar]

def mapgenes(alignfile, queryfile, reffiles, threads = 4):
	
	queryseqs, qnames =  readfasta(queryfile)
	qnames = [x for x in qnames]
	
	allalignedlocations, loadindex, loadrefnames, Ngaps = approlocation(alignfile, queryseqs)
	
	
	refseqs, rnames = readfasta(reffiles, loadrefnames)
	
	loadelements = {}
	lindex = 0
	with open(alignfile, mode = 'r') as rfile:
		
		for line in rfile:
			
			if len(line.strip()) == 0 or line[0]=="@":
				continue
			
			lindex += 1
			
			if lindex in loadindex:
				
				elements = line.strip().split()
				rname, strand, qname, qstart, cigar = elements[0], elements[1], elements[2], int(elements[3])-1, elements[5]
				strand = "+" if strand == '0' else '-'
				elements[1] = strand
				
				qname = qname.split("_group")
				if "NC_0609" in qname[1]:
					qname[1] = qname[1].replace("v", ".")
				qname = qname[0].replace("_",'-').replace("v",".") + "_group" +qname[1]
				
				if qname not in queryseqs:
					
					qname = [x for x in queryseqs.keys() if x.replace("_",'-').replace("v",".") == qname.replace("_",'-').replace("v",".")][0]
					
					
				elements[2] = qname
				
				
				
				loadelements[lindex] = elements
				
	results = []
	
	p=mul.Pool(processes=threads)
	
	for  (qname, rname), aligndata in allalignedlocations.items():
		
		for [alignpath, valid_aligns] in aligndata:
			
			valid_aligns = [x[:-1]+[loadelements[x[-1]]] for x in valid_aligns]
			
			result = p.apply_async(getaligndetails,(qname, rname, Ngaps[qname], alignpath, valid_aligns, queryseqs, refseqs))
			results.append(result)
			
	p.close()
	p.join()
	
	results = [x.get() for x in results]
	
	return results, qnames

def runmapgene(alignfile, queryfile, reffile, threads):
	
	iffile = runblastn(queryfile, reffile + "_db", alignfile)
	
	if not iffile:
		return {},{},{},{},{}
	
	maplocations, qnames = mapgenes(alignfile, queryfile, reffile, threads)
	
	mapped = cl.defaultdict(list)
	unmapped = cl.defaultdict(list)
	multimapped = cl.defaultdict(list)
	nomapped = [x for x in qnames]
	
	maplocations_query = cl.defaultdict(list)
	
	for  maplocation in maplocations:
		
		qname = maplocation[0]
		
		maplocations_query[qname].append(maplocation)
		
	for  qname in qnames:
		
		maplocations = maplocations_query[qname]
		
		maplocations = sorted(maplocations, key = lambda x:x[6], reverse = 1)
		
		maplocations_query[qname] = maplocations
		
		maplocations = [x for x in maplocations if x[6] > 0]
		
		contig = "_".join(qname.split('_group')[1].split("_")[1:-1])
		
		if contig in alterchrom:
			maplocations = [x for x in maplocations if x[2] not in alterchrom]
			
		if len(maplocations) == 0:
			continue
		
		nomapped.remove(qname)
		
		size = maplocations[0][1]
		
		highestscore = maplocations[0][6]
		
		valide_locations = [x for x in maplocations if x[6] > highestscore * 0.95]
		
		if highestscore < size * 0.5:
			
			unmapped[qname] = valide_locations
			continue
		
		
		hg38_loc = [x for x in valide_locations if x[2] in mainchrom and x[2] != 'chrY']
		
		chm13_loc = [x for x in valide_locations if x[2] in chm13chrom and x[2] != 'NC_060948.1']
		
		if (len(hg38_loc) > 1) or (len(chm13_loc) > 1):
			
			multimapped[qname] = valide_locations
			continue
		
		
		mapped[qname] = [maplocations[0]]
		
	return maplocations_query, mapped, unmapped, multimapped, nomapped


def main(args):
	
	allalignedlocations, mapped, unmapped,multimapped,nomapped = runmapgene(args.input, args.query, args.ref, args.threads)
	
	if len(allalignedlocations) == 0:
		return
	
	with open(args.output, mode = 'w') as f:
		
		for qname,aligns in allalignedlocations.items():
			
			for align in aligns:
				f.write("\t".join(list(map(str, align))) + "\n")
				
	with open(args.output+"_unmapped.txt", mode = 'w') as f:
		for qname in nomapped:
			f.write(qname + "\n")   
			
			
	with open(args.output+"_mapped.txt", mode = 'w') as f:
		
		for qname,aligns in mapped.items():
			
			for align in aligns:
				f.write("\t".join(list(map(str, align))) + "\n")
				
	with open(args.output+"_partimapped.txt", mode = 'w') as f:
		
		for qname,aligns in unmapped.items():
			
			for align in aligns:
				f.write("\t".join(list(map(str, align))) + "\n")
				
	with open(args.output+"_multimapped.txt", mode = 'w') as f:
		
		for qname,aligns in multimapped.items():
			
			for align in aligns:
				f.write("\t".join(list(map(str, align))) + "\n")
				
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	
	parser.add_argument("-r", "--ref", help="path to output file", dest="ref", type=str,required=True) 
	parser.add_argument("-q", "--query", help="path to output file", dest="query", type=str,required=True)  
	parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
	parser.add_argument("-t", "--threads", help="path to output file", dest="threads", type=int,default=4)
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
