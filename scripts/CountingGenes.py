#!/usr/bin/env python3

import pandas as pd
import collections as cl
import argparse
import pandas as pd


def determinescore(exons, thetype, exonnum):

        factor = 10 if thetype == "protein_coding" else 1

        if len(exons) >= exonnum:
                factor *= 2

        #penalty = 100 * exonnum -  100 * len(exons)

        score = ( sum([x[5] for x in exons]) ) 


        return factor * score 



def resolveoverlap(aligns_bycontig, exonstotran):


        aligns_bycontig_sort = sorted(aligns_bycontig, key = lambda x: x[-1]/x[1],reverse = 1)

        aligns_exons = cl.defaultdict(list)
        used_index = set()
        for index,row in enumerate(aligns_bycontig_sort):

                exonid, size, qstart0, qend0, strd0, rstart0, rend0,match0 = row
                tranid = exonstotran.get(exonid,None)


                if tranid is None or abs(qstart0 - qend0) < min(100, 0.5*size): 
                        continue


                overlap = 0
                ifoverlap = 0   
                for exon in aligns_exons[exonid]:

                        overlap = ( exon[1] - exon[0] ) + ( qend0 - qstart0 ) - max(qend0, exon[1]) + min(qstart0, exon[0])

                        if overlap > 0:
                                ifoverlap = 1
                                break

                if ifoverlap == 0:
                        aligns_exons[exonid].append((qstart0, qend0))
                        used_index.add(index)



        return sorted([x for i,x in enumerate(aligns_bycontig_sort) if i in used_index], key = lambda x: x[2])



def cleanoverlapexons(exons,aligns_exons):

        exons_cleaned = []
        exons_new = []
        for exon in exons:

                exonid, qstart0, qend0, strd0, match0, score0 = exon[:6]

                overlap = 0
                ifoverlap = 0
                for aligns_exon in aligns_exons:

                        overlap = ( aligns_exon[1] - aligns_exon[0] ) + ( qend0 - qstart0 ) - max(qend0, aligns_exon[1]) + min(qstart0, aligns_exon[0])

                        if overlap > 0:
                                ifoverlap = 1
                                break

                if ifoverlap == 0:
                        exons_new.append(exon)
                else:
                        exons_cleaned.append(exon)

        return exons_new, exons_cleaned

def resolveoverlap2(allmrnas, transtype, exons_infor, allowtruncate = dict()):



        allmrnas = sorted(allmrnas, key = lambda x:(x[1]*x[2],x[0],-x[-1][0][0]), reverse = 1)

        allmrnas_new = []
        aligns_exons = []
        used_index = set()
        for index in range(len(allmrnas)):

                row = allmrnas[index]
                if len(row[3]):

                        allmrnas_new.append(row)

                        for exon in row[3]:
                                exonid, qstart0, qend0, strd0, match0, score0,  tableindex = exon
                                aligns_exons.append((qstart0, qend0))

                        for index2, row2 in enumerate(allmrnas[(index+1):]):

                                exon_new , cleaned= cleanoverlapexons(row2[3], aligns_exons)


                                allmrnas[index2+index+1][3] = exon_new

                                allmrnas[index2+index+1][2] = sum([x[4] for x in exon_new])


                                allmrnas[index2+index+1][1] = determinescore(exon_new, transtype[row2[0]] ,len(allowtruncate.get(row2[0], exons_infor[row2[0]])))


                        allmrnas = allmrnas[:(index+1)] +sorted(allmrnas[(index+1):], key = lambda x: (x[1]*x[2],x[0]), reverse = 1)



        return allmrnas_new

"""
def getmrna(qranges):


        forward_orreverse = [1 if x[3] == '+' else -1 for x in qranges]

        if sum(forward_orreverse) < 0:
                qranges = qranges[::-1]

        mrnas = []
        for (index, qstart, qend, strd,match, sortindex) in qranges:

                alllinks = [ (i,index - mrna[-1][0] ,mrna) for i,mrna in enumerate(mrnas) if mrna[-1][0] < index]

                print(alllinks)

                if len(alllinks):
                        thelink = sorted(alllinks, key = lambda x: (x[1], -len(x[2])))[0][0]

                        mrnas[thelink].append([index, qstart, qend,strd, match, sortindex])

                else :
                        mrnas.append([[index, qstart, qend,strd, match, sortindex]])

        return mrnas
"""


def getmrna(qranges):

        forward_orreverse = [1 if x[3] == '+' else -1 for x in qranges]

        if sum(forward_orreverse) < 0:
                qranges = qranges[::-1]

        qranges_cleaned = dict()
        for i, (index, qstart, qend, strd, match, sortindex) in enumerate(qranges):

                ifoverlap = 0
                for j, (index2, qstart2, qend2, strd2,match2, sortindex2) in enumerate(qranges[:i]):

                        if j not in qranges_cleaned:
                                continue

                        if strd == strd2 and (qend - qstart) + (qend2- qstart2) > max(qend2,qend) - min(qstart,qstart2) + 20:

                                qranges_cleaned[j][0].append(index)
                                qranges_cleaned[j][1][index] = qstart
                                qranges_cleaned[j][2][index] = qend
                                qranges_cleaned[j][4][index] = match
                                qranges_cleaned[j][5][index] = sortindex

                                ifoverlap = 1
                                break

                if ifoverlap == 0:
                        qranges_cleaned[i] = ([index], {index:qstart}, {index:qend}, strd,{index:match},{index:sortindex})

        qranges_cleaned_ = sorted(list(qranges_cleaned.keys()))
        qranges_cleaned = [qranges_cleaned[x] for x in qranges_cleaned_]

        mrna_overlap = []   
        newmrnas = []
        mrnas = []
        for (index, qstart, qend, strd,match, sortindex) in qranges_cleaned:

                newmrnas = []
                ifadd = 0
                for i,mrna in enumerate(mrnas[::-1]):

                        thisadd = 0
                        #newmrna, cleaned = cleanoverlapexons(mrna , [[qstart, qend,strd, match, sortindex]])

                        lastnode =  mrna[-1][0]


                        if lastnode + 1 in index:

                                mrna.append([lastnode + 1, qstart[lastnode + 1], qend[lastnode + 1],strd, match[lastnode + 1], sortindex[lastnode + 1]])
                                ifadd = 1
                                thisadd =1

                        elif lastnode in index and mrna[-1][-2] < match[lastnode]:

                                #mrna[-1] = [lastnode , qstart[lastnode], qend[lastnode],strd, match[lastnode], sortindex[lastnode]]
                                newmrnas.append([x for x in mrna[:-1]]+[[lastnode , qstart[lastnode], qend[lastnode],strd, match[lastnode], sortindex[lastnode]]])

                        elif max(index) > lastnode + 1 and min([x-lastnode for x in index if x > lastnode]) < 10:

                                indexuse = min([x for x in index if x > lastnode])
                                mrna.append([indexuse, qstart[indexuse], qend[indexuse],strd, match[indexuse], sortindex[indexuse]])

                                ifadd = 1
                                thisadd  = 1

                        else:
                                allindex = set([x[0] for x in mrna])
                                missindex = set([x for x in range(mrna[0][0], mrna[-1][0]) if x not in allindex])

                                indexuses = set(index).intersection(missindex)

                                if len(indexuses):

                                        indexuse = min(indexuses)
                                        newmrna = [x for x in mrna  if x[0] < indexuse]+[[indexuse, qstart[indexuse], qend[indexuse],strd, match[indexuse], sortindex[indexuse]]]

                                        newmrnas.append(newmrna)
                                        ifadd = 1
                                        thisadd = 1

                if not ifadd:

                        newmrnas.append([[min(index), qstart[min(index)], qend[min(index)],strd, match[min(index)], sortindex[min(index)]]])

                mrnas += newmrnas

                exclude = [0 for x in mrnas]
                for i,mrna in enumerate(mrnas):
                        if exclude[i]:
                                continue

                        sorti = set([x[-1] for x in mrna])
                        indexi = set([x[0] for x in mrna])
                        for j,mrna2 in enumerate(mrnas[(i+1):]):
                                indexj = set([x[0] for x in mrna2])
                                sortj = set([x[-1] for x in mrna2])

                                if sortj.issubset(sorti) and max(indexj) >= max(indexi):
                                        exclude[j+i+1] = 1

                mrnas = [x for i,x in enumerate(mrnas) if exclude[i] == 0]




        mrnas = sorted(mrnas, key = lambda x:len(x), reverse = 1)


        return mrnas






def determinegene(aligns_bycontig, exonstotran, exons_infor,contig, transtogene, transtype, allowtruncate):

        size = aligns_bycontig[0][1]

        aligns_exons = cl.defaultdict(list)

        for i,row in enumerate(aligns_bycontig):

                exonid, size, qstart0, qend0, strd0, rstart0, rend0, match0 = row

                transids = exonstotran[exonid]

                for (transid, exonnum) in transids:

                        aligns_exons[transid].append((exonnum, qstart0, qend0,strd0, match0, i))

        allmrnas = []
        for transid, qranges in aligns_exons.items():


                mrnas = getmrna(qranges)


                for mrna in mrnas:

                        allexons = [x[-1] for x in  mrna]
                        allexons = [aligns_bycontig[index] for index in allexons]

                        exoninfor = exons_infor[transid]


                        matchsum = sum([( x[1] -  4*(x[1] - x[-1]) )  for x in allexons])

                        foundexon = set([x[0] for x in mrna])

                        #penalty = 100 * ( len(allowtruncate.get(transid, exons_infor[transid])) - len(foundexon) )

                        scores = [100* ( x[1] -  4*(x[1] - x[-1]) ) /x[1] for x in allexons]
                        score = sum(scores) 

                        factor = 1
                        if transtype[transid] == "protein_coding":
                                score *= 10
                                factor = 10
                        if  len(allowtruncate.get(transid, exons_infor[transid])) <= len(foundexon):
                                score *= 2


                        allexonnum = set([x[0] for x in mrna])

                        mrna = [x[:5]+[s]+x[5:] for x,s in zip(mrna,scores)]

                        allmrnas.append([transid,score,matchsum,mrna])

        allmrnas_filtered = resolveoverlap2(allmrnas,transtype, exons_infor, allowtruncate)


        allmrnas_filtered = [x for x in allmrnas_filtered if len(x[-1])]

        allmrna_new = []
        for row in allmrnas_filtered:

                transid, score, matchsum, mrna = row

                allexons =  [[x[0]]+aligns_bycontig[x[-1]] for x in mrna]

                #matchsum = sum([x[-1] for x in allexons])
                #score = determinescore(row[3], transtype[row[0]] ,len(exons_infor[row[0]]), int(row[0] in allowtruncate))

                allmrna_new.append([transid,score,matchsum ] +allexons)


        allmrna_new = sorted(allmrna_new, key = lambda x: (x[1]*x[2],x[0]), reverse = 1)

        contig_name = "_".join(contig.split("_")[:-2])
        contig_start = int(contig.split("_")[-2])
        contig_end = int(contig.split("_")[-1])

        rows = []
        for mrna in allmrna_new:

                transid = mrna[0]
                score = mrna[1]
                match = mrna[2]

                exoninfor = allowtruncate.get(transid,exons_infor[transid])

                exonnum = len(exoninfor)
                foundexonnum =  len(mrna[3:])

                allindex = ",".join([str(x[0]) for x in mrna[3:]])
                allstarts = ",".join([str(contig_start+x[3]) for x in mrna[3:]])
                allends = ",".join([str(contig_start+x[4]) for x in mrna[3:]])
                allexons = ",".join([str(x[1]) for x in mrna[3:]])
                allstrds = ",".join([str(x[5]) for x in mrna[3:]])

                simi = sum([100*x[-1]/x[-2] for x in mrna[3:]])/exonnum

                match0 = sum([( x[-1] )  for x in mrna[3:]])
                simi0 = sum([100* ( x[-1] ) /x[-2] for x in mrna[3:]])/exonnum

                start = min([contig_start+x[3] for x in mrna[3:]])
                end = max([contig_start+x[4] for x in mrna[3:]])
                totalsize = sum([abs(x[3]-x[4]) for x in mrna[3:]])


                tag = ""
                if exonnum != len(exons_infor[transid]):
                        trunexons = allowtruncate[transid]
                        tag = "(:{}-{})".format(min(trunexons)[0],max(trunexons)[0])


                row = [transid+tag, transtogene[transid][0], transtogene[transid][1], transtype[transid] , match,  simi, exonnum, exonnum - foundexonnum,allindex, contig_name, start ,end , allstarts, allends, allstrds,allexons]
                rows.append(row)

        return rows


def main(args):

        gff3file = args.gff
        inputfile = args.input
        region = args.region.replace('/rc','')
        truncateregion = args.truncate

        truncatecontig = ""
        truncatestart, truncateend = 0,0
        truncatesize = 0

        if truncateregion:
                truncatecontig,truncatelocus = truncateregion.split(":")
                truncatestart, truncateend = truncatelocus.split('-')
                truncatestart, truncateend = int(truncatestart), int(truncateend)
                truncatesize = truncateend - truncatestart

        outpath = args.output

        contig = region.split(':')[0]
        contigstart = int(region.split(':')[1].split('-')[0])
        contigend = int(region.split(':')[1].split('-')[1])
        contigsize = abs(contigend - contigstart)

        gff3 = pd.read_csv(gff3file, header = None, sep = '\t')
        gff3_exons  = gff3 [(gff3[2] == 'exon' )]
        #gff3_trans = gff3 [(gff3[2] == 'transcript' )]


        allowtruncate = dict()


        exons_texts = gff3_exons.values.tolist()
        exons_infor = cl.defaultdict(list)
        exonstotran = cl.defaultdict(list)
        transtype = cl.defaultdict(list)
        transtogene = dict()
        for row in exons_texts:
                text = row[-1]
                text = text.split(';')

                transcript_type  = [x for x in text if x.startswith("transcript_type=")][0][16:]
                exon_id= [x for x in text if x.startswith("exon_id=")][0][8:]
                if "." in exon_id:
                        exon_id = ".".join(exon_id.split(".")[:-1])
                exon_num= int([x for x in text if x.startswith("exon_number=")][0][12:])
                tran_id= [x for x in text if x.startswith("transcript_id=")][0][14:]
                if "." in tran_id:
                        tran_id = ".".join(tran_id.split(".")[:-1])
                gene_id= [x for x in text if x.startswith("gene_id=")][0][8:]

                if "." in gene_id:
                        gene_id = ".".join(gene_id.split(".")[:-1])

                gene_name = [x for x in text if x.startswith("gene_name=")][0][10:]


                transtogene[tran_id] = (gene_id, gene_name)

                exonstotran[exon_id].append((tran_id, len(exons_infor.get(tran_id, [])) + 1))

                exons_infor[tran_id].append([exon_id, abs(row[4] - row[3])])

                transtype[exon_id] = transcript_type
                transtype[tran_id] = transcript_type

                if len(truncateregion) and row[0] == truncatecontig and row[4] - row[3]  + truncatesize - max(row[4], truncateend ) + min(row[3], truncatestart) > 0:

                        if tran_id not in allowtruncate:
                                allowtruncate[tran_id] = []

                        allowtruncate[tran_id].append([exon_num , exon_id])



        table = pd.read_csv(inputfile, header = None, sep = '\t')
        table = table[(table[0].str.contains(contig)) & (table[5].str.contains('ENSE') )]
        table = table.sort_values(by=[9], ascending = 0)

        aligns = table[[0, 5, 6,2,3,4, 7,8,9,12]].values.tolist()

        aligns_bycontigs = cl.defaultdict(list)


        for row in aligns:

                querystart = min(int(row[0].split('_')[-2]), int(row[0].split('_')[-1]))

                contigstart_onquery = contigstart - querystart

                contigend_onquery = contigend - querystart

                if row[4] - row[3]  + contigsize - max(row[4],contigend_onquery ) + min(row[3], contigstart_onquery) <= 0:
                        continue
                if row[-2] < 0.9 * row[2]:
                        continue


                if '.' in row [1]:
                        row[1] = ".".join(row[1].split(".")[:-1])

                aligns_bycontigs[(row[0], row[5])].append(row[1:-1])


        outtable = []
        for (contig,strd), aligns_bycontig in aligns_bycontigs.items():

                aligns_bycontig = resolveoverlap(aligns_bycontig, exonstotran)
                aligns_bycontig = sorted(aligns_bycontig, key = lambda x: x[2])

                if len(aligns_bycontig) ==0:
                        continue

                rows = determinegene(aligns_bycontig, exonstotran,exons_infor,contig,transtogene, transtype, allowtruncate)

                #rows_inregion = [x for x in rows if x[11] - x[10] + contigsize - max(x[11],contigend ) + min(x[10], contigstart) > 0]

                outtable.extend(rows)


        outtable = sorted(outtable, key = lambda x: (x[5]*x[6]*x[4]), reverse = 1)

        pd.DataFrame.from_records(outtable).to_csv(outpath, index = False, header = None, sep = '\t', mode = 'w')

        return 

def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine psuedogene")
        parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
        parser.add_argument("-g", "--gff", help="path to input data file",dest="gff", type=str, required=True)
        parser.add_argument("-r", "--region", help="path to output file", dest="region",type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
        parser.add_argument("-t", "--truncate", help="truncate genes", dest="truncate",type=str, default = "")  

        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
