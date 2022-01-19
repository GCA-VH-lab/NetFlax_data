__author__		= "Chayan Kumar Saha, Gemma C. Atkinson"
__copyright__	= "MIT License: Copyright (c) 2020 Chayan Kumar Saha"
__email__		= "chayan.sust7@gmail.com"

import argparse
import os, sys, os.path, gzip
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC

usage= '''  Description: Filter prdicted toxins based on Abundance '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--input", required=True, help="*_TreeOrderOperonTA.txt")
parser.add_argument("-d", "--desc", required=True, help="*_outdesc.txt")
parser.add_argument("-l", "--list", required=True, help="List of input that was used for TAGs.py")
parser.add_argument("-t", "--threshold", help="threshold here is a number of different species that have same TA cluster, default=2")
parser.add_argument("-v","--version", action="version", version='%(prog)s 1.2.0.21')
parser.add_argument("-o", "--out_prefix", required= True, help=" Any Keyword to define your output eg. MyQuery ")
parser.add_argument("-lp", "--localProteomeDirectory", help=" Path for Local Proteome.faa Files, Default directory is './' which is the same directory where the script is located or running from. ")

args = parser.parse_args()
parser.parse_args()


if args.localProteomeDirectory:
	if os.path.isdir(args.localProteomeDirectory):
		if args.localProteomeDirectory[-1]=='/':
			localP=args.localProteomeDirectory
			#print('Local Proteome.faa Data path : ', localP, '\n')
		else:
			localP=args.localProteomeDirectory+'/'
			#print('Local Proteome.faa Data path : ', localP, '\n')
	else:
		print('No directory Found as : '+ args.localProteomeDirectory)
		sys.exit()
else:
	localP='./'

if args.threshold:
	thresh=args.threshold
else:
	thresh="2"

'''
#Query_Species  QueryAccession  TA_System	   AccessopnInPredictedTA-LikeRegion	   IntergenicSpace ClusterNumber_TA-likeSystem	 Number_Occurred Conserve(%)
WP_031529962.1#244_Dyadobacter_crusticola	   WP_031529962.1  2	   WP_031529962.1|WP_031529961.1   -4	  696|12  31/1062 2.92
WP_047819717.1#290_Croceicoccus_naphthovorans   WP_047819717.1  2	   WP_082134752.1|WP_047819717.1   53	  197|696 2/1062  0.19
WP_027360417.1#478_Desulforegula_conservatrix   WP_027360417.1  2	   WP_027360416.1|WP_027360417.1   45	  3|696   86/1062 8.1
WP_083580626.1#332_Bacteroides_ilei	 WP_083580626.1  2	   WP_083580626.1|WP_083580628.1   90	  696|12  31/1062 2.92
'''
'''
GCF_001639265.1 WP_066970962.1  Methanobrevibacter_filiformis_DSM_11501
'''

def remBadChar(item): #removing characters from species name
	import re
	return re.sub("[^a-zA-Z0-9]"," ",item).replace(" ","_")

taxaGCF={}
i=0
with open(args.list, 'r') as lIn:
	for line in lIn:
		Line=line.rstrip().split('\t')
		i+=1
		taxa=Line[1]+'#'+str(i)+'|'+remBadChar(Line[2])
		taxaGCF[taxa]=Line[0]

#print(taxaGCF)#'WP_004078543.1#908|Methanoplanus_limicola_DSM_2279': 'GCF_000243255.1'
#print(len(taxaGCF))#2281


def cognateFind(item1,item2):
	item1List=item1.split('|')
	#if item2 in item1List:
	item1List.remove(item2)
	return ''.join(map(str,item1List))
	#else:
	#return item1List, item1List.index(item2)

def speciesCheck(item1, item2):
	query=item1.split('#')[0]
	species=item1.split('#')[1][item1.split('#')[1].index('|')+1:]
	TaPairList=item2.split('|')
	IndexQuery=TaPairList.index(query)
	clustIndex=''
	if IndexQuery==0:
		clustIndex=1
	if IndexQuery==1:
		clustIndex=0
	ClustAccession=TaPairList[int(clustIndex)]
	return ClustAccession+'\t'+species

speciesDict={}
LineList=[]
cognateSet=set()

#j=0
with open(args.input, 'r') as TAin:
	for line in TAin:
		if line[0]!='#':
			Line=line.rstrip().split('\t')
			LineList.append(Line)
			#print(Line, cognateFind(Line[3],Line[1]))
			cognateSet.add(cognateFind(Line[3],Line[1])+'\t'+Line[0].split('#')[1].split('|')[0])
			#j+=1
			#print(j, Line,speciesCheck(Line[0],Line[3]))
			#speciesDict[taxaGCF[Line[0]]]=speciesCheck(Line[0],Line[3]).split('\t')[1]
			speciesDict[speciesCheck(Line[0],Line[3]).split('\t')[0]+'\t'+taxaGCF[Line[0]]]=taxaGCF[Line[0]]
			#lenCheckDict[cognateFind(Line[3],Line[1])]=taxaGCF[Line[0]]
#print(LineList) #['WP_052351501.1#784|Parageobacillus_genomosp._1', 'WP_052351501.1', '2', 'WP_052351501.1|WP_052351500.1', '-8', '696|2', '84/1062', '7.91'], ['WP_084455213.1#306_Novosphingobium_rosa', 'WP_084455213.1', '2', 'WP_084455213.1|WP_068083796.1', '61', '696|2', '84/1062', '7.91']
#print(speciesDict)#'WP_014547293.1\tGCF_000024665.1': 'GCF_000024665.1'
#print(len(cognateSet))#WP_014623656.1\t887'}
#print(lenCheckDict)
#print(len(speciesDict)) #1464 #1210
#print(len(cognateSet))#1496
#sys.exit()


def getTA(item):
	query=item[1]
	queryOrder=item[3].split('|').index(query)
	TaOrder=''
	if queryOrder==0:
		TaOrder=1
	elif queryOrder==1:
		TaOrder=0
	TAfam=item[5].split('|')[int(TaOrder)]
	occuranceNum=(item[6].split('/')[0])
	return str(TAfam)+'\t'+occuranceNum

occuranceSet=set()
TADict={}#TAGroup(key)=Clus\tOccur
for items in LineList:
	TADict[items[5]]=getTA(items)
	occuranceSet.add(int(getTA(items).split('\t')[1]))

#print(occuranceSet)
#print(TADict)#'696|13': '13\t10'
#print(len(TADict))
#print(occuranceSet)
##Change here -----
#sys.exit()

TopList=sorted(occuranceSet,reverse=True)

#print(TADict)#'696|13': '13\t10'
#if args.input.split('/')[-1]=='1_DUF4065_netOut_TreeOrderOperonTA.txt':
#	TopList=sorted(occuranceSet,reverse=True)#[:int(topTA)]
#else:
	#TopList=sorted(occuranceSet,reverse=True)[:int(topTA)]
#print(TopList)#[86, 84, 58, 54, 31]

clusterName=set()
for items in TopList:
	for ids in TADict:
		if int(TADict[ids].split('\t')[1])==items:
			#print(items, TADict[ids].split('\t')[1] )
			clusterName.add(int(TADict[ids].split('\t')[0]))

TopCluster=sorted(clusterName)
#print(TopCluster)#[2, 3, 7, 8, 9]
#print(clusterName, 'top')
clusList=[]

with open(args.desc, 'r') as desIn:
	for line in desIn:
		if '\t' in line:
			Line=line.rstrip().split('\t')
			#print(Line[0])
			for cluster in TopCluster:
				#print(cluster)
				if Line[0].split('(')[0]+'('==str(cluster)+'(':
					clusterACC=str(cluster)+'#'+Line[1]
					clusList.append(clusterACC)

def length_from_gwp(gcf_nr,accession_nr):
	faaFile=gcf_nr+'.faa.gz'
	LocalP=localP
	fastaSeq = gzip.open(LocalP+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==accession_nr:
			record.id = gcf_nr+'\t'+accession_nr
			return len(record.seq)



clustDict={}
for item in TopCluster:
	Clus_accList=[]
	for element in clusList:
		if element.split('#')[0]==str(item):
			Clus_accList.append(element.split('#')[1])
	clustDict[str(item)]=Clus_accList

#print(clustDict) #length of the accession wasn't checked
#print(clustDict) #'2': ['WP_054673950.1', 'WP_090147093.1', 'WP_113651984.1'


speciesClusDict={}
clusLenDict={} #similar to clustDict but sorted if protein <=400 aa
for ids in clustDict:
	speciesClustSet=set()
	clusLenList=[]
	for items in clustDict[ids]:
		for keys in speciesDict:
			if keys.split('\t')[0]==items:
				#print(items, speciesDict[keys], length_from_gwp(speciesDict[keys],items))
				if length_from_gwp(speciesDict[keys],items)<=450: #LengthLimit450AA
					clusLenList.append(items)
					speciesClustSet.add(speciesDict[keys])
	clusLenDict[ids]=clusLenList
	speciesClusDict[ids]=len(speciesClustSet)


cwd=os.getcwd()
RoundName=str(cwd).split('/')[-1].split('_')[0]
#print(cwd)

for ids in speciesClusDict:
	if int(speciesClusDict[ids])>=int(thresh):
		with open ('G'+ids+'#'+args.out_prefix+'.clus', 'w') as topOut, open('G'+ids+'#'+args.out_prefix+'.clusNum', 'w') as topNout:
			for items in clusLenDict[ids]:
				#print(items, clustDict[ids])
				for cogitems in cognateSet:
					if items==cogitems.split('\t')[0]:
						#print(items, cogitems)
						print(items, file=topOut)
						print(cogitems, file=topNout)
	if int(speciesClusDict[ids])<int(thresh) and int(speciesClusDict[ids])>=2:
		with open ('G'+ids+'#'+args.out_prefix+'.cutThresh', 'w') as cutThreshOut:
			for items in clusLenDict[ids]:
				#print(items, clustDict[ids])
				for cogitems in cognateSet:
					if items==cogitems.split('\t')[0]:
						print(items, file=cutThreshOut)

'''
for ids in speciesClusDict:
	if int(speciesClusDict[ids])>=int(thresh):
		with open ('G'+ids+'#'+args.out_prefix+'.clus', 'w') as topOut, open('G'+ids+'#'+args.out_prefix+'.clusNum', 'w') as topNout:
			for items in clustDict[ids]:
				#print(items, clustDict[ids])
				for cogitems in cognateSet:
					if items==cogitems.split('\t')[0]:
						#print(items, cogitems)
						print(items, file=topOut)
						print(cogitems, file=topNout)
	if int(speciesClusDict[ids])<int(thresh) and int(speciesClusDict[ids])>=2:
		with open ('G'+ids+'#'+args.out_prefix+'.cutThresh', 'w') as cutThreshOut:
			for items in clustDict[ids]:
				#print(items, clustDict[ids])
				for cogitems in cognateSet:
					if items==cogitems.split('\t')[0]:
						print(items, file=cutThreshOut)

'''














'''

python3 getTopTAFam.py -i big2281SpPanA150_TreeOrderOperonTA.txt -d big2281SpPanA150_flankgene.fasta_cluster_out_3_1e-10_outdesc.txt -t 2 -o big2281SpPanA150 > discradedList_PanT.txt
./discradedList_PanT.txt
'''
#	print(ids, len(clustDict[ids]))
