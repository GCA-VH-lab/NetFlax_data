#Author: Chayan Kumar Saha
import argparse
import gzip, time
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC
import os, sys, os.path

usage= '''  Description: Generate Filtered Output from hmmParsed File and extract sequences '''


parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--input", required=True, help="hmmParsed.txt")
parser.add_argument("-e", "--ethreshold", help=" E value threshold. Default = Top 5000 hits ")
parser.add_argument("-se", "--skipethreshold", action="store_true", help=" Skip E value threshold")
parser.add_argument("-l", "--limit", help=" rch filtered Limit. Default = 1000 ")
parser.add_argument("-lp", "--localProteomeDirectory", help=" Path for Local Proteome.faa Files, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-k", "--keep", action="store_true", help=" If you want to keep the intermediate files use [-k]. By default it will remove. ")
parser.add_argument("-v","--version", action="version", version='%(prog)s 1.5.0.21')
args = parser.parse_args()
parser.parse_args()

'''
hmmParsed.txt
GCF_902386925.1	WP_080800263.1	G6#DUF4065	1e-24	1.1e-24	84.4	7	123	1	123	toxin [Arabia massiliensis]
GCF_900130025.1	WP_200796512.1	G6#DUF4065	5.5e-15	6e-15	52.9	17	95	3	109	hypothetical protein
'''

if args.localProteomeDirectory:
	if os.path.isdir(args.localProteomeDirectory):
		if args.localProteomeDirectory[-1]=='/':
			localP=args.localProteomeDirectory
			print('Local Proteome.faa Data path : ', localP, '\n')
		else:
			localP=args.localProteomeDirectory+'/'
			print('Local Proteome.faa Data path : ', localP, '\n')
	else:
		print('No directory Found as : '+ args.localProteomeDirectory)
		sys.exit()
else:
	localP=os.getcwd()+'/'


numLimit=5000

if args.limit:
	ulimit=args.limit
else:
	ulimit="1000"


def isfloat(item):
	try:
		inS=float(item)
		return 'Yes'
	except ValueError:
		return 'No'

def length_from_gwp(gcf_nr,accession_nr):
	faaFile=gcf_nr+'.faa.gz'
	LocalP=localP
	fastaSeq = gzip.open(LocalP+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==accession_nr:
			record.id = gcf_nr+'\t'+accession_nr
			return len(record.seq)

currnetDir=os.getcwd()
sqlDB=currnetDir+'/'+args.input.split('.')[0]+'.sqlite'

import sqlite3 # Import the module
conn = sqlite3.connect(sqlDB) # Make a connection object
cursorObject = conn.cursor() # Create a cursor object
# The trailing comma is used to tell python that this is a tuple:
dropTableStatement = "DROP TABLE IF EXISTS hmmsearchP"
cursorObject.execute(dropTableStatement)
#GCF_002216145.1	WP_089035163.1	G19#DUF4065	2.3e-39	4.1e-39	131.6	4	191	1	197	GTP pyrophosphokinase [Neisseria chenwenguii]
#GCF_900080235.1	XP_020502886.1	eEF1A_linsi	1.2e-268	1.4e-268	892.0	1	458	1	459	elongation factor 1-alpha [Labrus bergylta
createTable = "CREATE TABLE hmmsearchP (NCBI_Accnr text, Target text, Query text, E_Value real, I_E_Value real) ;"
cursorObject.execute(createTable)
#cur.execute("INSERT INTO stocks VALUES ('2006-01-05','BUY','RHAT',100,35.14)")

eLineList=[]
iEvalueList=[]

with open(args.input, 'r') as fileIn:
	for line in fileIn:
		Line=line.rstrip().split('\t')
		#if length_from_gwp(Line[0],Line[1])<=400: #LengthLimit400AA
		addLine='INSERT INTO hmmsearchP VALUES ('+'"'+Line[0]+'"'+','+'"'+Line[1]+'"'+','+'"'+Line[2]+'"'+','+Line[3]+','+Line[4]+');'
		#print(addLine)
		cursorObject.execute(addLine)
		if Line and len(Line)>10:
			if isfloat(Line[4])=='Yes':
				eLineList.append(Line)
				iEvalueList.append(float(Line[4]))
			else:
				print(Line)
		else:
			print(Line)

conn.commit()

conn = sqlite3.connect(sqlDB) # Make a connection object
cursorObject = conn.cursor() # Create a cursor object

#queryList=[]
queryList=set()
# Prepare SQL query to retrieve a record from the database.
uniqQuery = 'SELECT DISTINCT Query FROM hmmsearchP;'
try:
	# Execute the SQL command
	cursorObject.execute(uniqQuery)
	# Fetch all the rows in a list of lists.
	queries = cursorObject.fetchall()
	for row in queries:
		#queryList.append("\t".join(map(str, row)))
		queryList.add("\t".join(map(str, row)))
except:
	print ("Error: unable to fetch data")

#spList=[]
spList=set()
# Prepare SQL query to retrieve a record from the database.
uniqSp = 'SELECT DISTINCT NCBI_Accnr FROM hmmsearchP;'
try:
	# Execute the SQL command
	cursorObject.execute(uniqSp)
	# Fetch all the rows in a list of lists.
	Sps = cursorObject.fetchall()
	for row in Sps:
		#spList.append("\t".join(map(str, row)))
		spList.add("\t".join(map(str, row)))
except:
	print ("Error: unable to fetch data")

#targetList=[]
targetList=set()
uniqTarget = 'SELECT DISTINCT Target FROM hmmsearchP;'
try:
	# Execute the SQL command
	cursorObject.execute(uniqTarget)
	# Fetch all the rows in a list of lists.
	targets = cursorObject.fetchall()
	for row in targets:
		#targetList.append("\t".join(map(str, row)))
		targetList.add("\t".join(map(str, row)))
except:
	print ("Error: unable to fetch data")

print('DatabaseMade', time.ctime())

percentileJack=0
FilteredHitList=[]
for item in sorted(spList):
	percentileJack+=1
	if percentileJack % 1000 == 0:
		print('\t'+'>>> '+str(round(int(percentileJack)*100/len(spList)))+'%'+' Completed...'+'('+str(percentileJack)+'/'+str(len(spList))+')')
	if percentileJack % len(spList) == 0:
		print('\t'+'>>> Completed ' +'\n')
	for element in queryList:
		bestHit = 'SELECT * FROM hmmsearchP WHERE NCBI_Accnr LIKE "%' + item + '%" AND Query LIKE "%' + element + '%"'
		cursorObject.execute(bestHit)
		results = cursorObject.fetchall()
		for ids in targetList:
			hitList=[]
			i_evalueList=[]
			for row in results:
				if ids == row[1]:
					hitList.append("\t".join(map(str, row)))
					i_evalueList.append(float(row[4]))
			if len(hitList)>0 and len(i_evalueList)>0:
				#print(hitList[i_evalueList.index(min(i_evalueList))])
				FilteredHitList.append(hitList[i_evalueList.index(min(i_evalueList))].split('\t'))
print('HitListMade', time.ctime())
cursorObject.close()

sp_set=set()
mysql_Dict={} #spIds+'\t'+proteinAccnr+'|'+domainHmm as keys
for hits in FilteredHitList:
	sp_set.add(hits[0])
	mysql_Dict[hits[0]+'\t'+hits[1]+'|'+hits[2]]= ("\t".join(map(str, hits)))

if args.ethreshold:
	evthresh=args.ethreshold
else:
	if len(eLineList)>numLimit:
		evthresh=str(sorted(iEvalueList)[numLimit])
	else:
		evthresh="1e-10"


netFlaxInputList=[] #['GCF_000008185.1', 'NP_972407.1']
if not args.skipethreshold:
	with open (args.input.split('.')[0]+'_'+str(evthresh)+'_filtered.txt', 'w') as out:
		for keys in mysql_Dict:
			#if keys[keys.index('|')+1:].upper()=='SYNTH':
			if float(mysql_Dict[keys].split('\t')[4])<float(evthresh) or float(mysql_Dict[keys].split('\t')[4])==float(evthresh):
				print(mysql_Dict[keys], file=out)
				netFlaxInputList.append(mysql_Dict[keys].split('\t')[:2])
else:
	with open (args.input.split('.')[0]+'_skipEval_filtered.txt', 'w') as sout:
		for keys in mysql_Dict:
			#if keys[keys.index('|')+1:].upper()=='SYNTH':
			#if float(mysql_Dict[keys].split('\t')[4])<float(evthresh) or float(mysql_Dict[keys].split('\t')[4])==float(evthresh):
			print(mysql_Dict[keys], file=sout)
			netFlaxInputList.append(mysql_Dict[keys].split('\t')[:2])
#print(netFlaxInputList[0])

def seq_from_wp_G(gcf_nr,accession_nr):
	faaFile=gcf_nr+'.faa.gz'
	#LocalP='/data/users/chayan/chayanBioPhD/netFlaX_AutoMation/LocalProteomes/'
	fastaSeq = gzip.open(localP+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==accession_nr:
			record.id = gcf_nr+'\t'+accession_nr
			record.description = str('')
			#handle.close()
			return record.format("fasta")

if not args.skipethreshold:
	if len(netFlaxInputList)>int(ulimit):
		netInFileName=args.input.split('.')[0]+'.netIn'
		fastaInFileName=args.input.split('.')[0]+'_'+str(evthresh)+'_filtered.fasta'
		with open (fastaInFileName, 'w') as fastaout:
			for item in netFlaxInputList:
				if length_from_gwp(item[0], item[1])<=450: #LengthLimit450AA
					print(seq_from_wp_G(item[0], item[1]), file=fastaout)
		outFileList=[]
		for i in range(1,10):
			identity=str('0.'+str(i))
			fastaOutFileName=args.input.split('.')[0]+'_'+str(evthresh)+'_filtered'+'_'+identity+'.fasta'
			outFileList.append(fastaOutFileName)
			command="usearch11.0.667_i86linux32 -cluster_fast %s -id %s -centroids %s" %(fastaInFileName, identity, fastaOutFileName)
			#command="usearch11.0.667_i86osx32 -cluster_fast %s -id %s -centroids %s" %(fastaInFileName, identity, fastaOutFileName)
			print(command)
			os.system(command)
		if not args.keep:
			if os.path.isfile(fastaInFileName):
				os.remove(fastaInFileName)
		fasDict={}
		for fastaFiles in outFileList:
			fileList=[]
			with open(fastaFiles, 'r') as fastaFileIn:
				for line in fastaFileIn:
					if line[0]=='>':
						fileList.append(line[1:].rstrip())
			fasDict[fastaFiles]=fileList

		maxlenList=[]
		minlenList=[]
		for item in fasDict:
			print(item, len(fasDict[item]))
			if len(fasDict[item])<int(ulimit)+1:
				maxlenList.append(len(fasDict[item]))
		if len(maxlenList)>0:
			for item in fasDict:
				if len(fasDict[item])==max(maxlenList):
					with open (netInFileName, 'w') as netout:
						print(("\n".join(map(str, fasDict[item]))), file=netout)

		else:
			for item in fasDict:
				minlenList.append(len(fasDict[item]))
			for item in fasDict:
				if len(fasDict[item])==min(minlenList):
					with open (netInFileName, 'w') as netout:
						print(("\n".join(map(str, fasDict[item]))), file=netout)

		if not args.keep:
			for item in fasDict:
				if os.path.isfile(item):
					os.remove(item)
	else:
		netInFileName=args.input.split('.')[0]+'.netIn'
		with open (netInFileName, 'w') as netout:
			for item in netFlaxInputList:
				if length_from_gwp(item[0], item[1])<=450: #LengthLimit450AA
					print(item[0], item[1], sep='\t', file=netout)
else:
	if len(netFlaxInputList)>int(ulimit):
		netInFileName=args.input.split('.')[0]+'.netIn'
		fastaInFileName=args.input.split('.')[0]+'_skipEval_filtered.fasta'
		with open (fastaInFileName, 'w') as fastaout:
			for item in netFlaxInputList:
				if length_from_gwp(item[0], item[1])<=450: #LengthLimit450AA
					print(seq_from_wp_G(item[0], item[1]), file=fastaout)
		outFileList=[]
		for i in range(1,10):
			identity=str('0.'+str(i))
			fastaOutFileName=args.input.split('.')[0]+'_skipEval_filtered'+'_'+identity+'.fasta'
			outFileList.append(fastaOutFileName)
			command="usearch11.0.667_i86linux32 -cluster_fast %s -id %s -centroids %s" %(fastaInFileName, identity, fastaOutFileName)
			#command="usearch11.0.667_i86osx32 -cluster_fast %s -id %s -centroids %s" %(fastaInFileName, identity, fastaOutFileName)
			print(command)
			os.system(command)
		if not args.keep:
			if os.path.isfile(fastaInFileName):
				os.remove(fastaInFileName)
		fasDict={}
		for fastaFiles in outFileList:
			fileList=[]
			with open(fastaFiles, 'r') as fastaFileIn:
				for line in fastaFileIn:
					if line[0]=='>':
						fileList.append(line[1:].rstrip())
			fasDict[fastaFiles]=fileList

		maxlenList=[]
		minlenList=[]
		for item in fasDict:
			print(item, len(fasDict[item]))
			if len(fasDict[item])<int(ulimit)+1:
				maxlenList.append(len(fasDict[item]))
		if len(maxlenList)>0:
			for item in fasDict:
				if len(fasDict[item])==max(maxlenList):
					with open (netInFileName, 'w') as netout:
						print(("\n".join(map(str, fasDict[item]))), file=netout)

		else:
			for item in fasDict:
				minlenList.append(len(fasDict[item]))
			for item in fasDict:
				if len(fasDict[item])==min(minlenList):
					with open (netInFileName, 'w') as netout:
						print(("\n".join(map(str, fasDict[item]))), file=netout)

		if not args.keep:
			for item in fasDict:
				if os.path.isfile(item):
					os.remove(item)
	else:
		netInFileName=args.input.split('.')[0]+'.netIn'
		with open (netInFileName, 'w') as netout:
			for item in netFlaxInputList:
				if length_from_gwp(item[0], item[1])<=450: #LengthLimit450AA
					print(item[0], item[1], sep='\t', file=netout)
