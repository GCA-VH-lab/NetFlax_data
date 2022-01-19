#Author:Chayan Kumar Saha

import time, math
import re
import urllib.request
import argparse
import glob
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC
import os, sys, os.path, subprocess
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

usage= '''  Description: NetFlaX Automation Script '''
#python3 NetFlaX_automationPipeline.py -lp proteomes/ -lh HMMs/ -ld LocalProteomes/ -c 40 -t 10  -n 5

parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-ld", "--localGFFDirectory", help=" Path for Local Gff3.gz and Proteome.gz Files, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-lp", "--localProteomeDirectory", help=" Path for Local Proteome.faa Files, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-lh", "--localHmmDirectory", help=" Path for Local Hmm Files, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-ls", "--localScriptDirectory", help=" Path for Local Scripts, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-c", "--cpu", help="Maximum number of parallel CPU workers to use for multithreads. ")
parser.add_argument("-e", "--ethreshold", help=" E value threshold. Default = 1e-10 ")
parser.add_argument("-t", "--threshold", help=" TA species threshold. Default = 10 ")
parser.add_argument("-nt", "--nthreshold", help=" NetFlaX query threshold. Default = 1000 ")
parser.add_argument("-n", "--number", help=" Minimum cut-off of accepted TA-like cluster per round for next round. Default = 5 ")
parser.add_argument("-v","--version", action="version", version='%(prog)s MT1.0.4')
args = parser.parse_args()
parser.parse_args()



if args.ethreshold:
	evthreshTADB=args.ethreshold
else:
	evthreshTADB="1e-10"

if args.nthreshold:
	nthresh=args.nthreshold
else:
	nthresh="1000"

if args.threshold:
	spthresh=args.threshold
else:
	spthresh="10"

if args.number:
	number=args.number
else:
	number="5"

if args.cpu:
	if int(args.cpu)>0:
		core=int(args.cpu)
	else:
		print('Please use number eg, 1,2...')
		sys.exit()

if args.localGFFDirectory:
	if os.path.isdir(args.localGFFDirectory):
		if args.localGFFDirectory[-1]=='/':
			localPG=args.localGFFDirectory
			print('Local GFF and Proteome Data path : ', localPG, '\n')
		else:
			localPG=args.localGFFDirectory+'/'
			print('Local GFF and Proteome Data path : ', localPG, '\n')
	else:
		print('No directory Found as : '+ args.localGFFDirectory)
		sys.exit()
else:
	localPG='./'

#beginDir=os.getcwd()
homeDir=os.getcwd()+'/'

def XConFind(listName):
	discardHmm=[]
	for item in listName: #G1#DUF4065.hmm G2#DUF4065
		if '#'.join(map(str,item[1].split('.')[0].split('#')[1:])) != '#'.join(map(str,item[2].split('#')[1:])):
			discardHmm.append(item[1])
	if discardHmm:
		if len(discardHmm)==0:
			return 0
		elif len(discardHmm)==1:
			return 1 #crossfound
		else:
			return len(discardHmm) #Multiple cross with prev hmm
	else:
		return 0 #no cross verified


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
	localP='./'

if args.localHmmDirectory:
	if os.path.isdir(args.localHmmDirectory):
		if args.localHmmDirectory[-1]=='/':
			localH=args.localHmmDirectory
			print('Local HMMs Data path : ', localH, '\n')
		else:
			localH=args.localHmmDirectory+'/'
			print('Local HMMs Data path : ', localH, '\n')
	else:
		print('No directory Found as : '+ args.localHmmDirectory)
		sys.exit()
else:
	localH='./'

if args.localScriptDirectory:
	if os.path.isdir(args.localScriptDirectory):
		if args.localScriptDirectory[-1]=='/':
			localS=args.localScriptDirectory
			print('Local Script path : ', localS, '\n')
		else:
			localS=args.localScriptDirectory+'/'
			print('Local Script path : ', localS, '\n')
	else:
		print('No directory Found as : '+ args.localScriptDirectory)
		sys.exit()
else:
	localS='./'

def checkHmmDirName(item):
	if '/' in item:
		return item.split('/')[1][:-4]
	else:
		return item[:-4]

def renameDir(item):
	if '#' in item:
		return item.split('#')[0]
	else:
		return item

def getPrevFromHMM(item): #G10R4#G2R3#G1R2#G16R1#DUF4065.hmm
	if '#' in item:
		splitItem= item.split('.')[0].split('#')
		return '#'.join(map(str,splitItem[1:]))



proteomeList=[]
for items in (glob.glob(localP+"*.faa")):
	proteomeList.append(items)

import queue as Queue

from queue import Queue,Empty
import threading, os

def worker_func():
	while not stopped.is_set():
		try:
			# use the get_nowait() method for retrieving a queued item to
			# prevent the thread from blocking when the queue is empty
			com = q.get_nowait()
		except Empty:
			continue
		try:
			os.system(com)
		except Exception as e:
			print('Error running command', str(e))
		finally:
			q.task_done()


for i in range(500):
	if (len(glob.glob(localH+"*.hmm")))>0:
		i+=1
		TAnetRes=set()
		print('\n')
		print('###TAsearch Round', i)
		print('\n')
		with open ('Round'+str(i)+'_netfile.txt', 'w') as deadOut:
			hmmList=[]
			hmmOutDirList=[]
			pCom1=[]
			pCom1B=[]
			for hmms in (glob.glob(localH+"*.hmm")):
				hmmList.append(hmms)
				HmmOutDir = str(i)+'_'+checkHmmDirName(hmms)
				HmmUsedDir = 'Round_'+str(i)+'_HMMs'
				hmmOutDirList.append(HmmOutDir)
				if not os.path.exists(HmmOutDir):
					os.makedirs(HmmOutDir)
				if not os.path.exists(HmmUsedDir):
					os.makedirs(HmmUsedDir)
				parsedTextFile=HmmOutDir+'/'+HmmOutDir+'_parsed.txt'
				#print(parsedTextFile)
				#for tabs in (glob.glob(outDir+'/'+"*_tab")):
					#print(tabs)
				scriptParse=localS+'hmmerParse.py'
					#prevVer#parseCommand="python3 ~/chayanBioPhD/netFlaX_AutoMation/hmmerParse.py -i %s >> %s" %(tabs, parsedTextFile)
				shList=[]
				for items in proteomeList:
					#print(items)

					tabFile=HmmOutDir+'/'+items.split('/')[1][:-4]+'_tab'
					outFile=HmmOutDir+'/'+items.split('/')[1][:-4]+'_out'
					if hmms.split('/')[1]=='DUF4065.hmm':
						#print(hmms.split('/')[1])
						command="hmmsearch --cpu 1 --domtblout %s --noali --cut_ga %s %s > %s" %(tabFile, hmms, items, outFile)
						parseCommand="python3 %s -i %s >> %s" %(scriptParse, tabFile, parsedTextFile)
						pCom1.append(command)
						shList.append(parseCommand)
					else:
						command="hmmsearch --cpu 1 --domtblout %s --noali --incE 1e-10 --incdomE 1e-10 -E 1e-10 %s %s > %s" %(tabFile, hmms, items, outFile)
						parseCommand="python3 %s -i %s >> %s" %(scriptParse, tabFile, parsedTextFile)
						pCom1.append(command)
						shList.append(parseCommand)
				with open(HmmOutDir+'/'+HmmOutDir+'_parsed.sh', 'w') as shout:
					for elem in shList:
						print(elem, file=shout)
				shcom="bash %s" %(HmmOutDir+'/'+HmmOutDir+'_parsed.sh')
				pCom1B.append(shcom)
					#with open(HmmOutDir+'/'+items.split('/')[1].replace('.faa','.sh'), 'w') as shout:
						#for elem in shList:
							#print(elem, file=shout)

					#shcom="bash %s" %(HmmOutDir+'/'+items.split('/')[1].replace('.faa','.sh'))
					#pCom1.append(shcom)


			thread_count = core # maximum parallel threads
			stopped = threading.Event()
			q = Queue()
			print('-- Processing : HmmSearch with '+ (str(len(pCom1))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')
			for item in pCom1:
				q.put(item)
			for x in range(thread_count):
				t = threading.Thread(target=worker_func)
				# t.daemon = True #Enable to run threads as daemons
				t.start()
			q.join()       # block until all tasks are done
			stopped.set()
			print('Process HmmSearch Done', '\n')
			#print(pCom1B)
			thread_count = core # maximum parallel threads
			stopped = threading.Event()
			q = Queue()
			print('-- Processing : HmmParse with '+ (str(len(pCom1B))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')
			for item in pCom1B:
				q.put(item)
			for x in range(thread_count):
				t = threading.Thread(target=worker_func)
				# t.daemon = True #Enable to run threads as daemons
				t.start()
			q.join()       # block until all tasks are done
			stopped.set()
			print('Process HmmParsing Done', '\n')


			for hmms in (glob.glob(localH+"*.hmm")):
				usedHMMFile= HmmUsedDir+'/'+checkHmmDirName(hmms)+'.hmm'
				renameHMMs= "mv %s %s" %(hmms, usedHMMFile)
				os.system(renameHMMs)
			rmCom3=[]
			for outDir in hmmOutDirList:
				#print(outDir)
				for GCFs in (glob.glob(outDir+'/'+"GC*_*")):
					removeFile="rm %s" %(GCFs)
					#print(removeFile)
					rmCom3.append(removeFile)
			thread_count = core # maximum parallel threads
			stopped = threading.Event()
			q = Queue()
			print('-- Processing : Removing HMMout with '+ (str(len(rmCom3))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')
			for item in rmCom3:
				q.put(item)
			for x in range(thread_count):
				t = threading.Thread(target=worker_func)
				# t.daemon = True #Enable to run threads as daemons
				t.start()
			q.join()       # block until all tasks are done
			stopped.set()
			print('Process Removing Done','\n')

			pCom3=[]
			#print('\n','Process Filtering starting..')
			outDirCount=set()
			pcount=set()
			for outDir in hmmOutDirList:
				outDirCount.add(outDir)
				parsedTextFile=outDir+'/'+outDir+'_parsed.txt'
				if outDir.split('_')[1]=='DUF4065':
					#FilterCommand="python3 ~/chayanBioPhD/netFlaX_AutoMation/filterHmmParsedFile.py -se -i %s -l %s" %(parsedTextFile, int(nthresh))  #LimitForFlaGsRunto 1000
					scriptFilterHMMParse=localS+'filterHmmParsedFile.py'
					FilterCommand="python3 %s -se -i %s -l %s -lp %s" %(scriptFilterHMMParse, parsedTextFile, int(nthresh), localPG)#LimitForFlaGsRunto 1000
					pCom3.append(FilterCommand)
				else:
					#FilterCommand="python3  ~/chayanBioPhD/netFlaX_AutoMation/filterHmmParsedFile.py -i %s -l %s" %(scriptFilterHMMParse, parsedTextFile, int(nthresh))
					scriptFilterHMMParse=localS+'filterHmmParsedFile.py'
					FilterCommand="python3 %s -i %s -l %s -lp %s" %(scriptFilterHMMParse, parsedTextFile, int(nthresh), localPG)
					pCom3.append(FilterCommand)
			thread_count = core # maximum parallel threads
			stopped = threading.Event()
			q = Queue()
			#print(pCom3)
			print('-- Processing : FilterHMMs with '+ (str(len(pCom3))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')
			for item in pCom3:
				q.put(item)
			for x in range(thread_count):
				t = threading.Thread(target=worker_func)
				# t.daemon = True #Enable to run threads as daemons
				t.start()
			q.join()       # block until all tasks are done
			stopped.set()
			for outDir in hmmOutDirList:
				for ftext in (glob.glob(outDir+'/'+"*_filtered.txt")):
					pcount.add(ftext)
			print('\t'+str(len(pcount))+'/'+str(len(outDirCount))+'...processed')
			print('\n'+'Process Filtering: Done'+'\n'+'\n')
					#pCom3.append(FilterCommand)

			outDirFCount=set()
			pFcount=set()
			if len(outDirCount)==len(pcount):
				def coreF(num):
					import math
					if num==1:
						return core
					if num>1:
						if (core//num)<2:
							return 2
						else:
							return (core//num)
				pCom4=[]
				for outDir in hmmOutDirList:
					outDirFCount.add(outDir)
					netflaxDir=outDir+'/'
					netflaxFile=outDir+'_parsed.netIn'
					if netflaxFile=='1_DUF4065_parsed.netIn':
						panA_filename=homeDir+"panA_TA.txt"
						panA_filenameOut=homeDir+outDir+'/'+outDir+'_parsedLab.net2In'
						labUnet2In=homeDir+outDir+'/'+outDir+'_parsedLabU.net2In'
						panA_NetCommand="cat %s | cut -f 1,2 > %s"%(panA_filename, panA_filenameOut)
						panA_NetCommand2="grep -vf %s %s > %s"%(homeDir+netflaxDir+netflaxFile, panA_filenameOut,labUnet2In)
						panA_NetCommand3="cat %s >> %s"%(labUnet2In, homeDir+netflaxDir+netflaxFile)
						os.system(panA_NetCommand)
						os.system(panA_NetCommand2)
						os.system(panA_NetCommand3)
						#print('PanA_TestQueryAddedInNet')
					netflaxFileOut=outDir+'_netOut'
					os.chdir(homeDir+netflaxDir)
					cwd=os.getcwd()
					#print(cwd)
					netflaxScript=homeDir+localS+'TAGs.py'
					coreTag=coreF(len(pcount))
					netFlaxCommand="python3 %s -l %s -ld %s -i %s -tl %s -n %s -o %s -c %s -u cksaha.mbio@gmail.com" %(netflaxScript, netflaxFile,homeDir+localPG,100, homeDir+'TaxaGCF.txt', 7, netflaxFileOut, coreTag) #Intergenic space 100
					#print(netFlaxCommand) #/data/users/chayan/chayanBioPhD/netFlaX_AutoMation/1_DUF4065
					#netFlaxCommand="python3 ~/chayanBioPhD/netFlaX_AutoMation/netFlaX.py -l %s -ld %s -i %s -t -to -o %s -c %s -vb" %(netflaxFile,'../'+localPG,100,netflaxFileOut,core) #Intergenic space 100
					##os.system(netFlaxCommand)
					commandLineFile=outDir+'_TAGsCommand.sh'
					with open (commandLineFile, 'w') as shOut:
						print('cd '+cwd,file=shOut)
						print('pwd', file=shOut)
						print(netFlaxCommand,file=shOut)
						print('cd ../',file=shOut)
						print('pwd', file=shOut)
					bashNetCom='/bin/bash '+homeDir+netflaxDir+'*_TAGsCommand.sh'
					pCom4.append(bashNetCom)
					cwdback='/'.join(map(str,cwd.split('/')[:-1]))
					os.chdir(cwdback)
				thread_count = core # maximum parallel threads
				stopped = threading.Event()
				q = Queue()
				#print(pCom3)
				print('-- Processing : FlaGs with '+ (str(len(pCom4))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')
				for item in pCom4:
					q.put(item)
				for x in range(thread_count):
					t = threading.Thread(target=worker_func)
					# t.daemon = True #Enable to run threads as daemons
					t.start()
				q.join()       # block until all tasks are done
				stopped.set()
				for outDir in hmmOutDirList:
					for outDesc in (glob.glob(homeDir+outDir+'/'+outDir+'_netOut_flankgene.fasta_cluster_out_*_outdesc.txt')):
						pFcount.add(outDesc)
			print(str(len(pFcount))+'/'+str(len(outDirFCount))+'...processed')
			print('\n'+'Process FlaGs: Done'+'\n'+'\n')
			if len(outDirFCount)==len(pFcount):
				pCom5=[]
				for outDir in hmmOutDirList:
					netflaxDir=outDir+'/'
					os.chdir(homeDir+netflaxDir)
					cwd=os.getcwd()
					clusFileExtensionName=str(cwd).split('/')[-1].split('_')[1]
					netFlax_TA=homeDir+outDir+'/'+outDir+'_netOut_*peronTA.txt'
					netFlax_TAtsv=homeDir+outDir+'/'+outDir+'_netOut_*peronTA.tsv'
					netFlax_taxaList=outDir+'_netOut_taxaList.txt'
					netFlax_outdesc=homeDir+outDir+'/'+outDir+'_netOut_flankgene.fasta_cluster_out_*_outdesc.txt'
					getClusScript=homeDir+localS+'getTAFam.py'
					topTAcommand="python3 %s -i %s -d %s -t %s -l %s -lp %s -o %s " %(getClusScript, netFlax_TA, netFlax_outdesc, int(spthresh), homeDir+outDir+'/'+netFlax_taxaList, homeDir+localPG, clusFileExtensionName) #threshold -t '2' needs to be 10
					tacommandLineFile=outDir+'_TAclus_Command.sh'
					with open (tacommandLineFile, 'w') as tashOut:
						print('cd '+cwd,file=tashOut)
						print('pwd', file=tashOut)
						print(topTAcommand,file=tashOut)
						print('cd ../',file=tashOut)
						print('pwd', file=tashOut)
					tabashNetCom='/bin/bash '+homeDir+netflaxDir+'*_TAclus_Command.sh'
					pCom5.append(tabashNetCom)
					cwdback='/'.join(map(str,cwd.split('/')[:-1]))
					os.chdir(cwdback)
				thread_count = core # maximum parallel threads
				stopped = threading.Event()
				q = Queue()
				#print(pCom3)
				print('-- Processing : TA-Like Cluster finding with '+ (str(len(pCom5))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')
				for item in pCom5:
					q.put(item)
				for x in range(thread_count):
					t = threading.Thread(target=worker_func)
					# t.daemon = True #Enable to run threads as daemons
					t.start()
				q.join()       # block until all tasks are done
				stopped.set()
			print('\t'+str(len(pCom5))+'/'+str(len(outDirFCount))+'...processed')
			print('\n'+'Process TA-Like Cluster finding: Done'+'\n')
			if len(outDirFCount)==len(pFcount): #not necessary but for checking maintenance
				'''
				###fileForCheck
				GCF_904420025.1	WP_191552224.1	DUF4065
				GCF_904423815.1	WP_191556125.1	DUF4065
				GCF_000006865.1	WP_010906205.1	G6#DUF4065
				GCF_000010045.1	WP_011461871.1	G6#DUF4065
				GCF_000011045.1	WP_101495023.1	G6#DUF4065
				'''
				crossCheckSet=set()
				#crosscheckAccGroupDict={}
				#crosscheckNameList=[]
				cOutDir=0
				for outDir in hmmOutDirList:
					print('\n'+'-- Cross Checking Directory: '+outDir+'\n')
					cOutDir+=1
					netflaxDir=outDir+'/'
					os.chdir(homeDir+netflaxDir)
					cwd=os.getcwd()
					netFlax_TAtsv=homeDir+outDir+'/'+outDir+'_netOut_*peronTA.tsv'
					fileForCheck=homeDir+outDir+'/'+outDir+'_parsed_*_filtered.txt'
					#fileForCheck=homeDir+outDir+'/'+outDir+'_parsed_*_filtered.txt' #rethink
					crosscheck_file=homeDir+"crossCheckList.txt"
					if '1_DUF4065'==str(outDir):
						panA_file=homeDir+"panA_TA.txt"
						panA_addCommand="cat %s > %s" %(panA_file, crosscheck_file)
						os.system(panA_addCommand)
					crossCheckFileCommand="cat %s | cut -f 1,2,3 >> %s" %(fileForCheck, crosscheck_file)
					os.system(crossCheckFileCommand)

					currnetDir=os.getcwd()
					sqlDB=currnetDir+'/crossCheckList.sqlite'

					import sqlite3 # Import the module
					conn = sqlite3.connect(sqlDB) # Make a connection object
					cursorObject = conn.cursor() # Create a cursor object
					# The trailing comma is used to tell python that this is a tuple:
					dropTableStatement = "DROP TABLE IF EXISTS sCheck"
					cursorObject.execute(dropTableStatement)
					#GCF_002216145.1	WP_089035163.1	G19#DUF4065
					createTable = "CREATE TABLE sCheck (NCBI_Accnr text, Target text, Query text) ;"
					cursorObject.execute(createTable)
					#cur.execute("INSERT INTO stocks VALUES ('2006-01-05','BUY','RHAT',100,35.14)")
					with open(homeDir+'crossCheckList.txt', 'r') as cIn:
						for line in cIn: #GCF_904420025.1	WP_191552224.1	DUF4065
							Line=line.rstrip().split('\t')
							#if length_from_gwp(Line[0],Line[1])<=400: #LengthLimit400AA
							addLine='INSERT INTO sCheck VALUES ('+'"'+Line[0]+'"'+','+'"'+Line[1]+'"'+','+'"'+Line[2]+'"'+');'
							#print(addLine)
							cursorObject.execute(addLine)
							crossCheckSet.add(Line[1])

					conn.commit()

					def crossDBchecker(accession):
						conn = sqlite3.connect(sqlDB) # Make a connection object
						cursorObject = conn.cursor() # Create a cursor object
						# Prepare SQL query to retrieve a record from the database.
						bestHit = 'SELECT * FROM sCheck WHERE Target LIKE "%' + accession +'%"'
						# Execute the SQL command
						cursorObject.execute(bestHit)
						# Fetch all the rows in a list of lists.
						results = cursorObject.fetchall()
						queryList=set()
						for row in results:
							if row[1]==accession:
								queryList.add(row[2])
						cursorObject.close()
						if queryList and len(queryList)>0:
							return queryList

					if crossCheckSet:
						ctime = time.asctime( time.localtime(time.time()) )
						print("-- Crosscheck DB updated at ", ctime, ' and size : ', len(crossCheckSet))

					'''
					with open(homeDir+'crossCheckList.txt', 'r') as cIn:
						for line in cIn: #GCF_904420025.1	WP_191552224.1	DUF4065
							Line=line.rstrip().split('\t')[1]
							Line2=line.rstrip().split('\t')[1]+'|'+line.rstrip().split('\t')[2]
							crossCheckSet.add(Line) #WP_191552224.1
							crosscheckNameList.append(Line2) #WP_191552224.1|DUF4065

					if crossCheckSet:
						ctime = time.asctime( time.localtime(time.time()) )
						print("-- Crosscheck DB updated at ", ctime, ' and size : ', len(crossCheckSet))
						for item in crossCheckSet: #WP_191552224.1
							itemsSet=set()
							for items in crosscheckNameList: #WP_191552224.1|DUF4065
								if item==items.split('|')[0]:
									itemsSet.add(items.split('|')[1])
							crosscheckAccGroupDict[item]=itemsSet
					'''
					print('-- Cross Checking: '+str(cOutDir)+'/'+str(len(hmmOutDirList)))
					#clusFileExtensionName=str(cwd).split('/')[-1].split('_')[1]
					#netFlax_TA=outDir+'_netOut_*peronTA.txt'
					#netFlax_TAtsv=outDir+'_netOut_*peronTA.tsv'
					#netFlax_taxaList=outDir+'_netOut_taxaList.txt'
					#netFlax_outdesc=outDir+'_netOut_flankgene.fasta_cluster_out_*_outdesc.txt'
					#getClusScript=homeDir+localS+'getTAFam.py'
					#topTAcommand="python3 %s -i %s -d %s -t %s -l %s -o %s " %(getClusScript, netFlax_TA, netFlax_outdesc, int(spthresh), homeDir+outDir+'/'+netFlax_taxaList, clusFileExtensionName) #threshold -t '2' needs to be 10
					###topTAcommand="python3 ~/chayanBioPhD/netFlaX_AutoMation/getTopTAFam.py -i %s -d %s -t %s -n %s -o %s " %(netFlax_TA, netFlax_outdesc, int(spthresh), int(number), clusFileExtensionName) #threshold -t '2' needs to be 10
					###print(topTAcommand)
					#os.system(topTAcommand) ##need a cross check if clus not in glob
					if (glob.glob("*clus")):
						clusDict={}
						for clus in (glob.glob("*clus")):
							clusIdList=[]
							with open (clus,'r') as clusIn:
								for line in clusIn:
									Line=line.rstrip()
									clusIdList.append(Line)
							clusIdList_set=set(clusIdList)
							#print(clusIdList_set & crossCheckSet)
							clusDict[clus]=clusIdList
						#{'G2R1#DUF4065.clus': ['WP_074138355.1', 'WP_095295403.1']}
						#print(clusDict)
						for ids in clusDict:
							print('-- Cluster FileName: ', ids, '\n')
							if len(set(clusDict[ids]) & crossCheckSet)>0:
								pIdsset=set()
								for item in set(clusDict[ids]) & crossCheckSet:
									if crossDBchecker(item):
										for pIds in crossDBchecker(item):
											pIdsset.add(pIds)
								if len(pIdsset)>0:
									for prevIds in pIdsset:
										print('\n'+'-- DeadEnd:AccFoundInPrevRound'+'\t'+str(i)+'\t'+prevIds+'\t'+ids[:-5]+'\n')
										TAnetRes.add('DeadEnd:AccFoundInPrevRound'+'\t'+str(i)+'\t'+prevIds+'\t'+ids[:-5])
							else:
								cwd=os.getcwd()
								#print('Creating HMM for', cwd)
								alnHmmScript=homeDir+localS+'alnHmm.py'
								alnHmmCommand="python3 %s -i %s -c %s -t %s -lp %s" %(alnHmmScript, ids, core, netFlax_TAtsv, homeDir+localPG)
								#print(alnHmmCommand)
								#alnHmmCommand="python3 ~/chayanBioPhD/netFlaX_AutoMation/alnHmm.py -i %s -c %s -t %s" %(ids, core, netFlax_TAtsv)
								os.system(alnHmmCommand)
								clusGroup=ids.split('.')[0]
								if (glob.glob(clusGroup+".hmm")):
									hmmName=clusGroup+".hmm"
									print('\n'+'-- predictedTAgroup:Found'+'\t'+str(i)+'\t'+getPrevFromHMM(hmmName)+'\t'+clusGroup+'\t'+'FilterPassed')
									TAnetRes.add('\n'+'predictedTAgroup:Found'+'\t'+str(i)+'\t'+getPrevFromHMM(hmmName)+'\t'+clusGroup+'\t'+'FilterPassed')
									copyHMMs="cp %s %s" %(clusGroup+".hmm", '../'+localH)
									os.system(copyHMMs)
								else:
									if (glob.glob(clusGroup+".ccout")):
										with open (clusGroup+".ccout",'r') as ccOIn:
											for line in ccOIn:
												Line=line.rstrip().split('\t')
												print('\n'+'-- DeadEnd:HmmHitFoundInPrevRound'+'\t'+str(i)+'\t'+Line[1].split('.')[0]+'\t'+Line[2]+'\t'+str(Line[0]))
												TAnetRes.add('DeadEnd:HmmHitFoundInPrevRound'+'\t'+str(i)+'\t'+Line[1].split('.')[0]+'\t'+Line[2]+'\t'+str(Line[0]))
							#crossCheckSet=set()
						#crosscheckAccGroupDict={}
						#crosscheckNameList=[]
						cwdback='/'.join(map(str,cwd.split('/')[:-1]))
						#print(cwdback)
						os.chdir(cwdback)
					else:
						#crossCheckSet=set()
						#crosscheckAccGroupDict={}
						#crosscheckNameList=[]
						print('\n'+'DeadEndFound')
						cwdback='/'.join(map(str,cwd.split('/')[:-1]))
						#print(cwdback)
						os.chdir(cwdback)
					#print('## Round '+ str(i) +' Finished.')
				localtime2 = time.asctime( time.localtime(time.time()) )
				print("-- Local current time after crosscheck done for ", str(cOutDir)+'/'+str(len(hmmOutDirList)), ' : ', localtime2)
				print("-- Crosscheck DB size : ", len(crossCheckSet))
				for item in TAnetRes:
					print(item, file=deadOut)
					time.sleep(1)
			else:
				print('not working'+'\n')
				sys.exit()
	else:
		print('#######Hopping Finished#######'+'\n')
		sys.exit()
#for hmms in hmms/*hmm;	 do		 mkdir ${hmms##hmms/}Out;		 for fasta in *faa;			 do hmmsearch --domtblout ${hmms##hmms/}Out/${fasta%%.faa}_${hmms##hmms/}_tab -A ${hmms##hmms/}Out/${fasta%%.faa}_${hmms##hmms/}_aln --incE 1e-10 --incdomE 1e-10 -E 1e-10 --cpu 4 $hmms $fasta > ${hmms##hmms/}Out/${fasta%%.faa}_${hmms##hmms/}_out;		 done;	 done


'''
source  link    target  colorestim      family  colorTAdb       Desc    TAId
DUF4065 pp      G28R1   pTox    -/RHH-MazF|T2   Tox      hypothetical protein #7        6252
DUF4065 pp      G7R1    pTox    G7R1    none     hypothetical protein #40
DUF4065 pp      G21R1   pTox    G21R1   none     hypothetical protein #13
DUF4065 pp      G13R1   pTox    -/(toxin)|T2    Tox      hypothetical protein #10       6093
DUF4065 pp      G6R1    pTox    G6R1    none     hypothetical protein #47
DUF4065 pp      G9R1    pTox    G9R1    none     hypothetical protein #25
DUF4065 pp      G8R1    pTox    G8R1    none     hypothetical protein #37
DUF4065 pp      G2R1    pTox    G2R1    none     hypothetical protein #61
DUF4065 pp      G12R1   pTox    G12R1   none     hypothetical protein #20

WP_102375551.1	155
WP_106007974.1	105
WP_123549286.1	119
WP_021168410.1	45
WP_087880394.1	18

'''
