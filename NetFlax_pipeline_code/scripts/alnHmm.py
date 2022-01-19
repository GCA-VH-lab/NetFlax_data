#Author: Chayan Kumar Saha
import argparse
import glob, gzip
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC
import os, sys, os.path

usage= '''  Description: Get .aln and .HMM file '''


parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--input", required=True, help="TopTA_Cluster.clus")
parser.add_argument("-t", "--tsv", required=True, help="_netOut_TreeOrderOperonTA.tsv")
parser.add_argument("-c", "--cpu", help="Maximum number of parallel CPU workers to use for multithreads. ")
parser.add_argument("-lp", "--localProteomeDirectory", help=" Path for Local Proteome.faa Files, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-v","--version", action="version", version='%(prog)s 2.0.21')
args = parser.parse_args()
parser.parse_args()

Entrez.email = 'chayan.sust7@gmail.com'
Entrez.api_key= 'eab0a98c3fa1495bc9532a9c9b137aa1e308'

if args.localProteomeDirectory:
	if os.path.isdir(args.localProteomeDirectory):
		if args.localProteomeDirectory[-1]=='/':
			localP=args.localProteomeDirectory
			#print('Local Proteome.faa Data path : ', localP, '\n')
		else:
			localPG=args.localProteomeDirectory+'/'
			#print('Local Proteome.faa Data path : ', localP, '\n')
	else:
		print('No directory Found as : '+ args.localProteomeDirectory)
		sys.exit()
else:
	localP='./'

homeDir='/'.join(map(str,localP.split('/')[:-2]))+'/'


if args.cpu:
	if int(args.cpu)>0:
		core=int(args.cpu)
	else:
		print('Please use number eg, 1,2...')
		sys.exit()

def seq_from_acc(accession_nr):
	"""
	:param accession_nr: NCBI protein accession
	:return: Protein Sequence
	"""
	try:
		handle = Entrez.efetch(db="protein", id=accession_nr, rettype="gbwithparts", retmode="text")
	except Exception as e:
		return False
	record = SeqIO.read(handle, "genbank")
	handle.close()
	return record.format("fasta")

tsvDict={}
#WP_006165776.1#236_Natrialba_chahannaoensis     318     +       +       129     -2795   -2478   35775   36092   WP_049896408.1#236      NZ_AOIN01000017.1       GCF_000337135.1
with open(args.tsv, 'r') as tsvIn:
	for line in tsvIn:
		if '#' in line:
			Line=line.rstrip().split('\t')
			tsvDict[Line[9].split('#')[0]]=Line[11]

#clusNum

idSet=set()
with open(args.input, 'r') as fileIn:
	for line in fileIn:
		Line=line.rstrip()
		idSet.add(Line)

taTXTlist=[]
with open(args.tsv[:-4]+'.txt', 'r') as txtIn:
	for line  in txtIn:
		if line[0]!='#':
			Line=line.rstrip().split('\t')
			taTXTlist.append(Line)


def seq_from_wp_G(gcf_nr,accession_nr):
	faaFile=gcf_nr+'.faa.gz'
	LocalP=localP
	fastaSeq = gzip.open(LocalP+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==accession_nr:
			#record.id = gcf_nr+'\t'+accession_nr
			#record.description = str('')
			#handle.close()
			return record.format("fasta")

def seq_from_wp_stat(gcf_nr,accession_nr):
	faaFile=gcf_nr+'.faa.gz'
	LocalP=localP
	fastaSeq = gzip.open(LocalP+faaFile, "rt")
	for record in SeqIO.parse(fastaSeq, "fasta"):
		if record.id==accession_nr:
			record.id = gcf_nr+'\t'+accession_nr
			return record.id+'\t'+str(len(record.seq))

if len(idSet)>0:
	with open(args.input.split('.')[0]+'.fas', 'w') as fasOut, open(args.input.split('.')[0]+'.stat', 'w') as statOut:
		for item in idSet:
			if item in tsvDict:
				print(seq_from_wp_G(tsvDict[item], item), file=fasOut)
				print(seq_from_wp_stat(tsvDict[item], item), file=statOut)



	fasFile=args.input.split('.')[0]+'.fas'

	HmmTrackList=args.input.split('.')[0].split('#')[1:]
	HMMListTobeChecked=[]
	for h in range(1,len(HmmTrackList)+1):
		roundHmms=homeDir+'Round_'+str(h)+'_HMMs/'
		#roundHmms='/data/users/chayan/chayanBioPhD/netFlaX_AutoMation/Round_'+str(h)+'_HMMs/'
		for hmms in (glob.glob(roundHmms+'*hmm')):
			#print(hmms)
			HMMListTobeChecked.append(hmms)
	recRound=homeDir+'HMMs/'
	for recentHmms in (glob.glob(recRound+'*hmm')):
		if recentHmms:
			HMMListTobeChecked.append(recentHmms)
	#print(HmmTrackList)
	print('## Checking: ', len(HMMListTobeChecked), 'HMMs')

	HMMcommandList=[]
	for items in HMMListTobeChecked:
		HMMs=items.split('/')[-1]
		#print(HMMs)
		if HMMs[0]=='D':
			#print(HMMs)
			command="A=`hmmsearch --cut_ga --cpu %s %s %s | grep '^Domain search space  (domZ):' | tr ':' '\t' | cut -f 2 | tr -d ' ' | cut -d '[' -f 1`; B=`echo %s`; C=`echo %s`; echo $A $B $C >> %s" %(core, items, fasFile, HMMs, fasFile.split('.')[0], fasFile.split('.')[0]+'.cc')
			#command="A=`hmmsearch --cut_ga --cpu %s %s %s | grep '^Domain search space  (domZ):' | tr ':' '\t' | cut -f 2 | tr -d ' ' | cut -d '[' -f 1`; B=`hmmsearch --cut_ga --cpu %s %s %s | grep '^Query:'`; echo $A $B >> %s" %(core, items, fasFile, core, items, fasFile, fasFile.split('.')[0]+'.cc')
			HMMcommandList.append(command)
		if HMMs[0]=='G':
			#print(HMMs)
			command="A=`hmmsearch --incE 1e-10 --incdomE 1e-10 -E 1e-10 --cpu %s %s %s | grep '^Domain search space  (domZ):' | tr ':' '\t' | cut -f 2 | tr -d ' ' | cut -d '[' -f 1`; B=`echo %s`; C=`echo %s`; echo $A $B $C >> %s" %(core, items, fasFile, HMMs, fasFile.split('.')[0], fasFile.split('.')[0]+'.cc')
			#command="A=`hmmsearch --incE 1e-10 --incdomE 1e-10 -E 1e-10 --cpu %s %s %s | grep '^Domain search space  (domZ):' | tr ':' '\t' | cut -f 2 | tr -d ' ' | cut -d '[' -f 1`; B=`hmmsearch --incE 1e-10 --incdomE 1e-10 -E 1e-10 --cpu %s %s %s | grep '^Query:'`; echo $A $B >> %s" %(core, items, fasFile, core, items, fasFile, fasFile.split('.')[0]+'.cc')
			HMMcommandList.append(command)

	#print(HMMcommandList)


	for commands in HMMcommandList:
		#print(commands)
		os.system(commands)

	HmmCheck=[]
	with open(fasFile.split('.')[0]+'.cc', 'r') as ccIn:
		for line in ccIn:
			line=line.rstrip().split(' ')
			if int(line[0])>0:
				HmmCheck.append(line)

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
				return len(discardHmm)
		else:
			return 0 #no cross verified
 #Multiple cross with prev hmm


	#print(HmmCheck)
	if len(HmmCheck)==0 or XConFind(HmmCheck)==0:
		alnFile=args.input.split('.')[0]+'.aln'
		HMMfile=args.input.split('.')[0]+'.hmm'
		command1="mafft --maxiterate 1000 --localpair --thread %s --quiet %s > %s" %(core, fasFile, alnFile)
		command2="hmmbuild %s %s" %(HMMfile, alnFile)
		#print(command1)
		os.system(command1)
		#print(command2)
		os.system(command2)
		if (glob.glob(HMMfile)):
			with open (HMMfile[:-4]+'.clusNum', 'r') as cNum, open (HMMfile[:-4]+'.TAinfo', 'w') as TAinfotxt: #WP_102375551.1	155
				for line in cNum:
					Line=line.rstrip().split('\t')
					for items in taTXTlist:
						#print(Line, items)
						if Line[0] in items[3] and int(items[0].split('#')[1].split('|')[0])==int(Line[1]):
							print('\t'.join(map(str,items)), file=TAinfotxt)
	else:
		if XConFind(HmmCheck)>=1:
			with open(fasFile.split('.')[0]+'.ccout', 'w') as ccXOut:
				for items in HmmCheck:
					print(items[0],items[1],items[2], sep='\t', file= ccXOut)
else:
	with open (args.input.split('.')[0]+'.cc', 'w') as ccOut:
		print('DeadEndFound : Filtering Criteria', file=ccOut)
