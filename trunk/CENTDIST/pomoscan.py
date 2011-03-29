import sys
import os
from time import time
from glob import *
grandstart=time()
params=dict()
params['genomedir']=None
params['matfile']="transfac.motifs"
params['inputfile']=None
params['w']=None
params['outputdir']=None
params['fastafile']=None
params['progdir']='.'
info=dict()
ERROR=0

try:
	i=1
	while i<len(sys.argv):
		if sys.argv[i]=='-f':
			i+=1
			params['fastafile']=sys.argv[i]
		elif sys.argv[i]=='-d':
			i+=1
			params['progdir']=sys.argv[i]
		elif sys.argv[i]=='-g':
			i+=1
			params['genomedir']=sys.argv[i]
		elif sys.argv[i]=='-i':
			i+=1
			params['inputfile']=sys.argv[i]
		elif sys.argv[i]=='-w':
			i+=1
			params['w']=int(sys.argv[i])
		elif sys.argv[i]=='-o':
			i+=1
			params['outputdir']=sys.argv[i]
		elif sys.argv[i]=='-m':
			i+=1
			params['matfile']=sys.argv[i]
		else:
			print "invalid parameter: "%sys.argv[i]
		i+=1
except Exception,e:
	print e
	ERROR=1
randbgfas='%s/randbg.fas'%params['progdir']

if params['fastafile']==None and (params['inputfile']==None or params['genomedir']==None or params['w']==None):
	ERROR=1

if params['outputdir']==None or params['matfile']==None:
	ERROR=1

if ERROR:
	print("ERROR: some parameters not set")
	print("arguments: -i peakfile  -g genomedir -w maxwindowsize -o outputdir -m motiffile")
	print("or -f fastafile -o outputdir -m motiffile")
	sys.exit(1)

for doonce in [1]:
	if glob(params['outputdir'])==[]:
		try:
			s=("mkdir '%s'"%params['outputdir'])
			print(s)
			if (os.system(s)):
				1/0
			s=("chmod 777 '%s'"%params['outputdir'])
			print(s)
			if (os.system(s)):
				1/0
		except:
			print("error creating output directory")

	if (params['fastafile']==None):
		s=("python extractfas.py '%s' '%s' '%s' > '%s/fasfile.fas'"%(params['genomedir'],params['inputfile'],6*params['w']+100,params['outputdir']))
		if (os.system(s)):
			print("ERROR OCCURRED")
			break
		fastafile='%s/fasfile.fas'%params['outputdir']
	else:
		fastafile=params['fastafile']

	params['N']=len(filter(lambda x:'>' in x,open(fastafile).readlines()))
	L=0
	first=1
	for line in open(fastafile):
		if first:
			first=0
		else:
			LL=len(line.strip())
			if LL and line[0]<>'>':
				L+=LL
			elif line[0]=='>':
				break
	params['L']=L

	start=time()
	#s=("%s/patternScan -i '%s' -o '%s' -f '%s' -fp '1' -bg %s"%(params['progdir'],params['matfile'],params['outputdir'],fastafile,randbgfas))
	s=("%s/patternScan -i '%s' -o '%s' -f '%s' -fp '1' "%(params['progdir'],params['matfile'],params['outputdir'],fastafile))
	print(s)
	if (os.system(s)):
		print("ERROR OCCURRED")
		break
	info['scan_time']=time()-start

for doonce in [1]:
	tempfile=params['outputdir']+'/temp.txt'
	os.system("rm %s >/dev/null 2>&1"%tempfile)
	for posfile in glob(params['outputdir']+'/*.pos'):
		s="%s/patternDistr -i '%s' -findThresh  >> '%s'"%(params['progdir'],posfile,tempfile)
		print s
		if (os.system(s)):
			None
			#print("ERROR OCCURRED")
			# break
	resultfile=open(params['outputdir']+'/results.txt','w')
	print>>resultfile,'NAME\tSCORE\tW\tTHRESHOLD\tZ0SCORE\tZ1SCORE'
	for line in open(tempfile):
		print>>resultfile,line.strip().replace(params['outputdir']+'/','').replace('.pos','')
	resultfile.close()
	#os.system("rm %s >/dev/null 2>&1"%tempfile)

for doonce in [1]:
	print>>open('%s/info.txt'%params['outputdir'],"w"),'\n'.join(map(lambda x:'%s\t%s'%x,params.items())+map(lambda x:'%s\t%s'%x,info.items()))
	start=time()
	s=("R %s  --no-save --args --slave < %s/plotgraphsnew.R"%(params['outputdir'],params['progdir']))
	print(s)
	if (os.system(s)):
		print("ERROR OCCURRED")
		break
	info['drawing_time']=time()-start
	info['total_time']=time()-grandstart
	print>>open('%s/info.txt'%params['outputdir'],"a"),'\n'.join(map(lambda x:'%s\t%s'%x,[('drawing_time',info['drawing_time']),('total_time',info['total_time'])]))
	#os.system('rm -rf %s/*.pos %s/*.fas'%(params['outputdir'],params['outputdir']))

