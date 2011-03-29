import sys
import os
from time import time
from glob import *
if os.getenv('PYTHONPATH')=='':
	sys.path.append('/home/chipseq/public_html/webseqtools/COMMONPY')
from common import *

#patternScan=sys.path[0]+'/patternScan.zzz.2010.11.26.exe'
#patternScan=sys.path[0]+'/patternScan.zzz.2010.12.01.exe'
patternScan=sys.path[0]+'/patternScan.zzz.2010.12.01.2.exe'
Rwrapper=format('python {sys.path[0]}/Rwrapper.py',locals())
plotgraphsR=format('{Rwrapper} {sys.path[0]}/plotgraphs.R',locals())
patternDistr=format('{sys.path[0]}/patternDistr',locals())
#patternDistr=format('python {sys.path[0]}/patternDistr.py',locals()) 
extractfas=format('python {sys.path[0]}/../extractfas.py',locals())

grandstart=time()
args=getargs()
try:
	configfile=args['configfile'][0]
	args=readconfig(configfile,args)
except:
	pass
info=dict()

genomedir=args.get('genomedir',[None])[0]
pwmfile=args.get('pwmfile')[0]
peakfile=args.get('peakfile',[None])[0]
fastafile=args.get('fastafile',[None])[0]
FP=float(args.get('FP',[0.0001])[0])*10000
w=int(args.get('w',[0])[0])
outputdir=args.get('outputdir')[0]
randbgfas=sys.path[0]+'/randbg.fas'

if fastafile==None and (peakfile==None and genomedir==None or w==None):
	raise Exception("fastafile or peakfile and genomedir and w required")

if glob(outputdir)==[]:
	command=format("mkdir {outputdir}",locals())
	if (ossystem(command)):
		raise Exception("mkdir failed")

if fastafile==None:
	seqlen=2*w
	fastafile=format('{outputdir}/fastafile.fa',locals())
	command=format("{extractfas} -genomedir {genomedir} -peakfile {peakfile} -w {seqlen} > {fastafile}",locals())
	if (ossystem(command)):
		raise Exception("extractfas failed")
N=len(filter(lambda x:x[0]=='>',open(fastafile)))
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

info['genomedir']=[genomedir]
info['pwmfile']=[pwmfile]
info['peakfile']=[peakfile]
info['fastafile']=[fastafile]
info['w']=[w]
info['outputdir']=[outputdir]
info['randbgfas']=[randbgfas]
info['N']=[N]
info['L']=[L]

start=time()
command=format("{patternScan} -i '{pwmfile}' -o '{outputdir}' -f '{fastafile}' -fp '{FP}' -one ",locals())
if (ossystem(command)):
	raise Exception("patternScan failed")
info['scan_time']=[time()-start]
resultfile=format("{outputdir}/results.txt",locals())
print >> open(format("{resultfile}",locals()),'w'),'\t'.join("NAME    SCORE   W       THRESHOLD       Z0SCORE Z1SCORE".split())
for posfile in glob(format('{outputdir}/*.pos',locals())):
	if ossystem(format("{patternDistr} -i '{posfile}' -findThresh >> '{resultfile}'",locals())):
	#if ossystem(format("{patternDistr} '{fastafile}' '{peakfile}' '{posfile}'  >> '{resultfile}'",locals())):
		#raise Exception("patternDistr failed")
		print "patternDistr failed"
		

start=time()
writeconfig(format("{outputdir}/info.txt",locals()),info)
if ossystem(format("{plotgraphsR} {outputdir}",locals())):
	raise Exception("plotgraphs failed")

info['draw_time']=[time()-start]
writeconfig(format("{outputdir}/info.txt",locals()),info)

ossystem(format('rm -f {outputdir}/*.pos {fastafile}',locals()))
