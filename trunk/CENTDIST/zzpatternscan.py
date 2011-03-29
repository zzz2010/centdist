import sys
import os

if len(sys.argv)<=1:
        print "arguments: pwmfile [genomedir peakfile w|fastafile] outputdir fp(specify 1 for 1 in 10000)"
	sys.exit(1)

from time import time
from glob import *
from common import *

#patternScan=sys.path[0]+'/patternScan.zzz.2010.11.26.exe'
patternScan=sys.path[0]+'/patternScan.zzz.2010.12.01.exe'
Rwrapper=format('python {sys.path[0]}/Rwrapper.py',locals())
plotgraphsR=format('{Rwrapper} {sys.path[0]}/plotgraphs.R',locals())
patternDistr=format('python {sys.path[0]}/Pos2Bed.py',locals())
extractfas=format('python {sys.path[0]}/../extractfas.py',locals())
plotmovingaverage=format('{Rwrapper} {sys.path[0]}/plotmovingaverage.R',locals())

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
w=int(args.get('w',[0])[0])
outputdir=args.get('outputdir')[0]
randbgfas=sys.path[0]+'/randbg.fas'
fp=float(args.get('fp')[0])
movingaveragebinsize=int(args.get('movingaveragebinsize',[1000])[0])
if fastafile==None and (peakfile==None and genomedir==None or w==0):
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
info['fp']=[fp]
with open('{outputdir}/info.txt'.format(outputdir=outputdir),'w') as infofile:
	for k,v in info.items():
		print>>infofile,k+'\t'+'\t'.join(map(str,v))

start=time()
#command=format("{patternScan} -i '{pwmfile}' -o '{outputdir}' -f '{fastafile}' -fp {fp} -bg '{randbgfas}'",locals())
command=format("{patternScan} -i '{pwmfile}' -o '{outputdir}' -f '{fastafile}' -fp {fp}",locals())
if (ossystem(command)):
	raise Exception("patternScan failed")

if ossystem(format("ls {outputdir}/*pos|xargs -n1 sh -c '{patternDistr} {fastafile} {peakfile} $0 ' ",locals())):
	  raise Exception("Pos2Bed failed")

command=format("ls {outputdir}/*pos|xargs -n1 sh -c '{plotmovingaverage} $0 {N} {movingaveragebinsize}'",locals())
if (ossystem(command)):
	raise Exception("patternScan failed")
	
