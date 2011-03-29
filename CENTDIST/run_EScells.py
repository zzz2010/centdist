import os
from glob import glob
from random import randint
srcfolder='../experiments/chipseq'
destfolder='EScells'
genomedir='/data2/genome/M_musculus_UCSC_mm8/orig'
while 1:
	left=list(set(map(lambda x:x.rsplit('/',1)[-1],glob('%s/loci*'%srcfolder)))-set(map(lambda x:x.rsplit('/',1)[-1],glob('%s/loci*'%destfolder))))
	if len(left)==0:
		break
	f=left[randint(0,len(left)-1)]
	os.system("python pomoscan.py -i %s/%s -w 500 -o %s/%s -m transfac.motifs -g %s"%(srcfolder,f,destfolder,f,genomedir))
	#os.system('python processcan.py %s/%s 3000 %s/%s %s'%(srcfolder,f,destfolder,f,genomedir))
	#os.system('rm %s/%s/*.pos'%(destfolder,f))
	#os.system('rm %s/%s/*.fas'%(destfolder,f))

