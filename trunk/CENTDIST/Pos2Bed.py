#!/usr/bin/python
import os
import sys
import commands

fastafile=sys.argv[1]
peakfile=sys.argv[2]
posfile=sys.argv[3]
extlen=500
motifbed=posfile+'.bed'
peakbed=peakfile+'.bed'
os.system('dos2unix '+peakfile)
os.system('awk \'BEGIN{OFS=\"\\t\"}{print $1,$2,$2+1}\' '+peakfile+' > '+peakbed)
motifbed=motifbed.replace('$','_')
##########################prepare motifbed file##########################################
lines=open(fastafile,'r').readlines()
seqname=list()
lastline='0'
row=0
for line in lines:
        row+=1
        if line[0]!='>' and lastline[0]=='>':
                comps=lastline[1:].split()
                seqname.append(comps[0])
                #seqname[row-1]=(comps[0])
        lastline=line
lastpos=('chr0',0,0)
lines=open(posfile,'r').readlines()
outf=open(motifbed,'w')
for line in lines:
        comps=line.split()
        seqid=int(comps[0])-1
        comps2=seqname[seqid].split(':')
        chr=comps2[0]
        comps3=comps2[1].split('-')
        pos1=int(comps3[0])+int(comps[1])
        pos2=int(comps3[0])+int(comps[1])+10 #motif len is 10 assume
#        if chr==lastpos[0] and (pos2-lastpos[1])*(pos2-lastpos[2])<0:
#                lastpos=(chr,pos1,pos2)
#                continue
#        else :
        outf.write(chr+'\t'+str(pos1)+'\t'+str(pos2)+'\t'+comps[2]+'\n')
        lastpos=(chr,pos1,pos2)
outf.close()
