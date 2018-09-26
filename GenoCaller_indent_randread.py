#!/usr/bin/env python
# -*- coding: ASCII -*-

###This program calls genotypes from bam files while taking into account post mortem damage in UDG-half treated file by ignoring the ends of reads, and only sampling one read randomly.
###One file is created:
###A emit all vcf file noting the call for all base pairs in the bed file
##GenoCaller_indent_randread.py <indexed bamfile> <bed file> <reference genome> <indent>


from sys import argv
import pysam
import math
import numpy as np
import string
import numpy.ma as ma
import random
import matplotlib
matplotlib.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import fmin
from scipy.optimize import fminbound
from scipy.optimize import curve_fit
import time
from random import randint



filenamein=argv[1]
filenameinB=argv[2]
ref_file=argv[3]
indent=int(argv[4])
min_RD=1
MQ=30 
BQ=30



if filenamein[-4:] == '.bam':
    filenameout=string.split(filenamein,'.bam')[0]
else:
    filenameout=filenamein[:]

    
if '/' in filenamein:
    filenameout=string.split(filenameout,'/')[-1]

filenameout1=filenameout+'.'+filenameinB+'.indent'+str(indent)+'.singleread.emit_all.vcf'


def phred2prob(x):
    return 10.0**(-x/10.0)

def prob2phred(x):
    return -10*math.log10(x)




###open up reference file
ref=pysam.FastaFile(ref_file)


###make a list of regions to be interogated
file = open(filenameinB)
data=file.read()
data=string.split(data,'\n')
file.close()

if data[-1]=='':
    del(data[-1])

SNPdir={}
SNPll={}
SNPlist=[]

count=0

for g in range(len(data)):
    k=string.split(data[g])
    key=k[0]+'_'+k[1]+'_'+k[2]
##    chromo=k[0]
##    start=int(k[1])
##    end=int(k[2])
##    seq=ref.fetch(chromo,start,end)
##    seq=seq.upper()
    SNPlist.append(key)
    SNPdir[key]=k

###open up bam file (must be indexed)
samfile = pysam.AlignmentFile(filenamein, "rb")
samp_name=samfile.header['RG'][0]['SM']
count=0


###set up output files
fileout=open(filenameout1,'w')

###set up output headers
out1='##fileformat=VCFv4.1\n'
outA='##Caller_arguments=<MappingQuality_filter='+str(MQ)+',BaseQuality_filter='+str(BQ)
outA=outA+',minRD='+str(min_RD)
outA=outA+',bam_in='+filenamein+',bedfile='+filenameinB
outA=outA+',indent='+str(indent)+'>\n'
outB='##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n'
outB=outB+'##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n'

outC='##Time created='+(time.strftime("%H:%M:%S"))+' '+(time.strftime("%d/%m/%Y"))+'\n'
for i in range(len(ref.lengths)):
    outC=outC+'##contig=<ID='+ref.references[i]+',length='+str(ref.lengths[i])+'>\n'
outC=outC+'##reference=file:'+ref_file+'\n'

out_e='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+samp_name+'\n'

fileout.write(out1+outA+outB+outC+out_e)


###Start running through regions in the bed
for gg in range(len(SNPlist)):
    #print 'Calculating genotype likelihoods for locus '+str(gg+1)+' : '+SNPlist[gg]
    chromo=SNPdir[SNPlist[gg]][0]
    pos_start=int(SNPdir[SNPlist[gg]][1])
    pos_end=int(SNPdir[SNPlist[gg]][2])
    use=0
    ###Go base by base through each region in the bam file
    for pileupcolumn in samfile.pileup(chromo,pos_start,pos_end,truncate=True,stepper='all'):
        nucl=ref.fetch(chromo,pileupcolumn.pos,pileupcolumn.pos+1)  ##grab reference sequence for region
        nucl=nucl.upper()
        var_list=[]  ##a list to store read bases
        index=pileupcolumn.pos-pos_start  ##tells us what position we are in for the arrays given the nucleotide position
        ###go read by read for each position
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:  ##ensure not an indel or duplicate
                #map_list.append(pileupread.alignment.mapping_quality)  ##record mapping quality of read
                if (pileupread.alignment.mapping_quality>=MQ) and (ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ):  ###Add base calls meeting the MQ and BQ filters
                    if (pileupread.query_position>indent) and ((len(pileupread.alignment.query_sequence)-1-pileupread.query_position)>indent): ###check base is not at end of read
                        var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33])

        if chromo[:3]=='chr':
            chromo_out=chromo[3:]
        else:
            chromo_out=chromo
        out=chromo_out+'\t'+str(pileupcolumn.pos+1)+'\t.\t'+nucl

        if len(var_list)>0:
            random.shuffle(var_list)
            call=var_list[0][0]
            if call<>nucl:
                out=out+'\t'+call+'\t99\tPASS\tAC=1;AF=1\tGT:AD:DP:GQ:PL\t1/1:0,1:1:99:255,255,0\n'
            else:
                out=out+'\t.\t99\tPASS\tAC=0;AF=0\tGT:DP:GQ:PL\t0/0:1:99:0,255,255\n'
        else:
            out=out+'\t.\t.\t.\t.\t.\t./.\n'
            
        fileout.write(out)
        use+=1

    if use==0:
        if chromo[:3]=='chr':
            chromo_out=chromo[3:]
        else:
            chromo_out=chromo
            
        nucl=ref.fetch(chromo,pos_start,pos_end)
        nucl=nucl.upper()
        out=chromo_out+'\t'+str(pos_end)+'\t.\t'+nucl+'\t.\t.\t.\t.\t.\t./.\n'
        fileout.write(out)

    if gg%10000==0:
        print SNPlist[gg]
        
fileout.close()

