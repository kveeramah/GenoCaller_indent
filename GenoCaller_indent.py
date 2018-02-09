#!/usr/bin/env python
# -*- coding: ASCII -*-

###This program calls genotypes from bam files while taking into account post mortem damage in UDG-half treated file by ignoring the ends of reads.
###Otherwise the algorithm is the same as GATK Unified Genotype for diploid calls
###Three files are created:
###An emit all vcf file noting the call for all base pairs in the bed file
###This version of the program will take any genotype with a low quality heterozygote call (default GQ=30) and convert to the next best homozygote call
###A vcf file with only those sites with evidence for at least one alternative allele and that passes a predetermined QUAL filter 
###A haploid emit all vcf, that gives you the most likely base under a haploid model. If two or more basepairs are tied, the reported allele is randomly chosen.
###to run type:
##GenoCaller_indent.py <indexed bamfile> <bed file> <reference genome> <indent>


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
MQ=15
BQ=15
GQ=30
QUALpass=30
theta=0.001


if filenamein[-4:] == '.bam':
    filenameout=string.split(filenamein,'.bam')[0]
else:
    filenameout=filenamein[:]

    
if '/' in filenamein:
    filenameout=string.split(filenameout,'/')[-1]

filenameout1=filenameout+'.'+filenameinB+'.indent'+str(indent)+'.emit_all.vcf'
filenameout2=filenameout+'.'+filenameinB+'.indent'+str(indent)+'.vcf'
filenameout3=filenameout+'.'+filenameinB+'.indent'+str(indent)+'.haploid.emit_all.vcf_like'



def phred2prob(x):
    return 10.0**(-x/10.0)

def prob2phred(x):
    return -10*math.log10(x)

#genotype called that incorporates damage for aDNA
def geno_caller_10GT_aDNA(X):
    GL=np.zeros(10)   #all 10 possible genotypes and order = AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
    hap=np.zeros((len(X),4))  #all 4 haploid possibilities, A,C,G,T
    all_dic={}
    all_dic['A']=0
    all_dic['C']=1
    all_dic['G']=2
    all_dic['T']=3
    count=0
    
    for g in range(len(X)):
        if all_dic.has_key(X[g][0])==False:
            continue
        err=phred2prob(X[g][1])
        hap[g]=err/3.0
        if X[g][0]=='A':
            hap[g][0]=1-err
            hap[g][2]=((1-X[g][3])*(err/3.0))+(X[g][3]*(1-err))
        elif X[g][0]=='C':
            hap[g][1]= ((1-X[g][2])*(1-err))+(X[g][2]*(err/3.0))
        elif X[g][0]=='G':
             hap[g][2]= ((1-X[g][3])*(1-err))+(X[g][3]*(err/3.0))           
        elif X[g][0]=='T':
            hap[g][3]=1-err
            hap[g][1]=((1-X[g][2])*(err/3.0))+(X[g][2]*(1-err))                   

        GL[0]=GL[0]+math.log10(hap[g][0])                      
        GL[1]=GL[1]+math.log10((hap[g][0]+hap[g][1])/2)       
        GL[2]=GL[2]+math.log10((hap[g][0]+hap[g][2])/2)
        GL[3]=GL[3]+math.log10((hap[g][0]+hap[g][3])/2)
        GL[4]=GL[4]+math.log10(hap[g][1])
        GL[5]=GL[5]+math.log10((hap[g][1]+hap[g][2])/2)
        GL[6]=GL[6]+math.log10((hap[g][1]+hap[g][3])/2)
        GL[7]=GL[7]+math.log10(hap[g][2])
        GL[8]=GL[8]+math.log10((hap[g][2]+hap[g][3])/2)
        GL[9]=GL[9]+math.log10(hap[g][3])
        count+=1

    if count==0:
        GL.fill(-9)
    return GL


#standard genotype caller
def geno_caller_10GT(X):
  
    GL=np.zeros(10)   #all 10 possible genotypes and order = AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
    hap=np.zeros((len(X),4))  #all 4 haploid possibilities, A,C,G,T

    all_dic={}
    all_dic['A']=0
    all_dic['C']=1
    all_dic['G']=2
    all_dic['T']=3
    
    count=0
    for g in range(len(X)):
        if all_dic.has_key(X[g][0])==False:
            continue
        hap[g]=phred2prob(X[g][1])*(1.0/3.0)
        hap[g][all_dic[X[g][0]]]=1-phred2prob(X[g][1])

        GL[0]=GL[0]+math.log10(hap[g][0])
        GL[1]=GL[1]+math.log10((hap[g][0]+hap[g][1])/2)
        GL[2]=GL[2]+math.log10((hap[g][0]+hap[g][2])/2)
        GL[3]=GL[3]+math.log10((hap[g][0]+hap[g][3])/2)

        GL[4]=GL[4]+math.log10(hap[g][1])
        GL[5]=GL[5]+math.log10((hap[g][1]+hap[g][2])/2)
        GL[6]=GL[6]+math.log10((hap[g][1]+hap[g][3])/2)

        GL[7]=GL[7]+math.log10(hap[g][2])
        GL[8]=GL[8]+math.log10((hap[g][2]+hap[g][3])/2)

        GL[9]=GL[9]+math.log10(hap[g][3])
        count+=1

    if count==0:
        GL.fill(-9.0)
    return GL



###open up reference file
ref=pysam.FastaFile(ref_file)


###set up various look up dictionaries
all_dic={}
all_dic['A']=0
all_dic['C']=1
all_dic['G']=2
all_dic['T']=3
all_dic['N']=-99
all_dic[0]='A'
all_dic[1]='C'
all_dic[2]='G'
all_dic[3]='T'
all_dic[-9]='.'
all_dic[-99]='N'

GL_dic={}

GL_dic['AA']=0
GL_dic['AC']=1
GL_dic['AG']=2
GL_dic['AT']=3
GL_dic['CC']=4
GL_dic['CG']=5
GL_dic['CT']=6
GL_dic['GG']=7
GL_dic['GT']=8
GL_dic['TT']=9

GL_dic[0]='AA'
GL_dic[1]='AC'
GL_dic[2]='AG'
GL_dic[3]='AT'
GL_dic[4]='CC'
GL_dic[5]='CG'
GL_dic[6]='CT'
GL_dic[7]='GG'
GL_dic[8]='GT'
GL_dic[9]='TT'
GL_dic[-9]='./.'

homs=[0,4,7,9]

alt_dic={}
alt_dic[0]=[0,0]
alt_dic[1]=[0,1]
alt_dic[2]=[0,2]
alt_dic[3]=[0,3]
alt_dic[4]=[1,1]
alt_dic[5]=[1,2]
alt_dic[6]=[1,3]
alt_dic[7]=[2,2]
alt_dic[8]=[2,3]
alt_dic[9]=[3,3]

tri_dic={}
tri_dic[1]='A,C'
tri_dic[2]='A,G'
tri_dic[3]='A,T'
tri_dic[5]='C,G'
tri_dic[6]='C,T'
tri_dic[8]='G,T'


LL0_map={}
LL0_map[0]=[0]
LL0_map[1]=[4]
LL0_map[2]=[7]
LL0_map[3]=[9]
            
LL1_map={}
LL1_map[0]=[1,2,3]
LL1_map[1]=[1,5,6]
LL1_map[2]=[2,5,8]
LL1_map[3]=[3,6,8]
            
LL2_map={}
LL2_map[0]=[4,5,6,7,8,9]
LL2_map[1]=[0,2,3,7,8,9]
LL2_map[2]=[0,1,3,4,6,9]
LL2_map[3]=[0,1,2,4,5,7]

homs=[0,4,7,9]

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
fileout2=open(filenameout2,'w')
fileout3=open(filenameout3,'w')

###set up output headers
out1='##fileformat=VCFv4.1\n'
out2='##fileformat=VCFv4.1\n'
out3='##fileformat=Custom variant caller with indented reads, emit all haploid\n'
outA='##Caller_arguments=<MappingQuality_filter='+str(MQ)+',BaseQuality_filter='+str(BQ)
outA=outA+',GenotypeQuality_filter='+str(GQ)+',minRD='+str(min_RD)
outA=outA+',bam_in='+filenamein+',bedfile='+filenameinB
outA=outA+',indent='+str(indent+'>\n'
outB='##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n'
outB=outB+'##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n'
outB=outB+'##INFO=<ID=MQ,Number=1,Type=Float,Description="mean Mapping Quality (not RMS)">\n'
outB=outB+'##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">\n'
outB=outB+'##INFO=<ID=SGC,Number=1,Type=Integer,Description="low GQ heterozygote switched to best homozygote">\n'
                         
outC='##Time created='+(time.strftime("%H:%M:%S"))+' '+(time.strftime("%d/%m/%Y"))+'\n'
for i in range(len(ref.lengths)):
    outC=outC+'##contig=<ID='+ref.references[i]+',length='+str(ref.lengths[i])+'>\n'
outC=outC+'##reference=file:'+ref_file+'\n'

out_e='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+samp_name+'\n'
out_f='#CHROM\tPOS\tREF\tCALL\tBEST_PROB\tALL_PROB\tAD\n'

fileout.write(out1+outA+outB+outC+out_e)
fileout2.write(out2+outA+outB+outC+out_e)
fileout3.write(out3+outA+outC+out_f)

###Precomputed prios for QUAL estimation
theta_prior=[0,theta/1,theta/2]
theta_prior[0]=1-sum(theta_prior)


###Start running through regions in the bed
for gg in range(len(SNPlist)):
    print 'Calculating genotype likelihoods for locus '+str(gg+1)+' : '+SNPlist[gg]
    chromo=SNPdir[SNPlist[gg]][0]
    pos_start=int(SNPdir[SNPlist[gg]][1])
    pos_end=int(SNPdir[SNPlist[gg]][2])

    ###set up arrays to store results for a given bed entry (probably not to efficient for SNP data, more suited for long regions)
    seq_len=pos_end-pos_start
    GLs=np.zeros((seq_len,10),dtype='float32')  ##Genotype likelihoods
    PLs=np.zeros((seq_len,10),dtype='float32')  ##Phred-scale normalized likelihoods
    RDs=np.zeros((seq_len,4),dtype='int32')     ##Read depth for each of the four bases
    GQs=np.zeros((seq_len),dtype='float64')     ##Estimated genotype qualitues
    REFs=np.zeros((seq_len),dtype='int32')      ##Reference allele index
    GTs=np.zeros((seq_len),dtype='int32')       ##Genotype index [order is AA,AC,AG,AT,CC,CG,CT,GG,GT,TT]
    POS=np.zeros((seq_len),dtype='int32')       ##1-based position
    ALTs=np.zeros((seq_len),dtype='int32')      ##Alternate allele index. If trinucleotide, store as 10
    QUALs=np.zeros((seq_len),dtype='float32')   ##Bayesian quality score
    MQ0=np.zeros((seq_len),dtype='int32')  ##number of reads with mapping 0
    MQm=np.zeros((seq_len),dtype='float32')  #mean mapping score
    SEG=np.zeros((seq_len),dtype='int32')  #hom ref, het with ref, hom alt, het without ref (tri)
    SWGQ=np.zeros((seq_len),dtype='int32')  #low GQhet switched to hom ref = 1


    ###make a list of positions
    count=0
    for ggg in range(pos_start,pos_end):
        POS[count]=ggg+1
        count+=1

    ###prepopulate arrays with -9s
    GTs.fill(-9)
    REFs.fill(-9)
    ALTs.fill(-9)
    SEG.fill(-9)

    ###Go base by base through each region in the bam file
    for pileupcolumn in samfile.pileup(chromo,pos_start,pos_end,truncate=True,stepper='all'):
        nucl=ref.fetch(chromo,pileupcolumn.pos,pileupcolumn.pos+1)  ##grab reference sequence for region
        nucl=nucl.upper()
        var_list=[]  ##a list to store read bases
        map_list=[]  ##a list to store mapping qualities
        index=pileupcolumn.pos-pos_start  ##tells us what position we are in for the arrays given the nucleotide position
        REFs[index]=all_dic[nucl]  ##recode the reference allele
        ###go read by read for each position
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:  ##ensure not an indel or duplicate
                map_list.append(pileupread.alignment.mapping_quality)  ##record mapping quality of read
                if (pileupread.alignment.mapping_quality>=MQ) and (ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ) and (REFs[index]<>-99):  ###Add base calls meeting the MQ and BQ filters
                    if (pileupread.query_position>indent) and ((len(pileupread.alignment.query_sequence)-1-pileupread.query_position)>indent): ###check base is not at end of read
                        var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33])

        ###rescale qualities that are greater than 40 to a max of 40.
        for ggg in range(len(var_list)):
            if var_list[ggg][1]>40:
                var_list[ggg][1]=40

        ###record read depth for each basetype
        if len(var_list)>0:
            all_list=list(zip(*var_list)[0])
            RDs[index]=[all_list.count('A'),all_list.count('C'),all_list.count('G'),all_list.count('T')]
        
        ###work out MQ0 and mean mapping quality for the snp (not RMSE)
        MQ0[index]=map_list.count(0)
        try:
            MQm[index]=sum(map_list)/float(len(map_list))
        except:
            MQm[index]=0.0

        ###if minimum read depth is met, try calling the genotype 
        if len(var_list)>=min_RD:
#            GLs[index]=geno_caller_10GT_aDNA(var_list)  ##ancient DNA aware calculation of genotype likelihoods
            GLs[index]=geno_caller_10GT(var_list)  ##indent aware calculation of genotype likelihoods
            
            GTs[index]=np.argmax(GLs[index])    ##record best genotype
            PLs[index]=(GLs[index]-np.max(GLs[index]))*-10   ###calculte Phred-scale values
            GQs[index]=np.msort(PLs[index])[1]  ###record genotype quality

            ###if GQ is less than threshold and site is heterozygous, swith to best homozygous genotype
            if (alt_dic[GTs[index]][0]<>alt_dic[GTs[index]][1]) and (GQs[index]<GQ):
                GTs[index]=homs[np.argmax(GLs[index][homs])]
                SWGQ[index]=1

        ###Work out basesian QUALity score
##            LL_0=np.sum(10**GLs[index][LL0_map[REFs[index]]])
##            LL_1=np.sum(10**GLs[index][LL1_map[REFs[index]]])
##            LL_2=np.sum(10**GLs[index][LL2_map[REFs[index]]])
##            norconst1=sum([LL_0,LL_1,LL_2])  ###Depristo says this, but it really should be multiplied by the prior
##            norconst2=sum([LL_0*theta_prior[0],LL_1*theta_prior[1],LL_2*theta_prior[2]])  ###Depristo says this, but it really should be multiplied by the prior
##
##            Pr0=(theta_prior[0]*LL_0)/norconst2
##            Pr1=(theta_prior[1]*LL_1)/norconst2
##            Pr2=(theta_prior[2]*LL_2)/norconst2
            ###this new version add constant to log likelihoods to avoid underflow
            LL_0=np.sum(10**(GLs[index][LL0_map[REFs[index]]]-np.max(GLs[index])))  ##total likelihood for q=0|X
            LL_1=np.sum(10**(GLs[index][LL1_map[REFs[index]]]-np.max(GLs[index])))  ##total likelihood for q=1|X, three different genotypes summed
            LL_2=np.sum(10**(GLs[index][LL2_map[REFs[index]]]-np.max(GLs[index])))  ##total likelihood for q=2|X, six different genotypes summed

            ###calculate normalizing constant for bayes formula
##            norconst1=sum([LL_0,LL_1,LL_2])  ###Depristo says this, but it really should be multiplied by the prior
            norconst2=sum([LL_0*theta_prior[0],LL_1*theta_prior[1],LL_2*theta_prior[2]])  ###Depristo says this, but it really should be multiplied by the prior

            ###calculate posterio probabilities for q=0,1 and 2
            Pr0=(theta_prior[0]*LL_0)/norconst2
            Pr1=(theta_prior[1]*LL_1)/norconst2
            Pr2=(theta_prior[2]*LL_2)/norconst2

            ###work out QUAL. If Pr=0 is greatest, work out qual from probability for q not equalling 0, otherwise, work out from qual=0
            if LL0_map[REFs[index]][0]==GTs[index]:
                try:
                    QUALs[index]=-10*math.log10(1-Pr0)
                except:
                     QUALs[index]=100.0 #add hoc solution to deal with extremely small ref homozygote probability                  
                SEG[index]=0
            else:
                try:
                    QUALs[index]=-10*math.log10(Pr0)
                except:
                    QUALs[index]=1000.0  #add hoc solution to deal with extremely small ref homozygote probability

                ###Work out what the reference allele is, and if there is in fact a trinucleotide to
                if REFs[index]==alt_dic[GTs[index]][0]:  ##homo,alt heterozygote
                    ALTs[index]=alt_dic[GTs[index]][1]
                    SEG[index]=1
                elif REFs[index]==alt_dic[GTs[index]][1]: ##homo,alt heterozygote
                    ALTs[index]=alt_dic[GTs[index]][0]
                    SEG[index]=1
                elif alt_dic[GTs[index]][0]==alt_dic[GTs[index]][1]: #homozygote alternate
                    ALTs[index]=alt_dic[GTs[index]][0]
                    SEG[index]=2
                else:  #must be trinucleotide
                    ALTs[index]=10
                    SEG[index]=3

    ###start writing calls to  vcf file                   
    for ggg in range(len(GLs)):
        LQ=0  ##indicator for a low quality site
        out=chromo+'\t'+str(POS[ggg])+'\t.\t'+all_dic[REFs[ggg]]+'\t'
        if ALTs[ggg]<10:
            out=out+all_dic[ALTs[ggg]]  #write none or single alternate allele
        else:
            out=out+tri_dic[GTs[ggg]]  #write trinucleotide alleles

        if SEG[ggg]<>-9:   ###this sites was callable, even if QUAL was low
            out=out+'\t'+str(round(QUALs[ggg],2))+'\t'

            ###indicate whether site had low QUAL score or not
            if QUALs[ggg]>=QUALpass:
                out=out+'PASS\t'
                LQ=1
            else:
                out=out+'low_quality\t'

            ###Info tag information, alternate allele count, frequency, as well as mapping qualities
            if SEG[ggg]==0:
                out=out+'AC=0;AF=0.0;MQ='+str(round(MQm[ggg],2))+';MQ0='+str(MQ0[ggg])
            elif SEG[ggg]==1:
                out=out+'AC=1;AF=0.5;MQ='+str(round(MQm[ggg],2))+';MQ0='+str(MQ0[ggg])
            elif SEG[ggg]==2:
                out=out+'AC=2;AF=1.0;MQ='+str(round(MQm[ggg],2))+';MQ0='+str(MQ0[ggg])
            elif SEG[ggg]==3:
                out=out+'AC=1,1;AF=0.5,0.5;MQ='+str(round(MQm[ggg],2))+';MQ0='+str(MQ0[ggg])

            ###Make a note that heterozygote call was converted to homozygous call
            if SWGQ[ggg]==1:
                out=out+';SGC='+GL_dic[np.argmax(GLs[ggg])]

            if SEG[ggg]==0:  #site is homozygous reference, just note, GT,DP and GQ
                out=out+'\tGT:DP:GQ:PL\t0/0:'+str(RDs[ggg][REFs[ggg]])+':'
                if GQs[ggg]>99:
                    out=out+'99:'
                else:
                    out=out+str(int(round(GQs[ggg])))+':'
                ##output all 10 PL values. Note the order is not the same as GATK
                for gggg in range(len(PLs[ggg])):
                    out=out+str(int(round(PLs[ggg][gggg])))+','
                out=out[:-1]+'\n'

            else:  #Site is heterozygous so note more information
                out=out+'\tGT:AD:DP:GQ:PL\t'
                
                if SEG[ggg]<3:  #there is only one alternate allele
                    if SEG[ggg]==1:
                        out=out+'0/1:'
                    elif SEG[ggg]==2:
                        out=out+'1/1:'
                    out=out+str(RDs[ggg][REFs[ggg]])+','+str(RDs[ggg][ALTs[ggg]])+':'+str(np.sum(RDs[ggg]))+':'

                else:   #there are two alternate alleles
                    out=out+'1/2:'
                    out=out+str(RDs[ggg][REFs[ggg]])+','+str(RDs[ggg][alt_dic[GTs[ggg]][0]])+','+str(RDs[ggg][alt_dic[GTs[ggg]][1]])+':'+str(np.sum(RDs[ggg]))+':'

                ##cap reported GQ score to 99
                if GQs[ggg]>99:
                    out=out+'99:'
                else:
                    out=out+str(int(round(GQs[ggg])))+':'

                ##output all 10 PL values. Note the order is not the same as GATK
                for gggg in range(len(PLs[ggg])):
                    out=out+str(int(round(PLs[ggg][gggg])))+','
                out=out[:-1]+'\n'
                           

        else:  ###Site had not data, i.e. usable read depth 0
            out=out+'\t.\t.\t.\t.\t./.\n'
            
        fileout.write(out) ##write to emit all
        if (SEG[ggg]>=1) and (LQ==1):
            fileout2.write(out)  ##write to variable only vcf that passes QUAL threshold

        ###bayesian best haploid call
        out_h=chromo+'\t'+str(POS[ggg])+'\t'+all_dic[REFs[ggg]]+'\t'

        if SEG[ggg]>=0:
            ###Use GLs for homozygous genotypes to calculate probablity of A,C,G,T. Reference unaware, so prior probability of a particular base is equal to the other three bases
            norm_const3=np.sum(10**(GLs[ggg][homs]-np.max(GLs[ggg][homs]))*0.25)
            post=(10**(GLs[ggg][homs]-np.max(GLs[ggg][homs]))*0.25)/norm_const3

            ##this section checks if the highest posterior is unique or found for multiple bases
            max_val=np.max(post)
            matches=[]
            post_str=''
            for gggg in range(len(post)):
                if post[gggg]==max_val:
                    matches.append(gggg)
                post_str=post_str+str(round(post[gggg],3))+','

            ##if max val is 0.25, there is not information at this site, so set as N
            if max_val<=0.25:
                out_h=out_h+'N\t0.25\t0.25,0.25,0.25,0.25'
            else:
                if len(matches)==1: ##get the base with the highest LL
                    out_h=out_h+all_dic[matches[0]]+'\t'+str(round(post[matches[0]],3))+'\t'+post_str[:-1]
                else:
                    rand_allele=randint(0,len(matches)-1)  ##if two or more bases are tied for best likelihood, randomly choose one
                    out_h=out_h+all_dic[matches[rand_allele]]+'\t'+str(round(post[matches[rand_allele]],3))+'\t'+post_str[:-1]
            out_h=out_h+'\t'+str(RDs[ggg][0])+','+str(RDs[ggg][1])+','+str(RDs[ggg][2])+','+str(RDs[ggg][3])+'\n'
        else:
            out_h=out_h+'N\t0.25\t0.25,0.25,0.25,0.25\t0,0,0,0\n'
        fileout3.write(out_h)

        
fileout.close()
fileout2.close()
fileout3.close()
        


