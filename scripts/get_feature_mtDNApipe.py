# -*- coding: UTF-8 -*- 
#coding:gbk
#编写人员：guojiaming
#修改人员：sulp 2021.11.11   原版：/data/data2/guojiaming/data/data_c_c_c/code/get_feature6_0.1.py
#用    途：特征提取(识别低频真突变)  
#备    注：可根据文件名称识别突变阈值。即根据out_true_vaf的突变阈值，输出out_feature_vaf
#          本代码在计算特征的过程本不考虑没有发生突变的位点，如某位点全为A或全为T。
#          故输入这些位点是计算不出特征的。

import pysam 
import os
import sys
import uuid
import numpy as np
import pandas as pd
import regex as re  
from collections import defaultdict
from pyfaidx import Fasta
import pysam   
#import pysamstats  #  这个不重要
import scipy.stats
from scipy.stats import mannwhitneyu
from scipy.special import beta   
from subprocess import *
from multiprocessing import Pool
import subprocess
import math
from scipy.special import comb, perm


if len(sys.argv) >= 4:
    reference_fasta = sys.argv[3]
else:
    reference_fasta = os.environ.get("MTDNAPIPE_REFERENCE", "/mnt/hdd/ref_data/modified_genome/human_mtDNA.fasta")
reference = Fasta(reference_fasta)
reference_fai =reference_fasta+".fai" 
genome=reference_fai

file0=open(genome)    #打开参考基因
chr_sizes=dict()
for line in file0:   #对fasta数据里面的每一行。    .fai 文件。
    line=line.rstrip()
    fields=line.split('\t')
    chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])   #   chr_sizes 为一个字典。   #  dict.get(key, default=None)  key -- 字典中要查找的键。default -- 如果指定键的值不存在时，返回该默认值。
	#   fields[0]为： 序列的名字。 fields[1]：序列的长度。
file0.close()
#print(chr_sizes)

#Mutation_threshold

#def process_line(line,label,seq_file_format,bam_dir,yangben_name,reference_fasta):  #raw
def process_line(line,seq_file_format,bam_dir,yangben_name,reference_fasta):  #revised by sulp 2021.11.02
        line=line.rstrip()
        #print(line)
        fields=line.split('\t')
        sample=fields[4]   #样本 此处为 walsh.
        Mutation_frequency=fields[5]
        chr=fields[0]
        pos=int(fields[1])
        pos1=max(0,int(pos)-1)
        pos2=min(int(chr_sizes[chr]),int(pos)+1)    #  chr_sizes: {"chrM":17119}
        major_allele=fields[2]       # major_allele ：主等位碱基
        minor_allele=fields[3]       #  minor_allele  ：次等位碱基
        name=str(sample)+'~'+str(chr)+'~'+str(pos)+"~"+str(major_allele)+"~"+str(minor_allele)
        chrom=str(chr)     #   DNA类型  。
        start=int(pos)-1     #序列开始位置 
        end=int(pos)
        #print(name)
        #print(chrom)
        #print(start)
        #print(end)
        #if seq_file_format=="cram":   #raw
        if seq_file_format=="bam":    #revised by sulp 2021.11.03 将cram文件换成bam文件
                #input_cram=bam_dir+"/"+str(yangben_name)+".cram"      #raw
                #crai_file=bam_dir+"/"+str(yangben_name)+".crai"       #raw
                #crai_file2=bam_dir+"/"+str(yangben_name)+".cram.crai" #raw
                #if not os.path.exists(input_cram): #如果文件不存在则执行如下操作。    #raw
                #  print("no sample.cram under the cram_dir")                          #raw
                #if not os.path.exists(crai_file) and not os.path.exists(crai_file2):  #raw
	              #  print("no cram index files under the cram_dir")                     #raw
                input_bam=bam_dir+"/"+str(yangben_name)+".mt.no.softclip.bam"                     #revised by sulp 2021.11.03 将cram文件换成bam文件
                bai_file=bam_dir+"/"+str(yangben_name)+".bai"                       #revised by sulp 2021.11.03 将cram文件换成bam文件
                bai_file2=bam_dir+"/"+str(yangben_name)+".mt.no.softclip.bam.bai" #revised by sulp 2021.11.03 将cram文件换成bam文件
                if not os.path.exists(input_bam): #如果文件不存在则执行如下操作。    #revised by sulp 2021.11.03 将cram文件换成bam文件
                  print("no sample.bam under the bam_dir")                           #revised by sulp 2021.11.03 将cram文件换成bam文件
                if not os.path.exists(bai_file) and not os.path.exists(bai_file2):  #revised by sulp 2021.11.03 将cram文件换成bam文件
	                print("no bam index files under the bam_dir")                      #revised by sulp 2021.11.03 将cram文件换成bam文件
                #print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@") #delete by sulp 2021.11.2
                #read_length=int(subprocess.check_output("samtools view "+ input_cram+"|head -1000|awk '{print length($10)}'|sort|uniq -c|awk '{OFS=\"\\t\";print $1,$2}'|sort -k1,1nr|head -1|cut -f2 ",shell=True).decode('ascii')) #raw
                read_length=int(subprocess.check_output("samtools view \""+input_bam+"\"|head -1000|awk '{print length($10)}'|sort|uniq -c|awk '{OFS=\"\\t\";print $1,$2}'|sort -k1,1nr|head -1|cut -f2 ",shell=True).decode('ascii'))  #revised by sulp 2021.11.03 将cram文件换成bam文件
                #print("read_length:",read_length) #delete by sulp 2021.11.2
				
                #a=pysam.AlignmentFile(input_cram, "rc",reference_filename=reference_fasta) #raw
                a=pysam.AlignmentFile(input_bam, "rb",reference_filename=reference_fasta)  #revised by sulp 2021.11.03 将cram文件换成bam文件
                #subprocess.run("mkdir -p "+path+"/"+file+"/"+fil+"/tmp", shell=True)       #raw
                subprocess.run(["mkdir","-p",path+"/"+file+"/tmp"])               #revised by sulp 2021.11.03  删除positive/test 文件夹
                #tmp1_localcram_filename=path+"/"+file+"/"+fil+"/tmp/"+sample+"_"+chrom+"_"+str(pos)+"_"+str(uuid.uuid4())+".cram"    #raw
                tmp1_localbam_filename=path+"/"+file+"/tmp/"+sample+"_"+chrom+"_"+str(pos)+"_"+str(uuid.uuid4())+".bam"  #revised by sulp 2021.11.03  删除positive/test 文件夹
                #a_local=pysam.AlignmentFile(tmp1_localcram_filename,'wc',template=a,reference_filename=reference_fasta)   #raw
                a_local=pysam.AlignmentFile(tmp1_localbam_filename,'wb',template=a,reference_filename=reference_fasta)    #revised by sulp 2021.11.03 将cram文件换成bam文件
                for read in a.fetch(chrom,start,end):
                          a_local.write(read)
                a_local.close()
                #print("$$$$$$$$$$$$$$$$$$$") #delete by sulp 2021.11.2
                #pysam.index(tmp1_localcram_filename,tmp1_localcram_filename+".crai")   #raw
                pysam.index(tmp1_localbam_filename,tmp1_localbam_filename+".bai")      #revised by sulp 2021.11.03 
        major_plus[name]=0
        minor_plus[name]=0
        major_minus[name]=0
        minor_minus[name]=0
        major_read1[name]=0
        minor_read1[name]=0
        major_read2[name]=0
        minor_read2[name]=0
        major_ids[name]=list()
        minor_ids[name]=list()
        conflict_num[name]=0
        context1_count=dict()
        context2_count=dict()
        context1_count[name]=context1_count.get(name,0)
        context2_count[name]=context2_count.get(name,0)
        mismatches_major[name]=list()
        mismatches_minor[name]=list()
        minor2_count=dict()
        #minor2_count['A']=0
        max_num_2ndallele=0
        indels_count=dict()
      
        if len(major_allele)==len(minor_allele) and not(major_allele==minor_allele):
                 #print("%%%%%%%%%%%%%%%%%%%%%%%%%%") #delete by sulp 2021.11.2
                 if len(major_allele)==1:
                     state="SNV"
                 elif len(major_allele)>1:
                     state="MNP"
                 length=0
                 major_softclippedreads=0
                 minor_softclippedreads=0
                 major_num=0
                 minor_num=0
                 indels_count[name]=0
                 context_20bp=str(reference[chrom][max(1,int(pos)-11):min(int(pos)+10,int(chr_sizes[chrom]))])
                 GCcontent=(context_20bp.count('G')+context_20bp.count('C'))/len(context_20bp)
                 #context1[name]=reference[chrom][int(pos)-2:int(pos)+1]  #revised by sulp 2021.11.02
                 #context2[name]=(base[str(reference[chrom][int(pos)-2:int(pos)-1])]+base[str(reference[chrom][int(pos)-1:int(pos)])]+base[str(reference[chrom][int(pos):int(pos)+1])])[::-1]   #revised by sulp 2021.11.02
                 #print(GCcontent)  #delete by sulp 2021.11.2
                 #print(context1)   #delete by sulp 2021.11.2
                 #print(context2)   #delete by sulp 2021.11.2
                 if pos == 1 :    #add by sulp 2021.11.02
                     context1[name]='GGA'
                     context2[name]='G'+minor_allele+'A'
                 if pos != 1 :    #add by sulp 2021.11.02
                     context1[name]=str(reference[chrom][int(pos)-2])+major_allele+str(reference[chrom][int(pos)])
                     context2[name]=str(reference[chrom][int(pos)-2])+minor_allele+str(reference[chrom][int(pos)]) 

                 #if seq_file_format=="cram":  #raw
                 if seq_file_format=="bam":    #revised by sulp 2021.11.03
                     #a=pysam.AlignmentFile(tmp1_localcram_filename, "rc",reference_filename=reference_fasta)  #raw
                     a=pysam.AlignmentFile(tmp1_localbam_filename, "rb",reference_filename=reference_fasta)   #revised by sulp 2021.11.03

                 for pileupcolumn in a.pileup(chrom, start, end, max_depth=8000):
                   for pileupread in pileupcolumn.pileups:
                        if pileupread.indel !=0:
                            indels_count[name]=indels_count.get(name,0)+1
                            #print("###########################pileupread.indel !=0##############")
                            #print(name,pileupread.indel,indels_count[name],pileupcolumn.pos)
                            #print("#############################")
                            continue
                        try:
                          querybase=pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+len(major_allele)]
                          #print("######################################################################################")
                          #print("querybase:",querybase)
                          #print("######################################################################################")
                          if pileupcolumn.pos==pos-1 and (not pileupread.alignment.flag & 256) and (not pileupread.alignment.flag & 1024):
                             #print("##################")
                              if querybase==major_allele: 
                                     major_num=major_num+1   #1132
                                     if (pileupread.alignment.get_cigar_stats()[0][4])>=10:
                                              major_softclippedreads=major_softclippedreads+1  #
                                     major_ids[name].append(pileupread.alignment.query_name)
                                     querypos_major[name].append(pileupread.query_position)
                                     mapq_major[name].append(pileupread.alignment.mapping_quality)
                                     baseq_major[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
                                     leftpos_major[name].append(pileupread.alignment.reference_start)
                                     mismatches_major[name].append(int(pileupread.alignment.get_tag('NM'))/read_length)
                                     if not pileupread.alignment.is_reverse:
                                           major_plus[name]=major_plus.get(name,0)+1
                                     elif pileupread.alignment.is_reverse:
                                           major_minus[name]=major_minus.get(name,0)+1
 
                                     if pileupread.alignment.flag & 64:
                                           major_read1[name]=major_read1.get(name,0)+1
                                     elif pileupread.alignment.flag &128:
                                           major_read2[name]=major_read2.get(name,0)+1
                                
                                     if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
                                          seqpos_major[name].append(pileupread.query_position)

                                          if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
                                               baseq_major_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+1])
                                          elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
                                               baseq_major_near1b[name].append("end")
                                     elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
                                          ##elif pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
                                          seqpos_major[name].append(len(pileupread.alignment.query_sequence)-pileupread.query_position)
                                          if pileupread.query_position >=1 :
                                             baseq_major_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-1])
                                                #print baseq_major_near1b[name]
                                          if pileupread.query_position ==0 :
                                             baseq_major_near1b[name].append("end")

                              elif querybase==minor_allele:
                                    minor_num=minor_num+1
                                    if int(pileupread.alignment.get_cigar_stats()[0][4])>10:
                                        minor_softclippedreads=minor_softclippedreads+1
                                    minor_ids[name].append(pileupread.alignment.query_name)
                                    querypos_minor[name].append(pileupread.query_position)
                                    mapq_minor[name].append(pileupread.alignment.mapping_quality)
                                    baseq_minor[name].append(pileupread.alignment.query_qualities[pileupread.query_position])
                                    leftpos_minor[name].append(pileupread.alignment.reference_start)
									                  #mismatches_minor[name].append(filterPick(pileupread.alignment.tags,'NM'))
									                  #mismatches_minor[name].append(pileupread.alignment.get_tag('NM'))
                                    mismatches_minor[name].append(int(pileupread.alignment.get_tag('NM'))/read_length)
                                    if not pileupread.alignment.is_reverse:
                                         minor_plus[name]=minor_plus.get(name,0)+1
                                    elif pileupread.alignment.is_reverse:
                                         minor_minus[name]=minor_minus.get(name,0)+1
									
                                    if pileupread.alignment.flag & 64:
                                         minor_read1[name]=minor_read1.get(name,0)+1
                                    elif pileupread.alignment.flag &128:
                                         minor_read2[name]=minor_read2.get(name,0)+1
									
                                    if pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
									                       #if pileupread.alignment.reference_start-pileupread.alignment.next_reference_start<0: 
                                         seqpos_minor[name].append(pileupread.query_position)
                                         if pileupread.query_position < len(pileupread.alignment.query_sequence)-1:
                                             baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position+1])
                                         elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
                                             baseq_minor_near1b[name].append("end")
				                                     #      elif pileupread.query_position == len(pileupread.alignment.query_sequence)-1:
				                                     #       baseq_minor_near1b[name].append("end")
                                    elif pileupread.alignment.is_proper_pair and pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
									                 #elif pileupread.alignment.reference_start-pileupread.alignment.next_reference_start>0: 
                                        seqpos_minor[name].append(len(pileupread.alignment.query_sequence)-pileupread.query_position)
                                        if pileupread.query_position >=1 :
                                                baseq_minor_near1b[name].append(pileupread.alignment.query_qualities[pileupread.query_position-1])
                                        elif pileupread.query_position ==0:
                                                baseq_minor_near1b[name].append("end")
                 
                                    if pileupread.alignment.is_proper_pair:
									                     ##http://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
                                        if pileupread.alignment.flag & 64 and (not pileupread.alignment.is_reverse):
                                             context1_count[name]=context1_count.get(name,int(0))+int(1)
                                        elif pileupread.alignment.flag & 128 and (pileupread.alignment.is_reverse):
                                             context2_count[name]=context2_count.get(name,int(0))+int(1)
                                        elif pileupread.alignment.flag & 64 and (pileupread.alignment.is_reverse):
                                             context2_count[name]=context2_count.get(name,int(0))+int(1)
                                        elif pileupread.alignment.flag & 128 and (not pileupread.alignment.is_reverse):
                                             context1_count[name]=context1_count.get(name,int(0))+int(1)

                         
                        except:
                              continue 
                 #print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7")  #delete by sulp 2021.11.02
				 
                 

                 conflict_reads=set(major_ids[name]) & set(minor_ids[name])
                 conflict_num[name]=len(conflict_reads)
                 #print(conflict_reads,conflict_num)
                 
                 if major_num != 0:   #add by sulp 2021.11.02
                     ref_softclip=str(major_softclippedreads/major_num)
                 else:
                     ref_softclip="NA" 
                 if minor_num !=0 :   #add by sulp 2021.11.02
                     alt_softclip=str(minor_softclippedreads/minor_num)
                 else:
                     alt_softclip="NA"   

                 #if seq_file_format=="cram":  #raw
                 if seq_file_format=="bam":   #revised by sulp 2021.11.03
                    #subprocess.run("rm "+tmp1_localcram_filename, shell=True)              #raw
                    subprocess.run("rm "+tmp1_localbam_filename, shell=True)               #revised by sulp 2021.11.03
                    #subprocess.run("rm "+tmp1_localcram_filename+".crai", shell=True)      #raw
                    subprocess.run("rm "+tmp1_localbam_filename+".bai", shell=True)        #revised by sulp 2021.11.03
		                #subprocess.run("rm "+tmp2_localcram_filename, shell=True)
		                #subprocess.run("rm "+tmp2_localcram_filename+".crai", shell=True)

                 #print("name:",name) #delete by sulp 2021.11.02

        #revised by sulp 2021.11.02
        #return label,name,Mutation_frequency,','.join(str(x) for x in querypos_major[name])+",",  ','.join(str(x) for x in querypos_minor[name])+",",  ','.join(str(x) for x in leftpos_major[name])+",",  ','.join(str(x) for x in leftpos_minor[name])+",",  ','.join(str(x) for x in seqpos_major[name])+",",  ','.join(str(x) for x in seqpos_minor[name])+",",  ','.join(str(x) for x in mapq_major[name])+",",  ','.join(str(x) for x in mapq_minor[name])+",", ','.join(str(x) for x in baseq_major[name])+",",  ','.join(str(x) for x in baseq_minor[name])+",",  ','.join(str(x) for x in baseq_major_near1b[name])+",", ','.join(str(x) for x in baseq_minor_near1b[name])+",", major_plus[name],major_minus[name],minor_plus[name],minor_minus[name], str(context1[name]), str(context2[name]), context1_count[name],context2_count[name],','.join(str(x) for x in mismatches_major[name])+",",','.join(str(x) for x in mismatches_minor[name])+",", major_read1[name],major_read2[name],minor_read1[name],minor_read2[name],conflict_num[name], state, str(length),str(GCcontent), str(major_softclippedreads/major_num), str(minor_softclippedreads/minor_num), str(indels_count[name]/(major_num+minor_num+indels_count[name])), max_num_2ndallele
        return name,Mutation_frequency,','.join(str(x) for x in querypos_major[name])+",",  ','.join(str(x) for x in querypos_minor[name])+",",  ','.join(str(x) for x in leftpos_major[name])+",",  ','.join(str(x) for x in leftpos_minor[name])+",",  ','.join(str(x) for x in seqpos_major[name])+",",  ','.join(str(x) for x in seqpos_minor[name])+",",  ','.join(str(x) for x in mapq_major[name])+",",  ','.join(str(x) for x in mapq_minor[name])+",", ','.join(str(x) for x in baseq_major[name])+",",  ','.join(str(x) for x in baseq_minor[name])+",",  ','.join(str(x) for x in baseq_major_near1b[name])+",", ','.join(str(x) for x in baseq_minor_near1b[name])+",", major_plus[name],major_minus[name],minor_plus[name],minor_minus[name], str(context1[name]), str(context2[name]), context1_count[name],context2_count[name],','.join(str(x) for x in mismatches_major[name])+",",','.join(str(x) for x in mismatches_minor[name])+",", major_read1[name],major_read2[name],minor_read1[name],minor_read2[name],conflict_num[name],str(GCcontent), ref_softclip, alt_softclip


path="/data/data2/guojiaming/data/data_c_c_c/FFPE_moqinqin"   #raw
path = sys.argv[1]   #revised by sulp 2021.11.02  传递参数

#data_all=os.listdir(path)                                    #raw
#data_all=[ i for i in os.listdir(path) if os.path.isdir(i) ] #revised by sulp 2021.11.02 只列出文件夹
data_all=[ sys.argv[2] ]                                     #revised by sulp 2021.11.11  传递参数  一个文件夹的样本名  为了便于shell的for循环并行运算
for file in data_all:     #file :cc1806
    data_little=os.listdir(path+"/"+file)
      #for fil in data_little:     # fil : postive    #delete by sulp 2021.11.03  删除positive/test的文件夹
      # if fil=="positive":    #delete by sulp 2021.11.02
      #    label=str(1)        #delete by sulp 2021.11.02
      # else:                  #delete by sulp 2021.11.02
      #    label=str(0)        #delete by sulp 2021.11.02
      #data_more=os.listdir(path+"/"+file+"/"+fil) #delete by sulp 2021.11.03
      #for fi in data_more:                        #raw
    for fi in data_little:                       #revised by sulp 2021.11.03
          #fff=fi.split(".")                       #delete by sulp 2021.11.02
          #if fff[-1]=="cram":                     #delete by sulp 2021.11.02
          #   yangben_name=fff[0]                  #delete by sulp 2021.11.02
          filename=fi.split("_")                                 #add by sulp 2021.11.02
          if filename[0]=='output' and filename[1]=='true' :    #add by sulp 2021.11.02
            vaf = filename[2]      #提取样本的突变阈值           #add by sulp 2021.11.02
          #if fi=="output_true_0.1":   # 需要修改                #delete by sulp 2021.11.02
            #print(file)
            #print(fi)
            #print(label)
            #input_pos=path+"/"+file+"/"+fil+"/"+fi   #raw
            input_pos=path+"/"+file+"/"+fi           #revised by sulp 2021.11.03  # 删除positive/test的文件夹
            #bam_dir=path+"/"+file+"/"+fil            #raw
            bam_dir=path+"/"+file                    #revised by sulp 2021.11.03  # 删除positive/test的文件夹
            #output=path+"/"+file+"/"+fil+"/output_feature_0.1"      # 需要修改   #raw
            output=path+"/"+file+"/output_feature_"+vaf                          #revised by sulp 2021.11.03   
            # NOTE: reference_fasta is resolved at the top of this script (argv[3] or MTDNAPIPE_REFERENCE)
            #seq_file_format="cram" #raw
            seq_file_format="bam"  #revised by sulp 2021.11.03
            n_jobs=1
            base=dict() #dict() 用于创建字典
            base['A']='T'
            base['T']='A'
            base['G']='C'
            base['C']='G'
			
            fo=open(output+".tmp","w")
            #revised by sulp 2021.11.02
            #header='label id Mutation_frequency querypos_major querypos_minor leftpos_major leftpos_minor seqpos_major seqpos_minor mapq_major mapq_minor baseq_major baseq_minor baseq_major_near1b baseq_minor_near1b major_plus major_minus minor_plus minor_minus context1 context2 context1_count context2_count mismatches_major mismatches_minor major_read1 major_read2 minor_read1 minor_read2 conflict_num type length GCcontent ref_softclip alt_softclip indel_proportion_SNPonly alt2_proportion_SNPonly'.split()
            header='id Mutation_frequency querypos_major querypos_minor leftpos_major leftpos_minor seqpos_major seqpos_minor mapq_major mapq_minor baseq_major baseq_minor baseq_major_near1b baseq_minor_near1b major_plus major_minus minor_plus minor_minus context1 context2 context1_count context2_count mismatches_major mismatches_minor major_read1 major_read2 minor_read1 minor_read2 conflict_num GCcontent ref_softclip alt_softclip'.split()
            #print (' '.join(header),file=fo)  #raw
            print ('\t'.join(header),file=fo) #revised by sulp 2021.11.11 修改分隔符
   
   
            querypos_major=defaultdict(list)
            mapq_major=defaultdict(list)
            baseq_major=defaultdict(list)  #都是字典
            leftpos_major=defaultdict(list)
            mismatches_major=defaultdict(list)
            major_plus=dict() #普通的字典
            major_minus=dict()
            major_read1=dict()
            major_read2=dict()
            major_ids=defaultdict(list)
            minor_ids=defaultdict(list)
            conflict_num=dict()
            baseq_major_near1b=defaultdict(list)
            seqpos_major=defaultdict(list)
            context1=dict()
            context2=dict()
            
            context2_count=dict()
            querypos_minor=defaultdict(list)
            dp_far=defaultdict(list)
            dp_near=defaultdict(list)
            mapq_minor=defaultdict(list)
            baseq_minor=defaultdict(list)
            leftpos_minor=defaultdict(list)
            mismatches_minor=defaultdict(list)
            minor_plus=dict()
            minor_minus=dict()
            minor_read1=dict()
            minor_read2=dict()
            baseq_minor_near1b=defaultdict(list)
            seqpos_minor=defaultdict(list)
       


            with open(input_pos,"r") as f:
                 for line in f.readlines():     
                   try:
                     if line.strip() == "":
                         continue

                     #result=process_line(line,label,seq_file_format,bam_dir,yangben_name,reference_fasta)  #raw
                     result=process_line(line,seq_file_format,bam_dir,file,reference_fasta)  #revised by sulp 2021.11.02

                     #print("result:",process_line(line,label,seq_file_format,bam_dir,yangben_name,reference_fasta))
                     #print(type(result))
                     for atuple in result:
                        #print(atuple)
                        #print(type(atuple))

                        fo.write(str(atuple))
                            #fo.write(''.join(str(x) for x in atuple))
							
                        #fo.write(" ") #raw
                        fo.write("\t") #revised by sulp 2021.11.11 修改分隔符
                            #print (''.join(str(x) for x in atuple),file=fo)
                     
                     fo.write("\n")
                   except:
                         continue
           

            fo.close()


            try:
                  #df=pd.read_csv(output+".tmp",sep=" ",index_col=False)  #output=sys.argv[2]  ：写入的方式     #raw 
                  df=pd.read_csv(output+".tmp",sep="\t",index_col=False)  #output=sys.argv[2]  ：写入的方式    #revised by sulp 2021.11.11 修改分隔符
                  #id
                  #print(df.label)  #delete by sulp 2021.11.02
                  df = df[df.querypos_minor != ',']
                  df = df[df.seqpos_minor != ',']
                  df = df[df.seqpos_major != ',']
                  df = df[df.baseq_minor_near1b != ',']
                  df = df[df.baseq_minor_near1b != 'end,']
                  df = df[df.leftpos_minor != ',']
                  df = df[df.baseq_major_near1b != ',']
                  df = df[df.baseq_major_near1b != 'end,']
                  
                  # print("*************************************") #delete by sulp 2021.11.02
                   
                  def my_mosaic_likelihood(a,b,c,d,e,f):
                  	depth=sum([int(a),int(b),int(c),int(d)])
                  	alt=sum([int(c),int(d)])
                  	r=0
                  	baseq_major=[float(i) for i in e.split(',')[:-1]]
                  	baseq_minor=[float(i) for i in f.split(',')[:-1]]
                  	r=sum([0.1**(float(i)/10) for i in baseq_major])
                  	r=r+sum([1-0.1**(float(i)/10) for i in baseq_minor])
                  	l=beta(r+1, depth-r+1)
                  	if l >0:
                  		l=math.log10(l)+math.log10(comb(depth,alt,exact=True))
                  		l=10**l
                  ##	return(float(beta(r+1, depth-r+1)))
                  	return(l)
                  
                  def my_het_likelihood(a,b,c,d):
                  	depth=sum([int(a),int(b),int(c),int(d)])
                  	alt=sum([int(c),int(d)])
                  	l=math.log10(comb(depth,alt,exact=True))+math.log10(0.5)*depth
                  	l=10**l
                  	return(l)
                  #	math.log10(comb(2000,1000,exact=True))+math.log10(0.5)*2000
                  #	return(0.5**depth)
                  
                  def my_refhom_likelihod(a,b):
                  	baseq_major=[float(i) for i in a.split(',')[:-1]]
                  	baseq_minor=[float(i) for i in b.split(',')[:-1]]
                  	depth=len(baseq_major)+len(baseq_minor)
                  	alt=len(baseq_minor)
                  	q=math.log10(1)
                  	q=sum(math.log10(1-0.1**(i/10)) for i in baseq_major)
                  	q=q+sum(math.log10(0.1**(i/10)) for i in baseq_minor)
                  	l=math.log10(comb(depth,alt,exact=True))+q
                  	#return(10**q)	
                  	return(10**l)	
                  
                  def my_althom_likelihod(a,b):
                  	baseq_major=[float(i) for i in a.split(',')[:-1]]
                  	baseq_minor=[float(i) for i in b.split(',')[:-1]]
                  	depth=len(baseq_major)+len(baseq_minor)
                  	alt=len(baseq_minor)
                  	q=math.log10(1)
                  	q=sum(math.log10(1-0.1**(i/10)) for i in baseq_minor)
                  	q=q+sum(math.log10(0.1**(i/10)) for i in baseq_major)
                  	l=math.log10(comb(depth,alt,exact=True))+q
                  	return(10**l)	
                  
                  def my_wilcox_pvalue(a, b):
                  	x1=[float(i) for i in a.split(',')[:-1]]
                  	x2=[float(i) for i in b.split(',')[:-1]]
                  	if x1!=x2:
                  		return (scipy.stats.ranksums(x1,x2)[1])
                  	elif x1==x2:
                  		return(float(1))
                  def my_wilcox_statistics(a, b):
                  	x1=[i for i in a.split(',')[:-1] if i!="end"]
                  	x2=[i for i in b.split(',')[:-1] if i!="end"]
                  	x1=[float(i) for i in x1]
                  	x2=[float(i) for i in x2]
                  	if x1!=x2:
                  		return (scipy.stats.ranksums(x1,x2)[0])
                  	elif x1==x2:
                  		return(float(0))
                  def my_wilcox_paired_pvalue(a, b):
                  	x2_index = [x for x in range(len(b.split(','))) if b.split(',')[x]!="end"][:-1]
                  	x2= [i for i in b.split(',') if i !="end"][:-1]
                  	x1=[float(a.split(',')[i]) for i in x2_index]
                  	x2=[float(i) for i in x2]
                  	if x1!=x2:
                  		return (scipy.stats.wilcoxon(x1,x2)[1])
                  	elif x1==x2:
                  		return(float(1))
                  def my_wilcox_paired_statistics(a, b):
                  	x2_index = [x for x in range(len(b.split(','))) if b.split(',')[x]!="end"][:-1]
                  	x2= [i for i in b.split(',') if i !="end"][:-1]
                  	x1=[float(a.split(',')[i]) for i in x2_index]
                  	x2=[float(i) for i in x2]
                  	if x1!=x2:
                  		return (scipy.stats.wilcoxon(x1,x2)[0])
                  	if x1==x2:
                  		return(float(0))
                  #	return (scipy.stats.mannwhitneyu(x1,x2)[0])
                  def my_ttest_statistics(a, b):
                  	if a!=b:
                  		x1=[float(i) for i in a.split(',')[:-1]]
                  		x2=[float(i) for i in b.split(',')[:-1]]
                  		return (scipy.stats.ttest_ind(x1,x2, equal_var = False)[0])
                  	elif a==b:
                  		return(float(0))
                  def my_fisher_pvalue(a,b,c,d):
                  	return (scipy.stats.fisher_exact([[int(a), int(b)], [int(c), int(d)]])[1])
                  def my_fisher_statistics(a,b,c,d):
                  	return (scipy.stats.fisher_exact([[int(a), int(b)], [int(c), int(d)]])[0])
                  def my_context_selection(a,b,c,d):
                  	if int(a)>=int(b):
                  		return(c)
                  	else:
                  		return(d)
                  def my_mean(a):
                  	x=[float(i) for i in a.split(',')[:-1]]	
                  	#return(sum(x)/len(x)/float(read_length))
                  	return(sum(x)/len(x))
                  def my_AF(a,b,c,d):
                  	depth=sum([int(a),int(b),int(c),int(d)])
                  	alt=sum([int(c),int(d)])
                  	return(float(alt)/float(depth))
                  def my_depth(a,b,c,d):
                  	depth=sum([int(a),int(b),int(c),int(d)])
                  	return(depth)
                  def my_mean_difference(a,b):
                  	x1=[float(i) for i in a.split(',')[:-1]]	
                  	x2=[float(i) for i in b.split(',')[:-1]]	
                  	return(sum(x1)/len(x1)-sum(x2)/len(x2))
                  def my_difference(a,b):
                  	return(float(a)-float(b))
                  	
                  	
                  	###  mappability :mappability[name]      np.mean(dp_near[name])  np.mean(dp_far[name])
                  	#dp_near dp_far mappability
                  	
                  df['querypos_p']=df.apply(lambda row: my_wilcox_pvalue(row['querypos_major'], row['querypos_minor']), axis=1)
                  df['leftpos_p']=df.apply(lambda row: my_wilcox_pvalue(row['leftpos_major'], row['leftpos_minor']), axis=1)
                  df['seqpos_p']=df.apply(lambda row: my_wilcox_pvalue(row['seqpos_major'], row['seqpos_minor']), axis=1)
                  df['mapq_p']=df.apply(lambda row: my_wilcox_pvalue(row['mapq_major'], row['mapq_minor']), axis=1)
                  df['baseq_p']=df.apply(lambda row: my_wilcox_pvalue(row['baseq_major'], row['baseq_minor']), axis=1)
                  df['baseq_t']=df.apply(lambda row: my_wilcox_statistics(row['baseq_major'], row['baseq_minor']), axis=1)
                  #df['baseq_t']=df.apply(lambda row: my_ttest_statistics(row['baseq_major'], row['baseq_minor']), axis=1)
                  #df['ref_baseq1b_p']=df.apply(lambda row: my_wilcox_pvalue(row['baseq_major'], row['baseq_major_near1b']), axis=1)
                  df['ref_baseq1b_t']=df.apply(lambda row: my_wilcox_statistics(row['baseq_major'], row['baseq_major_near1b']), axis=1)
                  #df['ref_baseq1b_t']=df.apply(lambda row: my_ttest_statistics(row['baseq_major'], row['baseq_major_near1b']), axis=1)
                  #df['alt_baseq1b_p']=df.apply(lambda row: my_wilcox_pvalue(row['baseq_minor'], row['baseq_minor_near1b']), axis=1)	
                  df['ref_baseq1b_p']=df.apply(lambda row: my_wilcox_paired_pvalue(row['baseq_major'], row['baseq_major_near1b']), axis=1)
                  #df['ref_baseq1b_t']=df.apply(lambda row: my_wilcox_paired_statistics(row['baseq_major'], row['baseq_major_near1b']), axis=1)
                  df['alt_baseq1b_p']=df.apply(lambda row: my_wilcox_paired_pvalue(row['baseq_minor'], row['baseq_minor_near1b']), axis=1)	
                  #df['alt_baseq1b_t']=df.apply(lambda row: my_wilcox_paired_statistics(row['baseq_minor'], row['baseq_minor_near1b']), axis=1)	
                  df['alt_baseq1b_t']=df.apply(lambda row: my_wilcox_statistics(row['baseq_minor'], row['baseq_minor_near1b']), axis=1)	
                  #df['alt_baseq1b_t']=df.apply(lambda row: my_ttest_statistics(row['baseq_minor'], row['baseq_minor_near1b']), axis=1)	
                  df['sb_p']=df.apply(lambda row: my_fisher_pvalue(row['major_plus'], row['major_minus'], row['minor_plus'], row['minor_minus']), axis=1)	
                  df['context']=df.apply(lambda row: my_context_selection(row['context1_count'], row['context2_count'], row['context1'], row['context2']), axis=1)	
                  df['major_mismatches_mean']=df.apply(lambda row: my_mean(row['mismatches_major']), axis=1)	
                  df['minor_mismatches_mean']=df.apply(lambda row: my_mean(row['mismatches_minor']), axis=1)	
                  df['mismatches_p']=df.apply(lambda row: my_wilcox_pvalue(row['mismatches_major'], row['mismatches_minor']), axis=1)	
                  df['AF']=df.apply(lambda row: my_AF(row['major_plus'], row['major_minus'], row['minor_plus'], row['minor_minus']), axis=1)	
                  df['dp']=df.apply(lambda row: my_depth(row['major_plus'], row['major_minus'], row['minor_plus'], row['minor_minus']), axis=1)	
                  df['mosaic_likelihood']=df.apply(lambda row: my_mosaic_likelihood(row['major_plus'], row['major_minus'], row['minor_plus'],row['minor_minus'],row['baseq_major'],row['baseq_minor']), axis=1)
                  df['het_likelihood']=df.apply(lambda row: my_het_likelihood(row['major_plus'], row['major_minus'], row['minor_plus'],row['minor_minus']), axis=1)
                  df['refhom_likelihood']=df.apply(lambda row: my_refhom_likelihod(row['baseq_major'],row['baseq_minor']), axis=1)
                  df['althom_likelihood']=df.apply(lambda row: my_althom_likelihod(row['baseq_major'],row['baseq_minor']), axis=1)
                  df['normalize']=df['mosaic_likelihood']+df['het_likelihood']+df['refhom_likelihood']+df['althom_likelihood']
                  df['mosaic_likelihood']=df['mosaic_likelihood']/df['normalize']
                  df['het_likelihood']=df['het_likelihood']/df['normalize']
                  df['refhom_likelihood']=df['refhom_likelihood']/df['normalize']
                  df['althom_likelihood']=df['althom_likelihood']/df['normalize']

                  
                  df['mapq_difference']=df.apply(lambda row: my_mean_difference(row['mapq_major'], row['mapq_minor']), axis=1)
                  df['sb_read12_p']=df.apply(lambda row: my_fisher_pvalue(row['major_read1'], row['major_read2'], row['minor_read1'], row['minor_read2']), axis=1)
                  
                  ## 没有计算df['dp_diff']=df.apply(lambda row: my_difference(row['dp_near'], row['dp_far']), axis=1)
                  #df['dp_p'].replace('', np.nan, inplace=True)
                  #df['alt_baseq1b_t'].replace('', np.nan, inplace=True)
                  
                  df['alt_baseq1b_t'].fillna(0,inplace=True)
                  #df['dp_p'].fillna(1,inplace=True)
                  
                  #  dp_diff ：没有计算
                  try:
                      #revised by sulp 2021.11.02
                      #df_new=df[['label','id','Mutation_frequency','conflict_num','type','length','GCcontent','ref_softclip','alt_softclip','querypos_p','leftpos_p','seqpos_p','mapq_p','baseq_p','baseq_t','ref_baseq1b_p','ref_baseq1b_t', 'alt_baseq1b_p','alt_baseq1b_t','sb_p','context','major_mismatches_mean','minor_mismatches_mean','mismatches_p','AF','dp','mosaic_likelihood','het_likelihood','refhom_likelihood','althom_likelihood', 'mapq_difference', 'sb_read12_p', 'indel_proportion_SNPonly', 'alt2_proportion_SNPonly']]
                      df_new=df[['id','Mutation_frequency','conflict_num','GCcontent','ref_softclip','alt_softclip','querypos_p','leftpos_p','seqpos_p','mapq_p','baseq_p','baseq_t','ref_baseq1b_p','ref_baseq1b_t', 'alt_baseq1b_p','alt_baseq1b_t','sb_p','context','major_mismatches_mean','minor_mismatches_mean','mismatches_p','AF','dp','mosaic_likelihood','het_likelihood','refhom_likelihood','althom_likelihood', 'mapq_difference', 'sb_read12_p']]
                      
                  #df_new=df[['type','id''conflict_num','mappability','type','length','GCcontent','ref_softclip','alt_softclip','querypos_p','leftpos_p','seqpos_p','mapq_p','baseq_p','baseq_t','ref_baseq1b_p','ref_baseq1b_t', 'alt_baseq1b_p','alt_baseq1b_t','sb_p','context','major_mismatches_mean','minor_mismatches_mean','mismatches_p','AF','dp','mosaic_likelihood','het_likelihood','refhom_likelihood','althom_likelihood', 'mapq_difference', 'sb_read12_p', 'indel_proportion_SNPonly', 'alt2_proportion_SNPonly']]
                  #id querypos_major querypos_minor leftpos_major leftpos_minor seqpos_major seqpos_minor mapq_major mapq_minor baseq_major baseq_minor baseq_major_near1b baseq_minor_near1b major_plus major_minus minor_plus minor_minus context1 context2 context1_count context2_count mismatches_major mismatches_minor major_read1 major_read2 minor_read1 minor_read2 dp_near dp_far dp_p conflict_num mappability type length GCcontent ref_softclip alt_softclip
                      fo2=open(output,"w")
                      df_new.to_csv(fo2, index=False,sep="\t")
                      fo2.close()
                  except:
                     continue
            except:
                continue
		
          else:
                #print(fi)
                #print("The corresponding output file was not found")
                continue
subprocess.run(["rm","-rf",path+"/"+file+"/tmp"])            #add by sulp 2021.11.03 #删除样本文件夹中的tmp文件夹
          
                        
                        
               
              
            
            
            



