#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:51:29 2019

@author: FinalZ
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:35:59 2019

@author: FinalZ
"""

from calculate_genotype_k_mer import Config_ComputeKmer,Genome_ComputeKmer
import getopt,sys,os

#python train_mod.py --infect human_infecting_virus --other Other_viruses --file mod_data --kmer 4
opts, args = getopt.getopt(sys.argv[1:], "hi:o:",["infect=","other=","file=","kmer="])

#get input parameter
for op, value in opts:
    if op == "--infect":
        infect = value
    elif op == "--other":
        other = value
    elif op == "--file":
        infile = value
    elif op == "--kmer":
        kmer = int(value)

#Determine 'infile' folder exists
isExists=os.path.exists('./'+infile+'')
if not isExists:
    os.makedirs(infile)

#for human infecting virus
in_file = open(''+infect+'','r')
in_files = in_file.readlines()
in_file.close()

#for other viruses
oth_file = open(''+other+'','r')
oth_files = oth_file.readlines()
oth_file.close()

#In order to reduce memory usage, we only calculate 1000 sequences at a time.
step = 1000
sub_infiles = [in_files[i:i+step] for i in range(0,len(in_files),step)]
sub_othfiles = [oth_files[i:i+step] for i in range(0,len(oth_files),step)]

piece = [500,1000,3000,5000,10000]
#get human infecting virus feature
isFirst = True
for sub in sub_infiles:
    #for make sure generate a new feature
    if isFirst:
        mode_type = 'w'
        isFirst = False
    else:
        mode_type = 'a'

    feature = Genome_ComputeKmer(sub,kmer,1)
    feature.to_csv('./'+infile+'/genome_'+str(kmer)+'_mer',mode = mode_type,sep='\t',header = None,index = False)
    #contig feature
    for p in piece:
        feature = Config_ComputeKmer(sub,kmer,p,1)
        feature.to_csv('./'+infile+'/'+str(p)+'bp_'+str(kmer)+'_mer',mode = mode_type,sep='\t',header = None,index = False)

#get other viruses feature
for sub in sub_othfiles:
    feature = Genome_ComputeKmer(sub,kmer,0)
    feature.to_csv('./'+infile+'/genome_'+str(kmer)+'_mer',mode = 'a',sep='\t',header = None,index = False)
    #contig feature
    for p in piece:
        feature = Config_ComputeKmer(sub,kmer,p,0)
        feature.to_csv('./'+infile+'/'+str(p)+'bp_'+str(kmer)+'_mer',mode = 'a',sep='\t',header = None,index = False)

#clear useless variable, Free memory
del feature,sub_infiles,sub_othfiles,in_files,oth_files
