#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:50:53 2019

@author: FinalZ
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 21:03:07 2017

@author: zhang
"""
def ComputeKmerVector(sequence,k,KmerCount):
    from  collections import OrderedDict
    n = len(sequence)
    for j in range(0,n-k+1):
        kmer = sequence[j:(j+k)]
        if kmer in KmerCount:
            KmerCount[kmer] += 1
        else:
            continue

    KmerVec = OrderedDict()
    for m in KmerCount:
        KmerVec[m] = float(KmerCount[m])/float(n-k+1)
    Vec = list(KmerVec.values())
    return Vec

def queryKmer(name,kmer):
    from  collections import OrderedDict
    import itertools
    import pandas as pd
    
    KmerCount = OrderedDict()
    s = itertools.product('ACGT',repeat = kmer)
    sl = list(s)
    for i in range(0,len(sl)):
        a0 = str(sl[i])
        a1 = a0.replace("(","")
        a2 = a1.replace(")","")
        a3 = a2.replace(",","")
        a4 = a3.replace("'","")
        a4 = a4.replace(" ", "")
        KmerCount[a4] = 0
    zero = [0 for _ in range(len(sl))]
    
    total_lst = list()
    for line in name:
        tt = line.strip('\n').split(' ')
        Seq = tt[1].upper()
        KmerCount = OrderedDict(zip(KmerCount.keys(),zero))
        temp = " ".join(map(str,ComputeKmerVector(Seq,kmer,KmerCount)))
        temp = tt[0]+' '+temp
        IVL = temp.split(' ')
        total_lst.append(IVL)
    ID_VEC_LABEL = pd.DataFrame(total_lst)
    return ID_VEC_LABEL


def Config_ComputeKmer(in_files,kmer,p,label):
    from  collections import OrderedDict
    import itertools,re
    import pandas as pd
    
    KmerCount = OrderedDict()
    s = itertools.product('ACGT',repeat = kmer)
    sl = list(s)
    for i in range(0,len(sl)):
        a0 = str(sl[i])
        a1 = a0.replace("(","")
        a2 = a1.replace(")","")
        a3 = a2.replace(",","")
        a4 = a3.replace("'","")
        a4 = a4.replace(" ", "")
        KmerCount[a4] = 0
    zero = [0 for _ in range(len(sl))]
    total_lst = list()
    
    for line in in_files:
        tt = line.strip('\n').split(' ')
        Seq = tt[1].upper()

        piece_Seq = re.findall(r'.{'+str(p)+'}',Seq)#sequence split into non-overlapping contigs of p
        if len(piece_Seq) == 0:
            continue
        count = 1
        for piece_seq in piece_Seq:
            KmerCount = OrderedDict(zip(KmerCount.keys(),zero))
            temp = " ".join(map(str,ComputeKmerVector(piece_seq,kmer,KmerCount)))
            temp = tt[0]+'_'+str(count)+' '+temp+' '+str(label)
            IVL = temp.split(' ')
            count += 1#mark contig
            total_lst.append(IVL)
    ID_VEC_LABEL = pd.DataFrame(total_lst)
    return ID_VEC_LABEL


def Genome_ComputeKmer(in_files,kmer,label):
    from  collections import OrderedDict
    import itertools
    import pandas as pd
    
    KmerCount = OrderedDict()
    s = itertools.product('ACGT',repeat = kmer)
    sl = list(s)
    for i in range(0,len(sl)):
        a0 = str(sl[i])
        a1 = a0.replace("(","")
        a2 = a1.replace(")","")
        a3 = a2.replace(",","")
        a4 = a3.replace("'","")
        a4 = a4.replace(" ", "")
        KmerCount[a4] = 0
    zero = [0 for _ in range(len(sl))]
    total_lst = list()
    
    for line in in_files:
        tt = line[0:-1].split(' ')
        Seq = tt[1].upper()
        KmerCount = OrderedDict(zip(KmerCount.keys(),zero))
        temp = " ".join(map(str,ComputeKmerVector(Seq,kmer,KmerCount)))
        temp = tt[0]+' '+temp+' '+str(label)
        IVL = temp.split(' ')
        total_lst.append(IVL)
    ID_VEC_LABEL = pd.DataFrame(total_lst)
    return ID_VEC_LABEL