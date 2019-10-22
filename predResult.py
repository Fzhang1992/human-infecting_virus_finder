#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 20:13:34 2019

@author: FinalZ
"""
from calculate_genotype_k_mer import queryKmer
import pandas as pd
from imblearn.ensemble import BalancedBaggingClassifier
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
import getopt,sys
import time

start_time = time.time()

#python predResult.py --query test --output predict_result --file mod_data --kmer 4 --kjob 1 --bjob 1
opts, args = getopt.getopt(sys.argv[1:], "hi:o:",["query=","output=","file=","kmer=","kjob=","bjob="])
for op, value in opts:
    if op == "--query":
        query = value
    elif op == "--output":
        output = value
    elif op == "--file":
        infile = value
    elif op == "--kmer":
        kmer = int(value)
    elif op == "--kjob":
        k_job = int(value)
    elif op == "--bjob":
        b_job = int(value)

###########################################################
#classifier query sequences by different length
test_500, test_1000,test_3000,test_5000,test_10000,test_genome = [],[],[],[],[],[]

file = open(''+query+'','r')
query_name = []
for line in file:
    line = line.strip('\n')
    if line[0] == '>':#save the query name
        name = line
        query_name.append(line)
        continue
    if len(line) < 1000:
        test_500.append(name+' '+line)
    elif len(line) < 3000:
        test_1000.append(name+' '+line)
    elif len(line) < 5000:
        test_3000.append(name+' '+line)
    elif len(line) < 10000:
        test_5000.append(name+' '+line)
    elif len(line) < 15000:
        test_10000.append(name+' '+line)
    else:
        test_genome.append(name+' '+line)
file.close()

result = pd.DataFrame({"name":query_name})
result.index = result.iloc[:,0]

###########################################################
lengths = [test_500,test_1000,test_3000,test_5000,test_10000,test_genome]
piece = ['500bp','1000bp','3000bp','5000bp','10000bp','genome']
pred_result = []
isFirst = True
for ll in range(len(lengths)):
    if lengths[ll] == []:
        continue

########################################
    data = pd.read_csv('./'+infile+'/'+piece[ll]+'_'+str(kmer)+'_mer',sep = '\t',header = None,index_col = False)
    X = data.iloc[:,1:(data.shape[1] - 1)]#get fearture
    y = data.iloc[:,data.shape[1] - 1]#get label
    del data
    X_train = np.array(X)
    y_train = np.array(y)
    del X,y#clear useless variable, Free memory
    
    #set model parameter
    model = BalancedBaggingClassifier(base_estimator = KNeighborsClassifier(n_neighbors = 1, n_jobs = k_job), n_estimators = 10, n_jobs = b_job)
    #training the model
    model.fit(X_train,y_train)
    
    #clear useless variable, Free memory
    del X_train,y_train
    
    print(''+piece[ll]+' model ready')
#########################################
#predict query sequences

    #In order to reduce memory usage, we only calculate 1000 sequences at a time.
    step = 1000
    sub_query = [lengths[ll][i:i+step] for i in range(0,len(lengths[ll]),step)]
    for sub in sub_query:
        data = queryKmer(sub,kmer)

        #get query sequence fearture
        X_test = np.array(data.iloc[:,1:])
        y_pred = np.array(model.predict(X_test))
        predict_prob_y = np.array(model.predict_proba(X_test)[:,1])
        pred_data = pd.DataFrame(data.iloc[:,0])
        pred_data.index = data.iloc[:,0]
        pred_data.drop([0],axis=1, inplace=True)
        pred_data['Model'] = piece[ll]
        pred_data['Label'] = y_pred
        pred_data['Probability'] = predict_prob_y

        if isFirst:
            isFirst = False
            pred_result = pred_data
        else:#append predict result
            pred_result = pred_result.append(pred_data)

    #clear useless variable, Free memory
    del data, X_test
#print the predict result
result = result.join(pred_result)
result.drop("name",axis=1, inplace=True)
result.index.name = 'Name'
result.to_csv(''+output+'',sep = '\t',header = True,index = True)

end_time = time.time()
total_time = end_time - start_time
print ('Total running time of the program: %.2f seconds' % (total_time))

#clear useless variable, Free memory
del pred_result, result,test_500,test_1000,test_3000,test_5000,test_10000,test_genome,lengths,model