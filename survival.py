# -*- coding: utf-8 -*-
"""
Created on Mon May  7 05:14:50 2018

@author: zchen
"""
import pandas as pd
import os
def extract_genes(file):
    Brst = pd.read_table(file,delim_whitespace=True,header = 1,names = ['Patient','Gene'],index_col=False)
    patients = Brst['Patient']
    pset = set(patients)
    
    pgdict = dict()
    f = open(file,'r')
    for line in f:
        newline = line.split()
        try:
            pgdict[newline[1]].append(newline[0])
        except KeyError:
            pgdict[newline[1]] = list()
            pgdict[newline[1]].append(newline[0])
            
    f.close()      
    statfile = 'stat_'+file
    st = pd.read_csv(statfile,delim_whitespace=True)
    temp =st.sort_values('HC',ascending= True)
    pairs = temp.iloc[:50,][['Gene1','Gene2']]
    gset = pairs.apply(tuple,axis = 1).tolist()
    
    return pset,gset


def surv_data(file):
    os.chdir('C:\\Users\\zchen\\Desktop\\research\\cancer')
    print("Working on ", file, "...")
    filename = file+"_MUT.txt"
    survname = file+"_surv.txt_clean"
    df = pd.read_csv(filename,delim_whitespace=True,header = 1,names = ['Patient','Gene'],index_col=False)
    surv = pd.read_csv(survname,delim_whitespace = True).dropna(subset = ["OS_MONTHS"])
    pset,gset = extract_genes(filename)
    dirname  = file+'_survdata'
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    for pair in gset:
        pair_str = str(pair)+'.txt'
        f = open(pair_str,'w')
        for patient in pset:
            patient_gene = set(df.loc[df['Patient']==patient]['Gene'])
            count = 0
            for gene in pair:
                #print(gene)
                if gene in patient_gene:
                    count+=1
            if count == 1:
                label = 1
            else:
                label = 2
            try:
                month = float(surv.loc[surv['ID']==patient]['OS_MONTHS'])
            except TypeError: #NA
                continue
            if month == 0: #dont care about OS_MONTH == 0
                continue
            status = surv.loc[surv['ID']==patient]['OS_STATUS']
            line = str(patient)+" "+str(month)+" "+str(status)+" "+str(label)+"\n"    
            f.write(line)
        f.close()
    os.chdir('C:\\Users\\zchen\\Desktop\\research\\cancer')
    
def main():
    os.chdir('cancer')
    f = open("tasks.txt","r")
    content = f.readlines()
    content = [x.strip() for x in content] 
    f.close()
    for task in content:
        surv_data(task)
    os.chdir('C:\\Users\\zchen\\Desktop\\research')
     
if __name__ == "__main__":
    main()