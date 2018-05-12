#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 18:04:06 2018

@author: DanielCzL315
"""

import pandas as pd
import numpy as np
from scipy import stats
import HC_estimator as hce
import sys
import os

def construct_fishertable(gene1,gene2):
    retval = np.array([[0,0],[0,0]])
    retval[0][0] = gene1.count(1)
    retval[0][1] = gene2.count(1)
    retval[1][0] = gene1.count(0)
    retval[1][1] = gene2.count(0)
    return retval

def construct_vector(file):
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
    vectors = dict()
    
    
    for key in pgdict.keys():
        vector = dict()
        for patient in pset:
            vector[patient] = 0
        for p in pgdict[key]:
            vector[p] = 1
        vectors[key] = vector
               
    vector_list = dict()
    
    for key in vectors.keys():
        curr = vectors[key]
        v = list()
        for p in curr.keys():
            v.append(curr[p])
        vector_list[key] = v
       
    more_than_ten = dict()
    
    for patient in vector_list.keys():
        v = vector_list[patient]
        if v.count(1) >= 10 :
            more_than_ten[patient] = v     
    
    return more_than_ten


def compute_fisher(table):
    #gene1 and gene2 are lists
    return stats.fisher_exact(table,alternative = 'less')

def compute_pearson(gene1,gene2):
    return stats.pearsonr(gene1,gene2)

def compute_HC(gene1,gene2):
    x = np.array(gene1)
    y = np.array(gene2)
    return hce.HC(x,y,0.53)

def paddle(string):
    while(len(string)<6):
        string+=' '
    return string

def paddle_number(number):
    while(len(number)<10):
        number+='0'
    return number

def calc(filename):
        try:
            os.chdir('cancer')
        except FileNotFoundError:
            print('',end='')
        genes_dict = construct_vector(filename)
        genes = list(genes_dict.keys())
        f = 'stat_'+filename
        file = open(f,'w')
        file.write("Gene1  Gene2  oddsratio  pvalue     pearson    p_pval     HC\n")
        for i in range(len(genes)):      
            g1 = genes[i]
            gene1 = genes_dict[g1]
            g1 = paddle(g1)
            for j in range(i+1,len(genes)):
                g2 = genes[j]
                gene2 = genes_dict[g2]
                g2 = paddle(g2)
                table = construct_fishertable(gene1,gene2)
                oddsratio, pval = compute_fisher(table)   
                oddsratio = round(oddsratio,7)
                pval = round(pval,7)
                if pval == 1.0:
                    pval = 0.9999999
                if pval < 0.3:
                    pearson = paddle_number(str(round(compute_pearson(gene1,gene2)[0],7)))
                    pearson_pval = paddle_number(str(round(compute_pearson(gene1,gene2)[1],7)))
                    try:
                        hc = (round(compute_HC(gene1,gene2),7))
                    except ValueError:
                        hc = 'NA        '                   
                else:
                    pearson = 'NA        '
                    pearson_pval = pearson
                    hc = 'NA       '
                entry = g1+" "+g2+" "+paddle_number(str(oddsratio))+" "+paddle_number(str(pval))+" "+str(pearson)+ " " +str(pearson_pval)+" "+str(hc)+"\n"
                file.write(entry)
        file.close()
        print('finished ',filename)


    