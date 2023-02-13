#!/usr/bin/env python
# coding: utf-8

# Copyright (c) 2023 Mitsuru Ohno  
# Use of this source code is governed by a BSD-3-style  
# license that can be found in the LICENSE file.  
#   
# kineticsform: set up reaction kinetics equations from reaction formulas.  


import numpy as np
import pandas as pd


#generate reaction kinetics equations  
def gen_eq(s):
    def m_lst(lst): #extract the chemical species to listing  
        ml = [lst[i+1] for i in range(0, len(lst), 2)]
        return ml  
    
    def dic_lst(lst): #convert list to dictionary, {even;chemical species = key : odd;coefficients = value}
        d = {lst[i+1]:lst[i] for i in range(0, len(lst), 2)}
        try:
            del d['1']
        except:
            pass
        return d  
    
    def rate_eqr(lst): #set up kinetics based on starting species 
        req = ''
        for i, j in enumerate(range(0, len(lst), 2)):
            if i != 0:
                if lst[j] != '1':
                    req = req+'*'+lst[j+1]+'**'+lst[j]
                else:
                    req = req+'*'+lst[j+1]
            else:
                if lst[j] != '1':
                    req = req+lst[j+1]+'**'+lst[j]
                else:
                    req = req+lst[j+1]
        return req
    
    row = []
    react = []
    prod = []
    eq_d = {}
    r = ''
    p = ''
    k = ''
    row = s.values.tolist()
    for i, elem in enumerate(row):
        if elem=='>':
            row.append(i)
    form = row[len(row)-2]
    lat =  row[len(row)-1]
    react=row[3:form]
    prod=row[lat+1:len(row)-2]
    k_n = row[0]
    k_eqr = 'k'+k_n+'*'+rate_eqr(react)
    prod_d = dic_lst(prod)
    react_m = m_lst(react)
    react_m2 = [e for e in react_m if e!='1']
    prod_m = m_lst(prod)
    prod_m2 = [e for e in prod_m if e!='1']
    for e in react_m2:
        eq_d['d'+e] = '-'+k_eqr
    for e in prod_m2:
        if prod_d[e] != '1':
            eq_d['d'+e] = '+(' + k_eqr + ')/' + str(prod_d[e])
        else:
            eq_d['d'+e] = '+' + k_eqr
    return eq_d


def react2kinetic(df):  #main body  
    def k2lst(x):
        temp_l = []
        temp_l = ['k'+str(x[0]), x[1]]
        return temp_l
    marged_eq = []
    k_lst=[]
    dfrev = df.copy()
    k_lst=list(dfrev.apply(k2lst, axis=1))
    dfrev = dfrev.fillna('1')
    dfrev = dfrev.astype('str')
    dfrev.iloc[:,1] = dfrev.iloc[:,1].astype(np.float64)
    dfrev.insert(1, 'eq_dic', np.NaN)
    dfrev['eq_dic'] = dfrev.apply(gen_eq, axis=1)
    uniq_m = set(sum([list(x.keys()) for x in dfrev['eq_dic']],[]))
    uniq_ml = sorted([e for e in list(uniq_m) if e!=None])
    uniq_ml2 = sorted([e[1:] for e in uniq_ml])
    for e in uniq_ml:
        l_eq = [e,]
        for d in dfrev['eq_dic']:
            l_eq.append(d.get(e))
        l_eq2 = [e for e in l_eq if e!=None]
        l_eq2.insert(1, '=')
        l_eq2.append('\n')
        marged_eq.append(''.join(l_eq2))
    uq_mn = 'number of the unique chemical species = '+str(len(uniq_ml))+'\n'+'list of unique chemical species\n'
    marged_eq2 = sorted(marged_eq)
    eq_txt = ''.join(marged_eq2)
    eq_txt = eq_txt.replace('=+', '=')
    kinetics_resl = uq_mn+str(uniq_ml2)+'\n\n'+eq_txt
    return [uniq_ml2, uniq_ml, k_lst, print(kinetics_resl)]

#end

