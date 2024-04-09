# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import random as rd
import numpy as np

#check



path_sheep_ramb1 = ''
path_sheep_QTL = r'C:\Users\h\Desktop\random\sheep-ram1.0-QTLdb-less5Mb.bed'
sheep_QTL = pd.read_csv(path_sheep_QTL,sep='\t',header=None)
QTLposition = sheep_QTL.iloc[:,0:3]

sheep_hot_position = pd.read_csv(r'C:\Users\h\Desktop\random\sheep_hot_position.txt',sep='\t',header=None)
sheep_hot_position_zhu = pd.read_csv(r'C:\Users\h\Desktop\random\sheep-svhotspot-QTL-breakpointCounts_fromzhu.txt',sep='\t',header=None)
sheep_hot_position.columns = ['Chr','s','t']

chr_size = pd.read_csv(r'C:\Users\h\Desktop\random\chrsize1.txt',sep='\t',header=None)

allgenomic = []
for i in range(len(chr_size)):
    if i == 0:
        allgenomic.append([chr_size.iloc[i,0],0,chr_size.iloc[i,1]])
        now = chr_size.iloc[i,1]
        continue
    
    allgenomic.append([chr_size.iloc[i,0],now,now+chr_size.iloc[i,1]])
    now = now + chr_size.iloc[i,1]
allgenomic = pd.DataFrame(allgenomic)
allgenomic.columns = ['Chr','s','t']
 
Range_size = sheep_hot_position['t']-sheep_hot_position['s']

all_cycle_qtl = []
for time in range(1000):
    all_hot = []
    for size in Range_size:
        center = rd.randint(1,2809021901)
        c = allgenomic.query('s<=@center').iloc[-1]
        Chr = c.iloc[0]
        s = c.iloc[1]
        t = c.iloc[2]
        if center-size < s:
            h_size = [Chr,s,s+size]
            all_hot.append(h_size)
            print('碰到边界')
            print(center,size)
            continue
        elif center+size > t:
            h_size = [Chr,t-size,t]
            all_hot.append(h_size)
            print('碰到边界')
            continue
        h_size = [Chr,center-size*0.5-allgenomic.iloc[Chr-1]['s'],center+size*0.5-allgenomic.iloc[Chr-1]['s']]
        #print(Chr,allgenomic.iloc[Chr-1]['s'],center-size*0.5)
        all_hot.append(h_size)
    all_hot_df = pd.DataFrame(all_hot)
    all_hot_df.columns = ['Chr','s','t']
    QTL = QTLposition
    QTL.columns = ['Chr','s','t']
    QTL = QTL.query('Chr!="Chr.27"')
    
    qtl_all = []
    for i in range(len(QTL)):
        qtl = QTL.iloc[i]
    
        qtl_s = qtl[1]
        qtl_t = qtl[2]
        qtl_size = qtl_t-qtl_s
        Chr = qtl.iloc[0]
        #if Chr == "Chr.X":
        #    Chr = "Chr.27"
        #Chr = Chr[4:]
        l = all_hot_df.query('Chr==@Chr')
        #扩大范围，有一半以上交集的都取
        center = qtl_s + qtl_size*0.5
        l['hot_center'] = l['s']*0.5 + l['t']*0.5
        #l['hot_size'] = l['t'] - l['s']
        #l['qtl_size'] = qtl_size
        #l['center_d'] = abs(l['hot_center'] - center)
        #targe_qtl = l.query('(center_d<hot_size*0.5) or (center_d<qtl_size*0.5)')
        targe_qtl = l.query('((s<=@center) and (t>=@center)) or (@qtl_t>=hot_center>=@qtl_s)')
        if len(targe_qtl)==1:
            qtl_all.append([targe_qtl.iloc[0,0],targe_qtl.iloc[0,1],targe_qtl.iloc[0,2]])
        elif len(targe_qtl)>1:
            print('targe_qtl>1')
            qtl_all.append([targe_qtl.iloc[0,0],targe_qtl.iloc[0,1],targe_qtl.iloc[0,2]])
    all_cycle_qtl.append(qtl_all)
    print(time)
    
    
count = []
for i in all_cycle_qtl:
    count.append(len(i))
count = np.array(count)
print(count)
     
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    