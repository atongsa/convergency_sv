# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:51:04 2022

@author: h
"""

import pandas as pd
import numpy as np
############分析sheep
path_sheep_QTL = r'C:\Users\h\Desktop\random\QTL\sheep\sheep-ram1.0-QTLdb-less5Mb.bed'
sheep_QTL = pd.read_csv(path_sheep_QTL,sep='\t',header=None)
QTLposition = sheep_QTL.iloc[:,0:4]
QTLposition.columns = ['Chr','s','t','qtl']
#导入来自师兄的overlapped文件
path_SV_QTL = r'C:\Users\h\Desktop\random\QTL\sheep\sheep-532-withoutTRA-SV-overlapped-with-less5MbQTL.txt'
SV_QTL = pd.read_csv(path_SV_QTL,sep='\t',header=None)
SV_QTL.columns = ['#CHROM','POS','t','d','e','f','g','QTL_name']
#导入QTL文件
QTL_path = r'C:\Users\h\Desktop\random\QTL\sheep\QTL.xlsx'
QTL_count = pd.read_excel(QTL_path)
#导入overlap到的vcf文件(该文件需要自己从总VCF文件里面提取SV位点得到)
VCF_path = r'C:\Users\h\Desktop\random\QTL\sheep\QTL_SV_overlapped.vcf.recode.vcf'
VCF = pd.read_csv(VCF_path,sep='\t',header=62)
#切片处理   
for i in range(9,len(VCF.columns)):
    VCF.iloc[:,i] = VCF.iloc[:,i].str[:3]
#融合，将VCF的位点信息放进QTL_SV信息中
SV_QTL.columns = ['#CHROM','POS','t','d','e','f','g','QTL_name']
#pd.merge(SV_QTL,VCF,on=['POS','#CHROM']).to_excel(r'C:\Users\h\Desktop\random\QTL\sheep\all_QTL_VCF.xlsx')
#保存重新读取融合QTL_VCF文件
all_QTL_VCF = pd.read_excel(r'C:\Users\h\Desktop\random\QTL\sheep\all_QTL_VCF.xlsx')

QTL_name_merge = all_QTL_VCF.QTL_name.str.split('QTL')
for i in range(len(QTL_name_merge)):
    all_QTL_VCF.iloc[i,6] = QTL_name_merge[i][0] + 'QTL'
    
#计算一下一共有多少个QTL，也就是画图文件的列
x = all_QTL_VCF.iloc[:,6]
plot = pd.DataFrame(x.drop_duplicates())
plot_name = plot.reset_index()['QTL_name']
n_data = pd.DataFrame(np.zeros((147,533)),columns = all_QTL_VCF.columns[7:])
n_data['ID'] = plot_name
#一共交集到了QTL的SV数目是2214 
SV_QTL_count = all_QTL_VCF.POS_s.drop_duplicates()

#计算df文件中一共存在多少个不同种类的基因型
def count_genotype(name,df):
    #output the count of (0/0,0/1,1/1,./.)
    try:
        count_00 = df[name].value_counts()['0/0']
    except KeyError:
        count_00 = 0
    try:
        count_01 = df[name].value_counts()['0/1']
    except KeyError:
        count_01 = 0
    try:
        count_11 = df[name].value_counts()['1/1']
    except KeyError:
        count_11 = 0
    try:
        count_none = df[name].value_counts()['./.']
    except KeyError:
        count_none = 0
    return count_00,count_01,count_11,count_none

#
all_QTL_VCF_drop = all_QTL_VCF.iloc[all_QTL_VCF.POS_s.drop_duplicates().index,:]

#公式(A/B)/(C/D) 开始计算
for i in range(len(n_data)):
    print(i/len(n_data))
    QTL_name = n_data.iloc[i,0]
    targe_SV = all_QTL_VCF.query('QTL_name==@QTL_name').iloc[:,8:]
    This_QTL_SV_count = len(targe_SV)  #所有跟这个QTL有交集的SV数目 *******   C
    all_SV_count = len(all_QTL_VCF.POS_s.drop_duplicates()) #所有跟QTL有交集的SV数目 ******* D
    for j in range(1,len(n_data.columns)):
        name = n_data.columns[j]
        individual_QTL_SV_typecount = count_genotype(name,targe_SV)
        individual_QTL_SV_count = individual_QTL_SV_typecount[1] + individual_QTL_SV_typecount[2]  ##跟该动物在这个QTL上有交集的SV数目 ******* A
        
        individual_allSV_typecount = count_genotype(name,all_QTL_VCF)
        individual_allSV_count = individual_allSV_typecount[1] + individual_allSV_typecount[2]  ##该动物的所有SV数目 ******* B
        if individual_QTL_SV_count == 0:
            n_data.iloc[i,j] = 'NA'
            continue
        
        select = individual_QTL_SV_count/individual_allSV_count
        Bg = This_QTL_SV_count/all_SV_count
        fold_enrichment = select/Bg
        n_data.iloc[i,j] = fold_enrichment    
    
#保存
n_data.to_excel(r'C:\Users\h\Desktop\result_fold_enrichment.xlsx')

n_data = pd.read_excel(r'C:\Users\h\Desktop\result_fold_enrichment.xlsx')
#建立一个对QTL进行分类的dict，用于*纵轴分类
x_path = r'C:\Users\h\Desktop\QLT_count_sort.xlsx'
x = pd.read_excel(x_path,header=None)
x = x.iloc[:,:5]
x.columns = ['Chr','cluster','s','t','name']
dict_name_cluster = x[['cluster','name']]
dict_name_cluster = dict_name_cluster.drop_duplicates()
dict_name = {}
for i in range(len(dict_name_cluster)):
    dict_name[dict_name_cluster.iloc[i,1]] = dict_name_cluster.iloc[i,0]
n_data['qtl_name'] = n_data.iloc[:,0]

#对所有纵轴QTL数据进行cluster分类
for i in range(len(n_data)):
    QTL_name = 'Name=' + n_data.iloc[i,0].split(' QTL')[0]
    cluster = dict_name[QTL_name]
    n_data.iloc[i,0] = cluster
    
id1 = n_data['ID']

#成功对列QTL进行了分类，保存
data = pd.read_excel(r'C:\Users\h\Desktop\result_fold_enrichment_QTLcluster.xlsx').iloc[:,1:]
for i in range(len(id1.str.split('_'))):
    n_data.iloc[i,0] = ' '.join(id1.str.split('_')[i][:-1])
#对所有fold-enrichment值取log2
n_data.iloc[:,1:-1] = np.log2(n_data.iloc[:,1:-1])
#所有Nan的值都用该表格的最小值进行填充
n_data.iloc[:,1:-1] = n_data.iloc[:,1:-1].fillna(n_data.iloc[:,1:-1].min().min())
#保存
n_data.to_excel(r'C:\Users\h\Desktop\result_fold_enrichment_QTLcluster_result.xlsx')
n_data.to_csv(r'C:\Users\h\Desktop\result_fold_enrichment_QTLcluster_result.csv',index=None,sep=' ')
#按照QTL的类别对纵轴进行排列，导出。
n_data.sort_values(by='ID').to_csv(r'C:\Users\h\Desktop\result_fold_sort.csv',sep=' ',index=None)
result_data = n_data.sort_values(by='ID')
x = result_data
oldindex = x.index
x = x.reset_index()
np.mean(x,axis=1).sort_values()
'''
x.query('ID=="Wool"').mean().mean()
Out[65]: -2.112416792969386
x.query('ID=="Reproduction"').mean().mean()
Out[67]: -2.646285368911443
x.query('ID=="Production"').mean().mean()
Out[69]: -2.422326227112234
x.query('ID=="Meat and Carcass"').mean().mean()
Out[70]: -2.434041172843377
x.query('ID=="Health"').mean().mean()
Out[71]: -2.7648975617935623
x.query('ID=="Exterior"').mean().mean()
Out[72]: -2.738802048804467
x.query('ID=="Milk"').mean().mean()
Out[112]: -3.4780209857007196
np.mean(x.query('ID=="Health"'),axis=1).max()
Out[118]: 0.17695985209296453
'''

#接下来对行（个体）进行分类
all_sample_data = pd.read_excel(r'C:\Users\h\Desktop\random\all_sample_data.xlsx')
all_sample_data = all_sample_data.iloc[2:]
sample_loc_dic = {}
for i in range(len(all_sample_data)):
    sample_loc_dic[all_sample_data.iloc[i,2]] = all_sample_data.iloc[i,0]
sample_loc_dic['6']
df = pd.DataFrame(sample_loc_dic,index=sample_loc_dic.keys()).iloc[:1].T
select_name = pd.DataFrame(n_data.columns)
select_name.columns = ['b']
df[1] = df.index
df.columns = ['a','b']
df = df.astype(str)
select_name = select_name.astype(str)
merge = pd.merge(select_name,df,on=1)
#单独导出分类文件
merge.to_csv(r'C:\Users\h\Desktop\name_columns.csv',sep=' ')


#plot
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot
data_path = r'C:\Users\h\Desktop\result_fold_enrichment.xlsx'

data = pd.read_excel(data_path).iloc[:,1:]
pyplot.figure(figsize=(15, 15))
plot=sns.heatmap(data)

data_log2 = np.log2(data)

min1 = data_log2.min().min()
#用整个数据框的最小值去填充缺失值NAN
data_log2.fillna(data_log2.min().min(),inplace=True)
pyplot.figure(figsize=(15, 15))
plot=sns.heatmap(data_log2)


for i in range(len(test)):
    l = test.iloc[i,:].value_counts()
    print(len(l))

############接下来分析goat
#先读取overlap文件，找到position用于vcf的位点提取
import pandas as pd
import numpy as np
############分析goat
path_goat_QTL = r'C:\Users\h\Desktop\random\QTL\goat\goats-ARS1-QTLdb-less5Mb.bed'
goat_QTL = pd.read_csv(path_goat_QTL,sep='\t',header=None)
QTLposition = goat_QTL.iloc[:,0:4]
QTLposition.columns = ['Chr','s','t','qtl']
#导入来自师兄的overlapped文件
#path_SV_QTL = r'C:\Users\h\Desktop\random\QTL\goat\capra-442-withoutTRA-SV-overlapped-with-less5MbQTL.txt'
path_SV_QTL = r'F:\svdata\capra_281_QTL_sv.txt'
SV_QTL = pd.read_csv(path_SV_QTL,sep='\t',header=None)
SV_QTL.columns = ['#CHROM','POS','t','d','e','f','g','QTL_name']
#导入QTL文件
QTL_path = r'C:\Users\h\Desktop\random\QTL\goat\QTL.xlsx'
QTL_count = pd.read_excel(QTL_path,header=None)
#导入overlap到的vcf文件(该文件需要自己从总VCF文件里面提取SV位点得到)
VCF_path = r'C:\Users\h\Desktop\random\QTL\goat\QTL_SV_overlapped.vcf.recode.vcf'
VCF = pd.read_csv(VCF_path,sep='\t',header=2208)
#切片处理
for i in range(9,len(VCF.columns)):
    VCF.iloc[:,i] = VCF.iloc[:,i].str[:3]
#融合，将VCF的位点信息放进QTL_SV信息中
    
#pd.merge(SV_QTL,VCF,on=['POS','#CHROM']).to_excel(r'C:\Users\h\Desktop\random\QTL\goat\all_QTL_VCF.xlsx')
#5 54214537
#3 107932280 #手动矫正

#pd.merge(SV_QTL,capra_281_plink,on=['#CHROM','POS']).to_excel(r'F:\svdata\all_QTL_VCF_20240117.xlsx')

#保存重新读取融合QTL_VCF文件
#删除不需要的列
#all_QTL_VCF = pd.read_excel(r'C:\Users\h\Desktop\random\QTL\goat\all_QTL_VCF.xlsx')
all_QTL_VCF = pd.read_excel(r'F:\svdata\all_QTL_VCF_20240118.xlsx')
#all_QTL_VCF['POS_t'] = all_QTL_VCF['end']
QTL_name_merge = all_QTL_VCF.QTL_name.str.split('QTL')

for i in range(len(QTL_name_merge)):
    all_QTL_VCF.iloc[i,6] = QTL_name_merge.iloc[i][0] + 'QTL'
    
#计算一下一共有多少个QTL，也就是画图文件的列
x = all_QTL_VCF.iloc[:,6]
plot = pd.DataFrame(x.drop_duplicates())
plot_name = plot.reset_index()['ID']
n_data = pd.DataFrame(np.zeros((len(plot),len(all_QTL_VCF.iloc[:,6:].columns))),columns=all_QTL_VCF.columns[6:])
n_data['ID'] = plot_name
#一共交集到了QTL的SV数目是342 
all_QTL_VCF 
'''
SV_QTL.QTL_name.value_counts()
Out[592]: 
Cannon bone circumference QTL (223248)    67
Rump length QTL (223251)                  65
Body width QTL (223246)                   59
Body length QTL (223247)                  56
Rump length QTL (223250)                  47
Withers height QTL (223252)               47
Teat number QTL (223240)                   1
Name: QTL_name, dtype: int64
'''
#公式(A/B)/(C/D) 开始计算
for i in range(len(n_data)):
    print(i/len(n_data))
    QTL_name = n_data.iloc[i,0]
    targe_SV = all_QTL_VCF.query('ID==@QTL_name').iloc[:,7:]
    This_QTL_SV_count = len(targe_SV)  #所有跟这个QTL有交集的SV数目 *******   C
    all_SV_count = len(all_QTL_VCF.POS_s.drop_duplicates()) #所有跟QTL有交集的SV数目 ******* D
    for j in range(1,len(n_data.columns)):
        name = n_data.columns[j]
        individual_QTL_SV_typecount = count_genotype(name,targe_SV)
        individual_QTL_SV_count = individual_QTL_SV_typecount[1] + individual_QTL_SV_typecount[2]  ##跟该动物在这个QTL上有交集的SV数目 ******* A
        
        individual_allSV_typecount = count_genotype(name,all_QTL_VCF)
        individual_allSV_count = individual_allSV_typecount[1] + individual_allSV_typecount[2]  ##该动物的所有SV数目 ******* B
        if individual_QTL_SV_count == 0:
            n_data.iloc[i,j] = 'NA'
            continue
        
        select = individual_QTL_SV_count/individual_allSV_count
        Bg = This_QTL_SV_count/all_SV_count
        fold_enrichment = select/Bg
        n_data.iloc[i,j] = fold_enrichment    
#n_data.to_excel(r'C:\Users\h\Desktop\random\QTL\goat\result_fold_enrichment.xlsx')   
#n_data = pd.read_excel(r'C:\Users\h\Desktop\random\QTL\goat\result_fold_enrichment.xlsx')
#n_data = pd.read_excel(r'F:\svdata\result_fold_enrichment_20240117.xlsx')
n_data = n_data.iloc[:,1:]
n_data.iloc[:,2:] = np.log2(n_data.iloc[:,1:])
#所有Nan的值都用该表格的最小值进行填充
n_data.iloc[:,2:] = n_data.iloc[:,1:].fillna(n_data.iloc[:,1:].min().min())
#保存

n_data.to_excel(r'C:\Users\h\Desktop\random\QTL\goat\result_fold_enrichment_log2_fillna.xlsx')   
n_data.to_csv(r'C:\Users\h\Desktop\random\QTL\goat\result_fold_enrichment_log2_fillna.csv',sep=' ',index=None)   
n_data = pd.read_excel(r'C:\Users\h\Desktop\random\QTL\goat\result_fold_enrichment_log2_fillna.xlsx')
#n_data.to_excel(r'F:\svdata\result_fold_enrichment_log2_fillna_20240117.xlsx')   
#n_data.to_csv(r'F:\svdata\result_fold_enrichment_log2_fillna.csv',sep=' ',index=None)
#n_data = pd.read_excel(r'F:\svdata\result_fold_enrichment_log2_fillna_20240117.xlsx')
#建立一个对QTL进行分类的dict，用于*纵轴分类


x = QTL_count
x.columns = ['Chr','cluster','s','t','id','name']
dict_name_cluster = x[['cluster','name']]
dict_name_cluster = dict_name_cluster.drop_duplicates()
dict_name = {}
for i in range(len(dict_name_cluster)):
    dict_name[dict_name_cluster.iloc[i,1]] = dict_name_cluster.iloc[i,0]

#对横轴的个体进行分类
all_sample_data = pd.read_excel(r'C:\Users\h\Desktop\random\all_sample_data.xlsx')
all_sample_data = all_sample_data.iloc[2:]
sample_loc_dic = {}
for i in range(len(all_sample_data)):
    sample_loc_dic[all_sample_data.iloc[i,3]] = all_sample_data.iloc[i,0]
    

df = pd.DataFrame(sample_loc_dic,index=sample_loc_dic.keys()).iloc[:1].T
select_name = pd.DataFrame(n_data.columns)
select_name.columns = ['b']
df[1] = df.index
df.columns = ['a','b']
df = df.astype(str)
select_name = select_name.astype(str)
merge = pd.merge(select_name,df,on='b')
select_name[select_name.b.str.contains('.',regex=False)]
#单独导出分类文件
merge.to_csv(r'C:\Users\h\Desktop\random\QTL\goat\name_columns.csv',sep=' ')
location = pd.read_csv(r'C:\Users\h\Desktop\random\QTL\goat\name_columns.csv',sep=',')

name_281 = pd.DataFrame(n_data.iloc[:,2:].columns)
name_281.columns = ['Sample ID']






