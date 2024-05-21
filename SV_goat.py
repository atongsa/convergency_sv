# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 16:37:47 2023

@author: h
"""
import numpy as np
import pandas as pd
import re
capra = pd.read_csv(r'F:\svdata\capra_sample281_flt_nohead.vcf',sep='\t')
#capra = pd.read_csv(r'F:\svdata\capra-442-have-sv-change.vcf',sep='\t',header=2208)
sheep = pd.read_csv(r'F:\download\jump_max\sheep-532-have-sv.vcf',sep='\t', header=62)

############************
sheep = pd.read_csv(r'F:\svdata\sheep_532_0516.vcf',sep='\t')
capra = pd.read_csv(r'F:\svdata\capra_281_0516.vcf',sep='\t')

location = pd.read_csv(r'F:\svdata\location.txt',sep=' ',header=None)

'''
indices = capra.loc[capra['INFO'].str.contains('=DEL'), :].index
# 然后，使用这些索引来对原始DataFrame进行赋值
capra.loc[indices, 'type'] = 'DEL'

indices = capra.loc[capra['INFO'].str.contains('=DUP'), :].index
# 然后，使用这些索引来对原始DataFrame进行赋值
capra.loc[indices, 'type'] = 'DUP'

indices = capra.loc[capra['INFO'].str.contains('=INV'), :].index
# 然后，使用这些索引来对原始DataFrame进行赋值
capra.loc[indices, 'type'] = 'INV'

indices = capra.loc[capra['INFO'].str.contains('=TRA'), :].index
# 然后，使用这些索引来对原始DataFrame进行赋值
capra.loc[indices, 'type'] = 'TRA'
indices = capra.loc[capra['INFO'].str.contains('=INS'), :].index
# 然后，使用这些索引来对原始DataFrame进行赋值
capra.loc[indices, 'type'] = 'INS'
'''
'''
INFO = capra[capra['INFO'].str.contains('=TRA')][['#CHROM','POS','INFO']]
#提取TRA信息
# Create a DataFrame
df = INFO

# Function to extract data between given start and end markers
def extract_data(s, start, end):
    pattern = re.compile(r'{}(.*?){}'.format(re.escape(start), re.escape(end)))
    matches = pattern.findall(s)
    return matches[0] if matches else None

# Extract 'CHR2' values
df['CHR2'] = df['INFO'].apply(lambda x: extract_data(x, 'CHR2=', ';'))

# Extract 'END' values, excluding the 'END=' and the trailing ';'
df['END2'] = df['INFO'].apply(lambda x: extract_data(x, 'END=', ';'))
# 提取 merged_df 中的 chr2 和 end2 列，并重置索引以便合并
updated_values = merged_df1[['chr2', 'end2']].reset_index()

# 将 updated_values DataFrame 与 capra_anno_havegene 合并，基于 index 列
# 注意这里使用的是 left join，以保留 capra_anno_havegene 中所有的行
capra_anno_havegene_updated = pd.merge(capra_anno_havegene.reset_index(), updated_values, on='index', how='left', suffixes=('', '_updated'))

# 用 updated_values 中的数据更新 capra_anno_havegene 的对应列
capra_anno_havegene_updated['chr2'] = capra_anno_havegene_updated['chr2_updated'].fillna(capra_anno_havegene_updated['chr2'])
capra_anno_havegene_updated['end2'] = capra_anno_havegene_updated['end2_updated'].fillna(capra_anno_havegene_updated['end2'])

# 删除临时列
capra_anno_havegene_updated.drop(columns=['chr2_updated', 'end2_updated'], inplace=True)

# 恢复原始的 index
capra_anno_havegene_updated.set_index('index', inplace=True)

# 查看更新后的 DataFrame
capra_anno_havegene_updated
'''



#提取breakpoint中的END数据
pattern = r'END=(\d+);'
# 使用 Pandas 的 str.extract 方法提取符合条件的字符串
capra['END_values'] = capra['INFO'].str.extract(pattern)
#
breakpoint_goat = capra[['#CHROM','POS','END_values']]

#capra_281_plink = pd.read_csv(r'F:\svdata\all_test.vcf.vcf',header=582, sep='\t')
capra_281_plink = pd.read_csv(r'F:\svdata\test.vcf',sep='\t',header=560)
#capra_281_input = pd.read_csv(r'F:\svdata\capra_sample281.annovar.input',header=None,sep='\t')

capra_281_plink['type'] = capra['type']


capra_281_plink.insert(2,'end',capra['END_values'])

#SV_type_f = capra_281_plink.iloc[:,3:6]
'''
# 创建第四列 'type'，根据条件填充
def match_types(row):
    types = re.findall(r'(DEL|INV|DUP|BND)', row)
    if not types:  # 如果列表为空，说明未匹配到任何类型
        return 'notag'
    elif len(set(types)) == 1:  # 如果列表中只有一个元素，说明匹配到了一种类型
        return types[0]
    else:
        return 'merge'  # 如果列表中有多个元素，说明匹配到了多种类型

# 对每一行打上SV的type标记
SV_type_f['type'] = 0
SV_type_f['type'] = SV_type_f.apply(lambda row: match_types(' '.join(map(str, row))), axis=1)

#先找TRA和BND
find_TRA = capra[capra.applymap(lambda x: pd.notna(x) and ('BND:' in str(x) or 'TRA:' in str(x)))]
TRA = find_TRA.dropna(how='all')
x = capra.iloc[TRA.index]
SV_type_f.iloc[TRA.index,-1] = 'BND'
#其他的notag为INV
SV_type_f.iloc[SV_type_f.query('type=="notag"').index,-1]='INV_1'

# 定义一个函数来从字符串中提取目标字符,获取INV的END
def extract_end(row):
    pattern = r'\[(.*?)\[|\](.*?)\]'
    ref_match = re.search(pattern, ['REF'])
    alt_match = re.search(pattern, row['ALT'])

    if ref_match and alt_match:
        return ''.join(part for part in ref_match.groups() + alt_match.groups() if part is not None)
    elif ref_match:
        return ''.join(part for part in ref_match.groups() if part is not None)
    elif alt_match:
        return ''.join(part for part in alt_match.groups() if part is not None)
    else:
        return None
INV_END = SV_type_f.query('type=="INV_1"').apply(extract_end, axis=1)
capra_281_plink['type'] = SV_type_f['type']
capra_281_plink.iloc[capra_281_plink.query('type=="INV_1"').index,2] =INV_END

BND_CHRB_END = SV_type_f.query('type=="BND"').apply(extract_end, axis=1)
capra_281_plink.iloc[capra_281_plink.query('type=="BND"').index,2] = BND_CHRB_END
#capra_281_plink中的end列的每一项如果出现有：，则将其值替换为：以后的字符
capra_281_plink.end = capra_281_plink.end.astype(str)
def replace_colon(value):
    """
    

    Parameters
    ----------
    value : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if ':' in value:
        return value.split(':', 1)[1]  # 将冒号后的部分提取出来
    else:
        return value
capra_281_plink['end'] = capra_281_plink['end'].apply(replace_colon)


'''
#capra_281_plink.insert(3,'len',capra_281_plink['end']-capra_281_plink['POS'])
#开始进行数目统计
#求长度diff
capra_281_plink_data = capra_281_plink.iloc[:,-282:-1]
#capra_281_plink.iloc[:,10:-3]

'''
#重新将分配给BND的INV找出来
find_INV = capra[capra.applymap(lambda x: pd.notna(x) and ('INV:' in str(x)))]
INV = find_INV.dropna(how='all')
capra_281_plink.iloc[INV.index,-1] = "INV"
'''
#将无法计算长度的TRA长度值赋予NA
capra_281_plink['diff'] = pd.to_numeric(capra_281_plink['end'], errors='coerce', downcast='float') - pd.to_numeric(capra_281_plink['POS'], errors='coerce', downcast='float')

'''
capra_281_plink_data['type'] = capra_281_plink['type']
capra_281_plink_data = capra_281_plink_data.replace('INV_1','INV')
capra_281_plink = capra_281_plink.replace('INV_1','INV')
'''
#sheep_532_plink_0110.to_csv(r'F:\svdata\sheep_532_plink_0110.vcf',sep='\t',index=None)

'''
#删除出现了重复的反向的inversion
capra_281_plink[capra_281_plink['diff'] <-1]
#先计算反向INV中end与其他所有POS相等的值，然后用end在反向INV的POS中验证重复的INV
same = capra_281_plink[capra_281_plink['POS'].isin(capra_281_plink[capra_281_plink['diff'] <-1]['end'].astype(int))]
s_f1 = capra_281_plink[capra_281_plink['diff'] <-1]['POS']
dup_INV = s_f1[s_f1.astype(int).isin(same['end'].astype(int))].index
#删除重复INV
capra_281_plink = capra_281_plink.drop(dup_INV)
#考虑删除merge了多种SV的行87069 --> 87040
capra_281_plink = capra_281_plink.query('type!="merge"') 
'''



'''
处理INV反向重复
——————————————————————————————————————————————————————————
'''
capra_281_plink.iloc[:,1] = capra_281_plink.iloc[:,1].astype(int)
capra_281_plink.iloc[:,2] = capra_281_plink.iloc[:,2].astype(int)
checked_pairs = {}

# 循环遍历 DataFrame 数据，检查反向重复且在同一染色体上
for idx, row in capra_281_plink.iterrows():
    pos = row['POS']
    end = row['end']
    chrom = row['#CHROM']

    # 构建反向键
    reverse_key = (chrom, end, pos)

    # 检查是否已存在反向键，并确认是否在同一条染色体上
    if reverse_key in checked_pairs:
        # 仅添加后出现的行索引到删除列表
        checked_pairs[reverse_key]['duplicates'].append(idx)
    else:
        # 否则，创建新键，记录首次出现的索引和后续重复的索引
        checked_pairs[(chrom, pos, end)] = {'first': idx, 'duplicates': []}

# 获取所有重复的索引，仅包括后出现的索引
to_remove = [index for pair in checked_pairs.values() for index in pair['duplicates']]
# 删除这些索引的行
capra_281_plink = capra_281_plink.drop(to_remove)
#capra = capra.drop(to_remove)
'''
处理TRA重复
——————————————————————————————————————————————————————————
'''
# 筛选出 type 为 "TRA" 的数据
tra_data = capra_281_plink[capra_281_plink['type'] == 'TRA']

# 查找具有反向重复的行的索引
# 使用字典存储检查过的组合，避免重复检查
checked_pairs = {}

# 循环遍历 DataFrame 数据，检查反向重复
for idx, row in tra_data.iterrows():
    pos = row['POS']
    end = row['end']
    
    # 构建反向键
    reverse_key = (end, pos)
    
    # 检查是否已存在反向键
    if reverse_key in checked_pairs:
        # 仅添加后出现的行索引到删除列表
        checked_pairs[reverse_key]['duplicates'].append(idx)
    else:
        # 否则，创建新键，记录首次出现的索引和后续重复的索引
        checked_pairs[(pos, end)] = {'first': idx, 'duplicates': []}

# 获取所有重复的索引，仅包括后出现的索引
to_remove = [index for pair in checked_pairs.values() for index in pair['duplicates']]
# 删除这些索引的行
capra_281_plink_clean = capra_281_plink.drop(to_remove)
#capra = capra.drop(to_remove)
#————————————————————————————————————————————————————————————————————————————————————————————
capra_281_plink = capra_281_plink_clean
#diff相减后要+1才是实际值
capra_281_plink['diff'] = capra_281_plink['diff'].abs()
capra_281_plink['diff'] = capra_281_plink['diff'] + 1
capra_281_plink = capra_281_plink.reset_index(drop=True)
capra_281_plink.iloc[capra_281_plink.query('type=="TRA"').index,-1]=None

'''###############################################
capra_281_plink.to_csv(r'F:\svdata\capra_281_plink_clean.vcf', index=None, sep='\t')
capra_281_plink = pd.read_csv(r'F:\svdata\capra_281_plink_clean.vcf',sep='\t')
capra_281_plink.type.value_counts()
Out[1058]: 
DEL    65731
TRA    11372
DUP     6596
INV     2583
INS        1
Name: type, dtype: int64
sheep_532_plink_0110.to_csv(r'F:\svdata\sheep_532_plink_0110_clean.vcf', index=None, sep='\t')
sheep_532_plink_0110.type.value_counts()
'''
#capra.iloc[capra_281_plink.index,:-1].to_csv(r'F:\svdata\capra_svout_87404_nohead.vcf',index=None,sep='\t')

capra_281_plink_data = capra_281_plink.filter(regex='RR')
capra_281_plink_data_except_type = capra_281_plink_data.iloc[:, :]
# 找出元素为 "1/1" 或 "0/1" 的行
filtered_rows = capra_281_plink_data_except_type[(capra_281_plink_data_except_type == '1/1') | (capra_281_plink_data_except_type == '0/1')]
# 合并 type 列
capra_281_plink_data['type'] = capra_281_plink['type']
merged_type = pd.concat([filtered_rows, capra_281_plink_data['type'], capra_281_plink['diff']], axis=1)
# 遍历每一列，统计非缺失值的 type 列的出现次数

all_SVtype_data = []
for col in merged_type.columns[:-2]:  # 排除最后一列 type
    ind_sv =  merged_type[[col, 'type']].dropna()['type']
    non_na_type_counts = ind_sv.value_counts()
    DEL_len = merged_type.iloc[ind_sv[ind_sv == 'DEL'].index]['diff'].abs().sum()
    DUP_len = merged_type.iloc[ind_sv[ind_sv == 'DUP'].index]['diff'].abs().sum()
    INV_len = merged_type.iloc[ind_sv[ind_sv == 'INV'].index]['diff'].abs().sum()
    INS_len = merged_type.iloc[ind_sv[ind_sv == 'INS'].index]['diff'].abs().sum()
    if INS_len !=0:
        all_SVtype_data.append([merged_type[[col]].columns[0], non_na_type_counts.DEL, non_na_type_counts.DUP, non_na_type_counts.TRA, non_na_type_counts.INV, non_na_type_counts.INS, DEL_len, DUP_len, INV_len, INS_len])
    else:
        all_SVtype_data.append([merged_type[[col]].columns[0], non_na_type_counts.DEL, non_na_type_counts.DUP, non_na_type_counts.TRA, non_na_type_counts.INV, 0, DEL_len, DUP_len, INV_len, 0])
    #print(f"Counts for {col}:")
    #print(non_na_type_counts)
    #print()

SVtype_data = pd.DataFrame(all_SVtype_data)
SVtype_data.columns = ['name', 'DEL', 'DUP', 'TRA', 'INV', 'INS', 'DEL_len', 'DUP_len', 'INV_len', 'INS_len']
#将vcf文件的sample中.1，.2这种删掉
#SVtype_data.iloc[SVtype_data[SVtype_data.name.str.contains('\.')].index,0] = [i[:-2] for i in SVtype_data[SVtype_data.name.str.contains('\.')].name]

SVtype_data['all_number'] = SVtype_data.DEL + SVtype_data.TRA + SVtype_data.DUP + SVtype_data.INV + SVtype_data.INS

SVtype_data['all_len'] = SVtype_data.DEL_len + SVtype_data.DUP_len + SVtype_data.INV_len + SVtype_data.INS_len

#根据名称信息将SRR转为对应的物种名
species = pd.read_csv(r'C:\Users\h\Desktop\SVdata\Species.txt',sep='\t')
species = species.iloc[:,[0,5,-1]]
species.columns = ['Species','code','name']
SV_result = pd.merge(species, SVtype_data, on='name')
#SVtype_data.to_excel(r'C:\Users\h\Desktop\GB图片修改\Table6.xlsx')
#SVtype_data.to_excel(r'C:\Users\h\Desktop\GB图片修改\Table6_sheep.xlsx') #20240514修改
SV_result.to_excel(r'C:\Users\h\Desktop\GB图片修改\Table6.xlsx')
#最终结果
capra_281_plink.to_csv(r'F:\svdata\SV_count_result.csv',index=None)

#统计table S7
#50-100bp	100bp-250bp	250bp-500bp	500bp-1kb	1kb-2kb	2kb-5kb	5kb-10kb	10kb-50kb	50kb-100kb	100kb-500kb	500kb-1Mb	Total
def count_length(SV):
    a = len(capra_281_plink.query('type==@SV and 50<diff<=100'))
    b = len(capra_281_plink.query('type==@SV and 100<diff<=250'))
    c = len(capra_281_plink.query('type==@SV and 250<diff<=500'))
    d = len(capra_281_plink.query('type==@SV and 500<diff<=1000'))
    e = len(capra_281_plink.query('type==@SV and 1000<diff<=2000'))
    f = len(capra_281_plink.query('type==@SV and 2000<diff<=5000'))
    g = len(capra_281_plink.query('type==@SV and 5000<diff<=10000'))
    h = len(capra_281_plink.query('type==@SV and 10000<diff<=50000'))
    i = len(capra_281_plink.query('type==@SV and 50000<diff<=100000'))
    j = len(capra_281_plink.query('type==@SV and 100000<diff<=500000'))
    k = len(capra_281_plink.query('type==@SV and 500000<diff<=1000000'))
    return [a,b,c,d,e,f,g,h,i,j,k]

capra_281_plink['diff'] = capra_281_plink['diff'].abs()
all_count_length = []
all_count_length.append(count_length('DEL'))
all_count_length.append(count_length('DUP'))
all_count_length.append(count_length('INV'))
all_count_length.append(count_length('INS'))
pd.DataFrame(all_count_length)
#capra_281_plink[['#CHROM','POS','end','type']].query('type=="DEL"').to_csv(r'F:\svdata\capra_281_DEL.bed',index=None,header=None,sep='\t')
'''
Out[411]: 
      0      1      2     3     4     5     6    7   8   9   10
0  20053  16619  10638  8010  5850  2933  1002  462  77  72  15
1    766   3078   1130   661   533   187    92  112  24  12   1
2     80    227    350   737   674   271   106   86  11  29  11
3      1      0      0     0     0     0     0    0   0   0   0

'''

"____________________________________________________________________________"
"注释结果与基因数目和分布统计_______________________________________________________________________"
capra_anno_input


#Capra
#capra_anno = pd.read_csv(r'F:\svdata\capra_87404.sv.vcf.anno.variant_function',sep='\t',header=None)
#capra_anno = pd.read_csv(r'F:\svdata\capra_87404.type.anno.variant_function',sep='\t',header=None)
capra_anno = pd.read_csv(r'F:\svdata\capra_281_0515.variant_function',sep='\t',header=None)
#sheep_anno = pd.read_csv(r'F:\svdata\sheep_532_type.vcf.variant_function',sep='\t',header=None)
#capra_281_0515.anno
#capra_anno = capra_anno.iloc[:,3:]
capra_anno.columns = [0,1,2,3,4,5,6,7,8,9,10,11]
capra_anno_havegene = capra_anno[capra_anno[0] != "intergenic"]
"""
capra_anno[0].value_counts()
Out[1465]: 
intergenic        49629
intronic          30384
upstream           1908
downstream         1808
exonic             1642
splicing           1073
ncRNA_intronic      846
UTR3                682
UTR5                375
ncRNA_exonic        291
ncRNA_splicing       31
Name: 0, dtype: int64
"""

#####################################################################################
# Apply the function to each row and expand any rows with multiple genes into separate rows
expanded_rows = []
for idx, row in capra_anno_havegene.iterrows():
    genes = extract_genes(row[1])
    for gene in genes:
        new_row = row.copy()
        new_row['gene'] = gene
        expanded_rows.append(new_row)
#capra_anno_havegene
# Create a new DataFrame from the expanded rows
expanded_df = pd.DataFrame(expanded_rows).reset_index(drop=True)

expanded_df




count_gene_occurrences = capra_anno_havegene[1].str.contains('gene-').sum()
#用正则表达式提取所有基因
pattern = 'gene-(.*?)(?=,|$|\()'
#提取基因
all_gene_location = capra_anno_havegene[1].str.extractall(pattern).astype(str)
all_gene_loc


INFO = capra[capra['INFO'].str.contains('=TRA')][['#CHROM','POS','INFO']]
df = INFO
df['CHR2'] = df['INFO'].apply(lambda x: extract_data(x, 'CHR2=', ';'))
df['END2'] = df['INFO'].apply(lambda x: extract_data(x, 'END=', ';'))
capra['chr2'] = capra['#CHROM'] 
capra.loc[capra['type'] == "TRA", 'chr2'] = df['CHR2'].values



"____________________________________________________________________________"
"计算有多少个不重复基因，基因被外显子注释、内含子注释、还是都有。_______________________________________________________________________"
#其中一共有12908个基因被SV注释到，其中1725个基因被外显子注释到，454个基因只被外显子注释到。
#df_exploded[df_exploded[0]=='exonic'][1].drop_duplicates()
#df_exploded[df_exploded[0]=='intronic'][1].drop_duplicates()
import pandas as pd

def clean_gene_name(gene):
    # 检查'gene-'前缀是否存在
    if 'gene-' in gene:
        # 如果存在，移除前缀和括号后的内容
        return gene.strip().split('gene-')[1].split('(')[0]
    else:
        # 如果'gene-'前缀不存在，返回原始字符串（或可以选择返回None或特定的占位符）
        return gene.strip()


"""
#capra_anno_havegene[[0, 10]] 计算DEL的不同类型,计算正文表2
#gene_SV_type = capra_anno_havegene[[0, 10]]
#gene_SV_type.columns = ['loc', 'type']
#pd.DataFrame(gene_SV_type.query('type=="DEL"')['loc']).value_counts()
import pandas as pd

# 假设 df 和 common_gene 已经被正确加载
import pandas as pd
import re

# 假设 df 和 common_gene 已经被正确加载

# 初始化一个空的DataFrame用于存放结果
result_df = pd.DataFrame(columns=df.columns)

# 更新正则表达式以匹配包括破折号和数字在内的基因名称
gene_pattern = r'-(\w[\w-]*\w)(?=[(,]|$)'

# 遍历df的每一行
for index, row in df.iterrows():
    # 使用正则表达式找出所有满足条件的基因名称
    genes = re.findall(gene_pattern, row[1])
    
    # 遍历找到的基因名称
    for gene in genes:
        # 如果基因名称在common_gene中，则创建一个新的行添加到result_df中
        if gene in common_gene[0].values:
            new_row = row.copy()
            new_row[1] = gene  # 将第二列设置为当前的基因名称
            result_df = result_df.append(new_row, ignore_index=True)

print(result_df)
result_df.query('type=="DUP"')['loc'].value_counts()

"""
'''
# Assuming capra_anno_havegene is similar to the sample_df defined earlier, with an additional requirement for handling multiple genes

def extract_genes(s):
    pattern = re.compile(r'gene-([^\(,]+)')
    return pattern.findall(s)



# Apply the function to each row and expand any rows with multiple genes into separate rows
expanded_rows = []
for idx, row in sample_df.iterrows():
    genes = extract_genes(row[1])
    for gene in genes:
        new_row = row.copy()
        new_row['gene'] = gene
        expanded_rows.append(new_row)
#capra_anno_havegene
# Create a new DataFrame from the expanded rows
expanded_df = pd.DataFrame(expanded_rows).reset_index(drop=True)

expanded_df


import pandas as pd

# 假设 capra_ref 已经加载到 DataFrame 中
# 示例提取基因名、染色体、起始位置和结束位置
capra_ref = capra_ref[capra_ref[2] == 'gene']
capra_ref['gene'] = capra_ref[8].str.extract('Name=([^;]+)')
capra_ref['chromosome'] = capra_ref[0]
capra_ref['start_position'] = capra_ref[3]
capra_ref['end_position'] = capra_ref[4]

# 现在创建一个新的 DataFrame 仅包含所需列
genes_info = capra_ref[['gene', 'chromosome', 'start_position', 'end_position']].dropna()



expanded_df_merge = pd.merge(expanded_df, genes_info, on='gene', how='left')
expanded_df_merge.to_excel(r'F:\svdata\S14.xlsx')

# 示例提取基因名、染色体、起始位置和结束位置

sheep_ref = pd.read_csv(r'F:\svdata\Oar_rambouillet_v1.0x.gff',sep='\t',comment='#',header=None)
sheep_ref = sheep_ref[sheep_ref[2] == 'gene']
sheep_ref['gene'] = sheep_ref[8].str.extract('Name=([^;]+)')
sheep_ref['chromosome'] = sheep_ref[0]
sheep_ref['start_position'] = sheep_ref[3]
sheep_ref['end_position'] = sheep_ref[4]

# 现在创建一个新的 DataFrame 仅包含所需列
genes_info = sheep_ref[['gene', 'chromosome', 'start_position', 'end_position']].dropna()



expanded_df_merge = pd.merge(expanded_df, genes_info, on='gene', how='left')
expanded_df_merge.to_excel(r'F:\svdata\S15.xlsx')





result = expanded_df[[0,2,3,'chr2','end2','gene']]
gff = pd.read_csv(r'F:\svdata\extracted_genes.txt',header=None,sep=' ')
#转换GFF
# Adjusted function for replacing chromosome identifiers
def replace_chromosome_v2(chromosome):
    prefix = "NC_0308"
    if chromosome.startswith(prefix):
        try:
            # Extract the numerical part after the prefix and convert to integer
            num = int(chromosome[len(prefix):].split('.')[0])
            # Check if the number is within the desired range
            if 8 <= num <= 36:
                return str(num - 7)  # Adjust the range to start from 1
        except ValueError:
            # In case the conversion to integer fails, return the original chromosome value
            return chromosome
    return chromosome

# Apply the adjusted function to the first column of the DataFrame
gff[0] = gff[0].apply(replace_chromosome_v2)
gff.columns = ['a','b','c','gene']
merge = pd.merge(expanded_df, gff, on='gene')
pd.merge(expanded_df, gff, on='gene').to_excel(r'F:\svdata\S14.xlsx')


# Adding new columns based on the conditions described by the user
merged_df['new_col_goat'] = np.where(merged_df['Chromosome_goat'] == '', '', merged_df['gene'])
merged_df['new_col_sheep'] = np.where(merged_df['Chromosome_sheep'] == '', '', merged_df['gene'])


######S15
sheep_result = pd.read_excel(r'C:\Users\h\Desktop\GB图片修改\Additional file 2.xlsx',sheet_name="Table S13 √",header=1)
goat_result = pd.read_excel(r'C:\Users\h\Desktop\GB图片修改\Additional file 2.xlsx',sheet_name="Table S14 √", header=1)

common_goat = pd.merge(common_genes, goat_result, on=['Gene ID'])
common_sheep = pd.merge(common_genes, sheep_result, on=['Gene ID'])
common_goat.columns = ['gene', 'Chromosome', 'Gene start', 'Gene end', 'SV chrA', 'SV start', 'SV chrB', 'SV end', 'SV type', 'Annotation']
common_sheep.columns = ['gene', 'Chromosome', 'Gene start', 'Gene end', 'SV chrA', 'SV start', 'SV chrB', 'SV end', 'SV type', 'Annotation']
# 为每个基因分配一个顺序号
common_sheep['order'] = common_sheep.groupby('gene').cumcount()
common_goat['order'] = common_goat.groupby('gene').cumcount()

# 找到每个基因在两个 DataFrame 中的最大出现次数
max_order_sheep = common_sheep.groupby('gene')['order'].max()
max_order_goat = common_goat.groupby('gene')['order'].max()
max_order = pd.concat([max_order_sheep, max_order_goat], axis=1).max(axis=1)

# 创建一个以基因和顺序号为组合键的 DataFrame
sheep_expanded = common_sheep.set_index(['gene', 'order']).reindex(
    pd.MultiIndex.from_tuples(
        [(gene, order) for gene in max_order.index for order in range(max_order[gene] + 1)],
        names=['gene', 'order']
    )
).reset_index()

goat_expanded = common_goat.set_index(['gene', 'order']).reindex(
    pd.MultiIndex.from_tuples(
        [(gene, order) for gene in max_order.index for order in range(max_order[gene] + 1)],
        names=['gene', 'order']
    )
).reset_index()

# 外部合并扩展后的 DataFrame
merged_df = pd.merge(sheep_expanded, goat_expanded, on=['gene', 'order'], how='outer', suffixes=('_sheep', '_goat'))

# 填充 NaN 值为空白
merged_df.fillna('', inplace=True)

# 删除辅助列 'order'
merged_df.drop('order', axis=1, inplace=True)
'''


# 应用更健壮的处理方法
df = capra_anno_havegene.iloc[:,:2]
split_genes = df[1].str.split(',').explode()
cleaned_genes = split_genes.apply(clean_gene_name)

# 重新统计不同基因的总数
unique_genes_count = cleaned_genes.nunique()

# 重新构建DataFrame，以便于筛选包含'exonic'的基因
df_exploded = pd.DataFrame({0: df[0].repeat(split_genes.groupby(split_genes.index).size()), 1: cleaned_genes})

# 筛选包含'exonic'的行，并提取基因名称
exonic_genes = df_exploded[1][df_exploded[0].str.contains("exonic")].unique()

# 重新计算包含'exonic'的不同基因数量
exonic_genes_count = len(exonic_genes)

unique_genes_count, exonic_genes_count

# 找到所有出现过的基因
all_genes = cleaned_genes.unique()

#goat_gene = all_genes
#sheep_gene = all_genes


# 初始化一个列表来保存只有'exonic'没有其他位置信息的基因
exonic_only_genes = []

# 对于每个基因，检查它是否只与'exonic'相关联
for gene in all_genes:
    # 找到这个基因的所有位置信息
    gene_positions = df_exploded[df_exploded[1] == gene][0].unique()
    
    # 如果这个基因的位置信息只有'exonic'，则添加到列表中
    if len(gene_positions) == 1 and 'exonic' in gene_positions:
        exonic_only_genes.append(gene)

# 计算只有'exonic'没有其他位置信息的基因数量
exonic_only_genes_count = len(exonic_only_genes)

exonic_only_genes_count

"____________________________________________________________________________"

#去重
all_gene = all_gene_location.drop_duplicates()
#all_gene = all_genes
result_str = ','.join(all_gene[0])
gene_df = pd.DataFrame(result_str.split(','))
#删除LOC开头的基因
gene_df = gene_df[~gene_df[0].str.startswith('LOC')]
#gene_df_sheep = pd.DataFrame(sheep_gene)
#gene_df_goat = pd.DataFrame(goat_gene)
#Sheep
sheep_anno = pd.read_csv(r'F:\svdata\sheep_532.type.anno1.variant_function',sep='\t',header=None)
#sheep_anno = pd.read_csv(r'F:\svdata\sheep-532-have-sv.annovar.variant_function',sep='\t',header=None)
sheep_anno_havegene = sheep_anno[sheep_anno[0] != "intergenic"]
all_gene_location = sheep_anno_havegene[1].str.extractall(pattern).astype(str)
"""
extracted_genes = sheep_anno_havegene[1].str.extractall(pattern).astype(str)

# 通过索引匹配将提取后的基因名称与原 DataFrame 合并
result = sheep_anno_havegene.loc[extracted_genes.index.get_level_values(0)].copy()
result['gene'] = extracted_genes[0].values
"""


all_gene = all_gene_location.drop_duplicates()
result_str = ','.join(all_gene[0])
gene_df_sheep = pd.DataFrame(result_str.split(','))
gene_df_sheep = gene_df_sheep[~gene_df[0].str.startswith('LOC')]
#计算SV注释到的merge基因
sheep_goat_merge = pd.merge(gene_df_sheep, gene_df, on='gene')

#计算山羊和绵羊SV注释到基因的交集 5948
pd.merge(gene_df, gene_df_sheep)
#common_gene.sort_values(by=[0])


"____________________________________________________________________________"
"Figure3A计算热点区域中有多少基因绘图 绘图网站：http://mg2c.iask.in/mg2c_v2.1/_______________________________________________________________________"

goat_hotspot_anno = pd.read_csv(r'F:\svdata\capra_hotspot_sv_anno.variant_function',sep='\t',header=None)
capra_anno_havegene = goat_hotspot_anno[goat_hotspot_anno[3] != "intergenic"]
count_gene_occurrences = capra_anno_havegene[4].str.contains('gene-').sum()
#用正则表达式提取所有基因
pattern = 'gene-(.*?)(?=,|$|\()'
#提取基因
all_gene_location = capra_anno_havegene[4].str.extractall(pattern).astype(str)
#去重
capra_all_gene = all_gene_location.drop_duplicates()
############ goat 1516个
sheep_hotspot_anno = pd.read_csv(r'F:\svdata\sheep_hotspot_sv_anno.variant_function',sep='\t',header=None)
sheep_anno_havegene = sheep_hotspot_anno[sheep_hotspot_anno[3] != "intergenic"]
count_gene_occurrences = sheep_anno_havegene[4].str.contains('gene-').sum()
#用正则表达式提取所有基因
pattern = 'gene-(.*?)(?=,|$|\()'
#提取基因
all_gene_location = sheep_anno_havegene[4].str.extractall(pattern).astype(str)
#去重
sheep_all_gene = all_gene_location.drop_duplicates()


#df_combined = df.groupby(['chr', 's', 't'])['gene'].apply(','.join).reset_index()



"____________________________________________________________________________"
"Figure3a数据计算,绘图文件准备_______________________________________________________________________"
#输出基因交集，用于绘制figure 3a
pd.merge(capra_all_gene,sheep_all_gene).to_csv(r'F:\svdata\hotspot_gene_goat-sheep.txt', index=None,header=None,sep='\t')
goat_sheep_hotspotgene = pd.merge(capra_all_gene,sheep_all_gene)
#输出sheep的热点区域与其注释基因的交集，在此处读取
sheep_hotspot_gene_old = pd.read_csv(r'F:\svdata\sheep_hotspot_sv.txt',sep='\t',header=None)
sheep_hotspot_gene = sheep_hotspot_gene_old.iloc[:,:7]
sheep_hotspot_gene = sheep_hotspot_gene[sheep_hotspot_gene[6] != "intergenic"]
pattern = 'gene-(.*?)(?=,|$|\()'
#提取基因
sheep_hotspot_gene = sheep_hotspot_gene.loc[sheep_hotspot_gene[7].drop_duplicates().index]
sheep_hotspot_gene[7] = sheep_hotspot_gene[7].replace(to_replace=r'gene-|\([^)]*\)', value='', regex=True)
sheep_hotspot_gene.columns = ['chr1','s1','t1','chr2','s2','t2','type','gene']
#去重
sheep_hotspot_gene = sheep_hotspot_gene.loc[sheep_hotspot_gene.gene.drop_duplicates().index]
goat_sheep_hotspotgene
goat_sheep_hotspotgene.columns = ['gene']
all_merge = pd.merge(sheep_hotspot_gene, goat_sheep_hotspotgene,on='gene')
all_merge.to_excel(r'F:\svdata\sheep_goat_hotspotgene.xlsx')
#进行处理#GPR146
all_merge = pd.read_excel(r'F:\svdata\sheep_goat_hotspotgene.xlsx')
#合并行,准备作为绘图数据
all_merge_plot = all_merge.groupby(['chr1', 's1', 't1'])['gene'].apply(','.join).reset_index()

"____________________________________________________________________________"
"Figure2a数据计算_______________________________________________________________________"
type1 = pd.read_csv(r'F:\svdata\type.txt','\t', header=None)
#type1 = pd.read_csv(r'F:\svdata\sheep_type1.txt','\t', header=None)
"""
######################对sheep的SV进行分类
find_DEL = sheep.loc[sheep['INFO'].str.contains('=DEL'), :].index
find_INV = sheep.loc[sheep['INFO'].str.contains('=INV'), :].index
find_DUP = sheep.loc[sheep['INFO'].str.contains('=DUP'), :].index
find_TRA = sheep.loc[sheep['INFO'].str.contains('=TRA'), :].index
find_INS = sheep.loc[sheep['INFO'].str.contains('=INS'), :].index
sheep_532_plink_0110.loc[find_DEL, 'type'] = 'DEL'
sheep_532_plink_0110.loc[find_INV, 'type'] = 'INV'
sheep_532_plink_0110.loc[find_DUP, 'type'] = 'DUP'
sheep_532_plink_0110.loc[find_TRA, 'type'] = 'TRA'
sheep_532_plink_0110.loc[find_INS, 'type'] = 'INS'
#sheep_532_plink_0110.to_csv(r'F:\svdata\sheep_532_plink_0110.vcf',sep='\t',index=None)

"""
# 从type1中获取Wild goat和Native goat的ID列表
wild_goat_ids = type1[type1[0] == 'Wild goat'][1].tolist()
native_goat_ids = type1[type1[0] == 'Native goat'][1].tolist()
improved_goat_ids = type1[type1[0] == 'Improved goat'][1].tolist()
# 从capra_281_plink中筛选出Wild goat和Native goat个体对应的列
wild_goat_df = capra_281_plink[wild_goat_ids]
native_goat_df = capra_281_plink[native_goat_ids]
improved_goat_df = capra_281_plink[improved_goat_ids]

wild_goat_df[['type','diff']] = capra_281_plink[['type','diff']]
native_goat_df[['type','diff']] = capra_281_plink[['type','diff']]
improved_goat_df[['type','diff']] = capra_281_plink[['type','diff']]
# 检查DataFrame中每行是否包含"0/1"或"1/1"
wild_goat_df_filtered = wild_goat_df[wild_goat_df.apply(lambda row: ('0/1' in row.values) or ('1/1' in row.values), axis=1)]
# 检查DataFrame中每行是否包含"0/1"或"1/1"
native_goat_df_filtered = native_goat_df[native_goat_df.apply(lambda row: ('0/1' in row.values) or ('1/1' in row.values), axis=1)]
improved_goat_df_filtered = improved_goat_df[improved_goat_df.apply(lambda row: ('0/1' in row.values) or ('1/1' in row.values), axis=1)]

# 获取每个DataFrame的索引集合
native_indices = set(native_goat_df_filtered.index)
wild_indices = set(wild_goat_df_filtered.index)
improved_indices = set(improved_goat_df_filtered.index)
# 找出只在native_goat_df_filtered中存在的索引
len(wild_indices - improved_indices - native_indices)
len(native_indices - wild_indices - improved_indices)
len(improved_indices - native_indices -wild_indices)
# 计算三个集合的交集
len(native_indices.intersection(wild_indices, improved_indices))
# 计算两两交集
len(native_indices.intersection(wild_indices))
len(native_indices.intersection(improved_indices))
len(wild_indices.intersection(improved_indices))

"""

len(wild_indices - improved_indices)
Out[1191]: 35816

len(wild_indices - improved_indices - native_indices)
Out[1192]: 31407

len(native_indices - wild_indices - improved_indices)
Out[1193]: 14323

len(improved_indices - native_indices -wild_indices)
Out[1194]: 1539
"""


capra_281_plink_plt = capra_281_plink.query('type!="BND"')

capra_281_plink_plt['diff'] = capra_281_plink_plt['diff'].abs()
capra_281_plink_plt = capra_281_plink_plt.query('50<diff<=1000').sort_values(by=['diff'])


#如果diff<0,则将起点和终点交换 additional file 2 S16
#capra_281_noTRA_st.loc[capra_281_noTRA_st['diff'] < 0, ['POS', 'end']] = capra_281_noTRA_st.loc[capra_281_noTRA_st['diff'] < 0, ['end', 'POS']].values

"____________________________________________________________________________"
"Figure2c绘图_______________________________________________________________________"
#读取0/1格式的绵羊和山羊数据
sheep_532_plink_0110 = pd.read_csv('F:\svdata\sheep_532_plink_0110.vcf', sep='\t')
sheep_532_plink_0110['animal'] = "sheep"
sheep_532_plink_0110 = sheep_532_plink_0110.drop('END_values',axis=1)
capra_281_plink = pd.read_csv(r'F:\svdata\SV_count_result.csv', sep=',')
#打上标签后合并长度数据
"""
x = sheep_anno.iloc[:,:5]
x.columns = ['loc','gene','chr','s','t']
sheep_532 = sheep_532_plink_0110[['#CHROM','POS','end','type']]
sheep_532.columns = ['chr','s','t','type']
x = x.query('loc!="intergenic"')
"""







capra_281_plink['animal'] = "goat"
sheep_diff = sheep_532_plink_0110[['type','diff','animal']]
goat_diff = capra_281_plink[['type','diff','animal']]
all_diff = pd.concat([sheep_diff, goat_diff])
#去NA、取绝对值、过滤到50-1000
all_diff = all_diff.dropna()
all_diff['diff'] = all_diff['diff'].abs()
#50-1000、1000-10000、>10000一共分3个梯度
diff1 = all_diff.query('50<diff<=1000')
diff2 = all_diff.query('1000<diff<=10000')
diff3 = all_diff.query('10000<diff<100000')
diff4 = all_diff.query('100000<diff')
import matplotlib.pyplot as plt

import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess
bin1 = np.arange(50, 1001, 10)
bin2 = np.arange(1000, 10001, 100)
bin3 = np.arange(10000, 100000, 1000)
bin4 = np.arange(100000, 1000000, 10000)
def diff_binplot(all_diff, bins, number):
    # 计算diff值的计数，并对计数取对数
    counts = all_diff['diff'].value_counts().sort_index()
    log_counts = np.log10(counts)
    
    #bins = np.arange(50, 1001, 10)
    all_diff['bin'] = pd.cut(all_diff['diff'], bins=bins)
    binned_data = all_diff.groupby(['bin', 'type', 'animal'])['diff'].agg(['count']).reset_index()
    # 取对数
    binned_data['log_count'] = np.log10(binned_data['count'])
    # 创建图表
    plt.figure(figsize=(14, 8))
    
    # 分别绘制sheep和goat的曲线
    for animal in ['sheep', 'goat']:
        animal_data = binned_data[binned_data['animal'].str.lower() == animal]
        for sv_type in animal_data['type'].unique():
            type_data = animal_data[animal_data['type'] == sv_type]
            linestyle = '-' if animal == 'sheep' else '--'  # sheep为实线，goat为虚线
            plt.plot(type_data['bin'].apply(lambda x: x.mid), type_data['log_count'], linestyle=linestyle, label=f'{animal.capitalize()}-{sv_type}')
    
    # 添加图例
    plt.legend()
    
    # 添加标题和坐标轴标签
    plt.title('Log10 Counts of SV Types by Animal Type')
    plt.xlabel('SV Size (10 BP bins)')
    plt.ylabel('Log10(Count)')
    
    # 调整图表布局并展示
    plt.tight_layout()
    plt.savefig(r"F:\svdata\svplot\Figure2cS{}.pdf".format(number))
    plt.show()

"____________________________________________________________________________"
"Figure2d绘图_______________________________________________________________________"
def calculate_maf(row):
    # 计算每个等位基因的计数
    allele_counts = row.str.cat(sep='').count('0'), row.str.cat(sep='').count('1')

    # 计算总的有效等位基因数（排除了'./.'）
    total_alleles = sum(allele_counts)

    # 计算次等位基因频率
    if total_alleles == 0:  # 防止除以零
        return None
    maf = min(allele_counts) / total_alleles
    return maf

# 应用函数计算每个位点的MAF
maf_sheep = sheep_532_plink_0110.iloc[:, 10:-3].apply(calculate_maf, axis=1)
maf_goat = capra_281_plink.iloc[:, 10:-3].apply(calculate_maf, axis=1)
sheep_532_plink_0110['MAF'] = maf_sheep
capra_281_plink['MAF'] = maf_goat

all_maf = pd.concat([sheep_532_plink_0110[['type','MAF','animal']], capra_281_plink[['type','MAF','animal']]])

# 计算maf值的计数，并对计数取对数
counts = all_maf['MAF'].value_counts().sort_index()
log_counts = np.log10(counts)

bins = np.arange(-0.000001, 0.51, 0.01)
all_maf['bin'] = pd.cut(all_maf['MAF'], bins=bins)
binned_data = all_maf.groupby(['bin', 'type', 'animal'])['MAF'].agg(['count']).reset_index()
# 取对数
binned_data['log_count'] = np.log10(binned_data['count'])
# 创建图表
plt.figure(figsize=(14, 8))

# 分别绘制sheep和goat的曲线
for animal in ['sheep', 'goat']:
    animal_data = binned_data[binned_data['animal'].str.lower() == animal]
    for sv_type in animal_data['type'].unique():
        type_data = animal_data[animal_data['type'] == sv_type]
        linestyle = '-' if animal == 'sheep' else '--'  # sheep为实线，goat为虚线
        plt.plot(type_data['bin'].apply(lambda x: x.mid), type_data['log_count'], linestyle=linestyle, label=f'{animal.capitalize()}-{sv_type}')

# 添加图例
plt.legend()

# 添加标题和坐标轴标签
plt.title('Log10 Counts of SV Types by Animal Type')
plt.xlabel('SV Size (10 BP bins)')
plt.ylabel('Log10(Count)')

# 调整图表布局并展示
plt.tight_layout()
plt.savefig(r"F:\svdata\svplot\Figure2d.pdf")
plt.show()
"____________________________________________________________________________"
"Figure2F绘图_______________________________________________________________________"
data = {
    'type': ['A', 'B', 'A', 'C', 'B', 'A', 'C', 'A', 'B', 'C'],
    'BQ2602': ['0/0', '0/1', '1/1', '0/0', '1/1', '0/0', '1/1', '0/1', '1/1', '0/0'],
    'BQ2620': ['0/1', '0/0', '0/1', '1/1', '0/0', '1/1', '0/0', '0/1', '1/1', '0/0'],
    # ... 其他列
}

# 选取需要分析的列
genotypes = capra_281_plink.iloc[:, 10:-3]

# 初始化一个空的DataFrame来存储结果
type_counts_per_individual = pd.DataFrame(index=capra_281_plink['type'].unique(),columns=genotypes.columns)
# 对于每个个体（即DataFrame的每一列），统计各种type出现的次数
for individual in genotypes.columns:
    # 筛选出符合条件的行索引
    valid_indices = genotypes[individual].isin(['0/1', '1/1'])
    
    # 根据这些索引，统计type的出现次数
    counts = capra_281_plink.loc[valid_indices, 'type'].value_counts()
    
    # 将统计结果添加到结果DataFrame中
    type_counts_per_individual[individual] = counts

type_counts_per_individual = type_counts_per_individual.T
name_data = pd.read_csv(r'F:\svdata\Fig2Edata.txt', sep='\t', header=None)
name_data.columns = ['breed', 'type', 'name', 'srr']
#type_counts_per_individual = type_counts_per_individual.query('type != "Unknown" and type != "Hybrid"')
# 展示结果
type_counts_per_individual


type_counts_per_individual['breed'] = name_data['breed'].to_list()
type_counts_per_individual['type'] = name_data['type'].to_list()


merged_data = type_counts_per_individual.groupby(['breed', 'type']).mean().reset_index()
merged_data['count'] = merged_data.DEL + merged_data.INV + merged_data.DUP + merged_data.TRA
merged_data = merged_data.sort_values(by=['count'], ascending=False)
# 设置颜色映射
colors = {'DEL': 'blue', 'DUP': 'purple', 'INV': 'green', 'TRA': 'orange'}

# 绘制堆叠图
fig, ax = plt.subplots(figsize=(10, 6))

# 用于存储每个品种的底部位置的变量
bottom = np.zeros(len(merged_data))

# 循环遍历每种类型的结构变异，将它们堆叠起来
for sv_type, color in colors.items():
    ax.bar(merged_data['breed'], merged_data[sv_type], bottom=bottom, label=sv_type, color=color)
    bottom += merged_data[sv_type]

# 添加图例
ax.legend()

# 设置x轴和y轴标签
ax.set_xlabel('Breed')
ax.set_ylabel('Counts')

# 设置标题
ax.set_title('Stacked Counts of Structural Variants by Breed')
plt.savefig(r"F:\svdata\svplot\Figure2e.pdf")
# 显示图表
plt.show()
"____________________________________________________________________________"
"Figure2e_1绘图_______________________________________________________________________"
type_colors = {
    'Wild goat': 'red',
    'Native goat': 'blue',
    'Improved goat': 'green'
}

new_data = merged_data.loc[merged_data.breed.drop_duplicates().index]
# Create a new figure
plt.figure(figsize=(15, 1))  # Width is larger to spread out the points, height is small as we only need one row

# Plot a row of points with colors based on the type
for i, t in enumerate(new_data['type']):
    plt.scatter(i, 0, color=type_colors[t], label=t if t not in plt.gca().get_legend_handles_labels()[1] else "")

# Customize the plot
plt.gca().axes.get_yaxis().set_visible(False)  # Hide the y-axis as it's not needed
plt.gca().axes.get_xaxis().set_visible(False)  # Optionally hide the x-axis for cleaner look
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Type')  # Place the legend to the right of the figure
plt.savefig(r"F:\svdata\svplot\Figure2e_1.pdf")
# Show the plot
plt.show()



#用bedtools取QTL交集


r'F:\svdata\capra_sample281.anno.variant_function'
####

intersection = pd.read_csv(r'F:\svdata\intersection_over_50_percent.bed',header=None,sep='\t')
sheep_merge_annovar = pd.read_csv(r'F:\svdata\sheep_sv.merged.vcf.annovar.input',header=None,sep='\t')
peak_sv = intersection[[3,4,5]]
peak_sv = peak_sv.drop_duplicates()
peak_sv.columns = ['chr','s','t']
sheep_merge_annovar.columns = ['chr','s','t','a','b','c','e','h']

diff = pd.merge(sheep_merge_annovar, peak_sv, on=['chr', 's', 't'], how='outer', indicator=True).query('_merge != "both"').drop('_merge', axis=1)


#peak_sv.to_csv(r'F:\svdata\peak_sv.bed',index=None,header=None,sep='\t')
pd.merge(sheep_merge_annovar, peak_sv, on=['chr', 's', 't']).to_csv(r'F:\svdata\peak_sv.bed',index=None,header=None,sep='\t')
diff.to_csv(r'F:\svdata\nopeak_sv.bed',index=None,header=None,sep='\t')
#linux系统上annovar注释结果
nopeak_anno = pd.read_csv(r'F:\svdata\sheep_sv.nopeak.anno.variant_function',header=None,sep='\t')
peak_anno = pd.read_csv(r'F:\svdata\sheep_sv.peak.anno.variant_function',header=None,sep='\t')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 定义数据
nopeak_data = pd.Series({
    "downstream": 607,
    "exonic": 983,
    "intronic": 6948,
    "upstream": 563,
    "3'UTR": 726,
    "5'UTR": 216,
    "intergenic": 16040
})

peak_data = pd.Series({
    "downstream": 114,
    "exonic": 2099,
    "intronic": 975,
    "upstream": 137,
    "3'UTR": 26,
    "5'UTR": 153,
    "intergenic": 2567
})

# 计算比例
nopeak_ratios = nopeak_data / nopeak_data.sum()
peak_ratios = peak_data / peak_data.sum()
all_sv_ratios = (nopeak_data + peak_data) / (nopeak_data + peak_data).sum()
"____________________________________________________________________________"
"Figure7堆积图绘图_______________________________________________________________________"
# 计算堆积数据的起始位置
def get_starts(ratios):
    starts = [0]
    for i in range(len(ratios) - 1):
        starts.append(starts[i] + ratios[i])
    return starts

nopeak_starts = get_starts(nopeak_ratios)
peak_starts = get_starts(peak_ratios)
all_sv_starts = get_starts(all_sv_ratios)

# 绘制堆积图
fig, ax = plt.subplots(figsize=(10, 3))

categories = nopeak_data.index
bar_width = 0.4

ax.barh('nonPeak-SVs', nopeak_ratios, left=nopeak_starts, color=['grey', 'orange', 'green', 'yellow', 'red', 'purple', 'blue'], edgecolor='white', height=bar_width)
ax.barh('peak-SVs', peak_ratios, left=peak_starts, color=['grey', 'orange', 'green', 'yellow', 'red', 'purple', 'blue'], edgecolor='white', height=bar_width)
ax.barh('all-SVs', all_sv_ratios, left=all_sv_starts, color=['grey', 'orange', 'green', 'yellow', 'red', 'purple', 'blue'], edgecolor='white', height=bar_width)

ax.set_xlabel('Percentage')
ax.set_title('SVs Overlap with Peaks')
ax.set_yticks([])
ax.set_xlim(0, 1)

# 添加分类标签
for i, cat in enumerate(categories):
    ax.text(-0.05, 1, cat, ha='right', va='center', fontsize=8, color='black', transform=ax.transAxes)

# 添加百分比标签
for i, (nopeak_ratio, peak_ratio, all_sv_ratio) in enumerate(zip(nopeak_ratios, peak_ratios, all_sv_ratios)):
    ax.text(nopeak_starts[i] + nopeak_ratio / 2, 1, f"{nopeak_ratio:.2%}", ha='center', va='center', fontsize=8, color='black', transform=ax.transAxes)
    ax.text(peak_starts[i] + peak_ratio / 2, 0.5, f"{peak_ratio:.2%}", ha='center', va='center', fontsize=8, color='black', transform=ax.transAxes)
    ax.text(all_sv_starts[i] + all_sv_ratio / 2, 0, f"{all_sv_ratio:.2%}", ha='center', va='center', fontsize=8, color='black', transform=ax.transAxes)

plt.tight_layout()
plt.show()

######################common gene
gene_list = 'ENSOARG00000006800,SLAMF7,SLAMF1,KITLG,HMGI-C,HERC5,HERC6,SLC34A2,POP1,NBEA,CRYL1,RNF213,U1,HBE1,TRIP13,SLC12aA7,SUPT3H,EXOC2,DUSP22,HBM,LUC7L,MTMR7'

'''
my_list = 'AGBL1,ANO3,ATP8A1,ATRNL1,BMPR2,CCBE1,CCSER1,CD226,CDH2,CNTNAP2,CPQ,CSMD1,CTNNA3,GABRB3,GALNTL6,GPC6,GRID2,GRM8,GTDC1,HCN1,KCNB2,KCNIP4,KLHL1,MAGI2,MECOM,MGAT5,PARD3B,PCDH15,PIGK,PLCB1,PRKAR1B,PRKG1,PTPRN2,RASGRP1,RFX3,RNASEH2B,SH3D19,SH3GL2,SHC3,SVEP1,THSD7A,THSD7B,TMEM117,VWA8'
'''

all_peak_anno = pd.read_csv('F:\\svdata\\peak.txt',sep='\t')

BMPR2 = all_peak_anno.query('geneId=="BMPR2"')
BMPR1B = all_peak_anno.query('geneId=="BMPR1B"')



#######################附图4
"____________________________________________________________________________"
"附图S4_______________________________________________________________________"
# 假设你的VCF文件已经加载到DataFrame中，列名包括'CHROM', 'POS'等基本信息
# 这里使用一个模拟的DataFrame来代表VCF文件的内容


# 假设的VCF数据
#第二圈
vcf_data = capra_281_plink.iloc[:,:2]
#vcf_data = sheep_532_plink_0110.iloc[:,:2]
# 计算每个染色体的最大长度（Mb）
max_pos_per_chrom = vcf_data.groupby('#CHROM')['POS'].max()
max_pos_per_mb = (max_pos_per_chrom / 1e6).apply(np.ceil).astype(int)

# 创建包含所有染色体、所有可能Mb区间的DataFrame
full_index = pd.MultiIndex.from_product([max_pos_per_mb.index, range(max_pos_per_mb.max() + 1)], names=['#CHROM', 'Mb_bin'])
full_df = pd.DataFrame(index=full_index).reset_index()

# 计算VCF中每个SNP所在的Mb区间
vcf_data['Mb_bin'] = (vcf_data['POS'] / 1e6).astype(int)

# 统计每个Mb区间的SNP数目
snp_counts = vcf_data.groupby(['#CHROM', 'Mb_bin']).size().reset_index(name='SNP_Count')

# 合并统计结果到全区间DataFrame中
result_df = pd.merge(full_df, snp_counts, on=['#CHROM', 'Mb_bin'], how='left')
result_df['SNP_Count'] = result_df['SNP_Count'].fillna(0).astype(int)

result_df['Start_pos'] = result_df['Mb_bin'] * 1e6  # 起始位置
result_df['End_pos'] = (result_df['Mb_bin'] + 1) * 1e6 - 1  # 结束位置，减1是因为区间是闭合的

# 选择需要的列并重新排列，去掉'Mb_bin'列
result_df = result_df[['#CHROM', 'Start_pos', 'End_pos', 'SNP_Count']]

result_df
capra_MB_result_df = result_df[~result_df['#CHROM'].astype(str).str.contains('NW')]
capra_MB_result_df.to_csv(r'F:\svdata\capra_MB_result_df.bed',index=None,sep='\t')

sheep_MB_result_df = result_df[~result_df['#CHROM'].astype(str).str.contains('NW')]
sheep_MB_result_df.to_csv(r'F:\svdata\sheep_MB_result_df.bed',index=None,sep='\t')

sheep_MB_result_df['#CHROM'] = 'a' + sheep_MB_result_df['#CHROM'].astype(str)
capra_MB_result_df['CHROM'] = 'a' + capra_MB_result_df['CHROM'].astype(str)

#第三圈
#F:\svdata\capra_87404_hotspot_FigS4

#第四圈
sheep_havegene_sv = sheep_anno[~sheep_anno[0].str.contains('intergenic')]
capra_anno = pd.read_csv(r'F:\svdata\capra_87404.sv.vcf.anno.variant_function',sep='\t',header=None)
sheep_havegene_sv = sheep_havegene_sv[[2,3,4,1]]

goat_havegene_sv = capra_anno[capra_anno[3] != "intergenic"]
goat_havegene_sv = goat_havegene_sv[~goat_havegene_sv[0].astype(str).str.contains('NW')]
goat_havegene_sv = goat_havegene_sv[[0,1,2,4]]

goat_havegene_sv.columns = ['Chromosome', 'Start', 'End', 'Gene_Info']
sheep_havegene_sv.columns = ['Chromosome', 'Start', 'End', 'Gene_Info']

def extract_genes(gene_info):
    genes = []
    for part in gene_info.split(','):
        if 'gene-' in part:
            gene = part.split('gene-')[-1].split('(')[0]
            genes.append(gene)
    return genes

# 应用函数提取基因，每个基因作为一行
goat_havegene_sv['Genes'] = goat_havegene_sv['Gene_Info'].apply(extract_genes)
exploded_goat_havegene_sv = goat_havegene_sv.explode('Genes')


# 确保没有基因的行被去除
exploded_goat_havegene_sv = exploded_goat_havegene_sv[exploded_goat_havegene_sv['Genes'].notna()]

# 计算每个区间的Mb，用于后续的分组
exploded_goat_havegene_sv['Start_Mb'] = (exploded_goat_havegene_sv['Start'] // 1e6) * 1e6
exploded_goat_havegene_sv['End_Mb'] = ((exploded_goat_havegene_sv['End'] // 1e6) + 1) * 1e6 - 1

# 按染色体和1Mb区间进行分组，计算每个区间的唯一基因数
gene_counts_per_mb = exploded_goat_havegene_sv.groupby(['Chromosome', 'Start_Mb', 'End_Mb'])['Genes'].nunique().reset_index(name='Unique_Gene_Count')

# 最终结果，包括染色体、区间起始和结束位置以及该区间的唯一基因数
gene_counts_per_mb = gene_counts_per_mb[['Chromosome', 'Start_Mb', 'End_Mb', 'Unique_Gene_Count']]
gene_counts_per_mb.rename(columns={'Start_Mb': 'Start', 'End_Mb': 'End'}, inplace=True)

gene_counts_per_mb.to_csv(r'F:\svdata\goat_MB_gene.bed',index=None,sep='\t')

#########################################
# 确保没有基因的行被去除
sheep_havegene_sv['Genes'] = goat_havegene_sv['Gene_Info'].apply(extract_genes)
exploded_sheep_havegene_sv = sheep_havegene_sv.explode('Genes')
exploded_sheep_havegene_sv = exploded_sheep_havegene_sv[exploded_sheep_havegene_sv['Genes'].notna()]

# 计算每个区间的Mb，用于后续的分组
exploded_sheep_havegene_sv['Start_Mb'] = (exploded_sheep_havegene_sv['Start'] // 1e6) * 1e6
exploded_sheep_havegene_sv['End_Mb'] = ((exploded_sheep_havegene_sv['End'] // 1e6) + 1) * 1e6 - 1



# 按染色体和1Mb区间进行分组，计算每个区间的唯一基因数
gene_counts_per_mb = exploded_sheep_havegene_sv.groupby(['Chromosome', 'Start_Mb', 'End_Mb'])['Genes'].nunique().reset_index(name='Unique_Gene_Count')

# 最终结果，包括染色体、区间起始和结束位置以及该区间的唯一基因数
gene_counts_per_mb = gene_counts_per_mb[['Chromosome', 'Start_Mb', 'End_Mb', 'Unique_Gene_Count']]
gene_counts_per_mb.rename(columns={'Start_Mb': 'Start', 'End_Mb': 'End'}, inplace=True)
gene_counts_per_mb.to_csv(r'F:\svdata\sheep_MB_gene.bed',index=None,sep='\t')

#共线性link圈


link = pd.read_csv(r'F:\svdata\link.xls',sep='\t',header=None)
link.columns = ['a','s1','t1','b','s2','t2']








link = link[~link.a.str.contains('NW')]
trans = pd.read_csv(r'F:\svdata\trans.txt',sep='\t',header=None)

replacement_dict = trans.set_index(0)[1].to_dict()
# 然后，我们将这个映射字典应用于 'link' DataFrame的相关列
link['a'] = link['a'].replace(replacement_dict)
link['b'] = link['b'].replace(replacement_dict)
link = link[~link.a.str.contains('NW')]
link[['b','s2','t2','a','s1','t1']].to_csv(r'F:\svdata\link',index=None,sep='\t',header=None)

####Fig S############
species = pd.read_excel('F:\svdata\species.xlsx',header=None)
#species = pd.read_csv('F:\svdata\sheep_species.txt',header=None,sep='\t')
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_unweighted
from itertools import combinations
result_df['type'] = capra_281_plink['type']
df2 = species
df2.columns = ['dq','id']
df1 = capra_281_plink
groups = defaultdict(list)
for index, row in df2.iterrows():
    groups[row['dq']].append(row[1])

individuals = df2['id'].tolist()

# 检查第一个DataFrame的columns中是否包含这些个体名
individual_columns = [col for col in df1.columns if col in individuals]

# 提取包含这些个体名的列
result_df = df1[individual_columns]

# 显示结果
print(result_df)
df1 = result_df

# 合并df1和df2以获取地理信息
# 将 df1 索引转换为列，用作 SNP 位置或编号
df1_reset = df1.reset_index().rename(columns={'index': 'SNP_Position'})

# 将 df1 转换为长格式，并包含 SNP 位置信息
df1_melted_corrected = df1_reset.melt(id_vars="SNP_Position", var_name="id", value_name="Genotype")

# 使用 df2 的地理信息合并
df_merged_corrected = pd.merge(df1_melted_corrected, df2, on="id")

# 筛选出非空（即 '0/1' 或 '1/1'）的 SNP 记录
df_valid_snps_corrected = df_merged_corrected[df_merged_corrected["Genotype"].isin(["0/1", "1/1"])]

# 按地区和 SNP 位点分组，然后计算每个地区的每个 SNP 位点上的非空记录数量
valid_snp_counts_corrected = df_valid_snps_corrected.groupby(["dq", "SNP_Position"]).size().reset_index(name="count")

# 对于每个地区，计算有效 SNP 的总数（即在至少一个个体中非空的 SNP 位点数量）
snp_counts_by_region_corrected = valid_snp_counts_corrected.groupby("dq")["SNP_Position"].nunique().reset_index(name="valid_snp_count")

snp_counts_by_region_corrected

valid_snp_counts_corrected.iloc[:,:2].query('dq=="Wild goat"').SNP_Position.to_csv(r'F:\svdata\Wild goat.txt',index=None,header=None,sep='\t')
valid_snp_counts_corrected.iloc[:,:2].query('dq=="Europe"').SNP_Position.to_csv(r'F:\svdata\Europe.txt',index=None,header=None,sep='\t')
valid_snp_counts_corrected.iloc[:,:2].query('dq=="Africa"').SNP_Position.to_csv(r'F:\svdata\Africa.txt',index=None,header=None,sep='\t')
valid_snp_counts_corrected.iloc[:,:2].query('dq=="Asia"').SNP_Position.to_csv(r'F:\svdata\Asia.txt',index=None,header=None,sep='\t')

##################################表格1统计
capra_281_plink = pd.read_csv(r'F:\svdata\SV_count_result.csv',sep=',')
df = capra_281_plink.iloc[:,10:-1]
df.columns = df.columns.map(lambda x: x.replace('.1', '').replace('.2', ''))
goat_breed = pd.read_csv(r'F:\svdata\goat_breed.txt',sep='\t',header=None)
#sheep_breed = pd.read_csv(r'F:\svdata\sheep_breed.txt',sep='\t',header=None)
sample_to_breed = goat_breed.set_index(1)[0].to_dict()
df_mapped = df.rename(columns=sample_to_breed)
# Function to check if a variant exists in a group of samples at a given position
# Initializing a DataFrame to store the mutation counts for each breed
breed_mutation_counts_corrected = pd.DataFrame(columns=["Breed", "DEL", "DUP", "INV", "TRA", "INS"])


df_mapped_T = df_mapped.T
df_mapped_T = df_mapped_T.reset_index()

def count_SV(df_mapped_T):
    new_row = pd.Series(0, index=df_mapped_T.columns)
    
    # 由于第一列是描述性的列（比如可能是样本名），我们将其保留为空或其他特定值
    new_row[0] = None  # 或者你可以赋予它一个特定的标识，如 'Contains 0/1 or 1/1'
    
    # 遍历DataFrame的列（除了第一列）
    for col in df_mapped_T.columns[1:]:
        # 检查每列是否包含 "0/1" 或 "1/1"
        if df_mapped_T[col].isin(['0/1', '1/1']).any():
            # 如果包含，设置对应的值为1
            new_row[col] = 1
    
    # 将新的行添加到DataFrame的末尾
    df_mapped_T = df_mapped_T.append(new_row, ignore_index=True)
    return df_mapped_T


domestication = df_mapped_T.query('index=="domestication"')
Alpine_Ibex = df_mapped_T.query('index=="Alpine Ibex"')
Bezoar = df_mapped_T.query('index=="Bezoar"')
Siberian_Ibex = df_mapped_T.query('index=="Siberian Ibex"')
Iberian_Ibex = df_mapped_T.query('index=="Iberian Ibex"')
Nubian_Ibex = df_mapped_T.query('index=="Nubian Ibex"')
Markhor = df_mapped_T.query('index=="Markhor"')
'''
sheep
domestication = df_mapped_T.query('index=="domestication"')
Ovis_musimon = df_mapped_T.query('index=="Ovis musimon"')
Ovis_orientalis = df_mapped_T.query('index=="Ovis orientalis"')
Ovis_vignei = df_mapped_T.query('index=="Ovis vignei"')
Ovis_ammon = df_mapped_T.query('index=="Ovis ammon"')
Ovis_nivicola = df_mapped_T.query('index=="Ovis nivicola"')
Ovis_dalli = df_mapped_T.query('index=="Ovis dalli"')
Ovis_canadensis = df_mapped_T.query('index=="Ovis canadensis"')
'''



name_list = [domestication,Alpine_Ibex,Bezoar,Siberian_Ibex,Iberian_Ibex,Nubian_Ibex,Markhor]
#name_list = [domestication,Ovis_musimon,Ovis_orientalis,Ovis_vignei,Ovis_ammon,Ovis_nivicola,Ovis_dalli,Ovis_canadensis]
domestication = count_SV(domestication)
Alpine_Ibex = count_SV(Alpine_Ibex)
Bezoar = count_SV(Bezoar)
Siberian_Ibex = count_SV(Siberian_Ibex)
Iberian_Ibex = count_SV(Iberian_Ibex)
Nubian_Ibex = count_SV(Nubian_Ibex)
Markhor = count_SV(Markhor)
domestication.iloc[-1,0] = 'tag'
Alpine_Ibex.iloc[-1,0] = 'tag'
Bezoar.iloc[-1,0] = 'tag'
Siberian_Ibex.iloc[-1,0] = 'tag'
Iberian_Ibex.iloc[-1,0] = 'tag'
Nubian_Ibex.iloc[-1,0] = 'tag'
Markhor.iloc[-1,0] = 'tag'

'''
sheep
domestication = count_SV(domestication)
Ovis_musimon = count_SV(Ovis_musimon)
Ovis_orientalis = count_SV(Ovis_orientalis)
Ovis_vignei = count_SV(Ovis_vignei)
Ovis_ammon = count_SV(Ovis_ammon)
Ovis_nivicola = count_SV(Ovis_nivicola)
Ovis_dalli = count_SV(Ovis_dalli)
Ovis_canadensis = count_SV(Ovis_canadensis)

domestication.iloc[-1,0] = 'tag'
Ovis_musimon.iloc[-1,0] = 'tag'
Ovis_orientalis.iloc[-1,0] = 'tag'
Ovis_vignei.iloc[-1,0] = 'tag'
Ovis_ammon.iloc[-1,0] = 'tag'
Ovis_nivicola.iloc[-1,0] = 'tag'
Ovis_dalli.iloc[-1,0] = 'tag'
Ovis_canadensis.iloc[-1,0] = 'tag'
'''



def count_2(df, type1):
    df = df.T
    df.columns = df.iloc[0,:]
    df = df.iloc[1:,:]
    df['type'] = type1
    df = df.query('tag==1')
    return df.type.value_counts()
# Iterating over each breed in the goat_breed DataFrame
for breed in goat_breed[0].unique():
    # Filter the columns for the current breed
    breed_columns = [col for col in df_mapped.columns if col == breed]
    
    # If the breed is not found among the columns, continue to the next iteration
    if not breed_columns:
        continue
    
    # Filter the DataFrame for the current breed and check for variants
    breed_df = df_mapped[breed_columns].apply(lambda row: row.str.contains('1/1|0/1').any(), axis=1)
    
    # Combine the filtered DataFrame with the 'type' column
    breed_df_combined = pd.concat([breed_df, df['type']], axis=1, keys=['Exists', 'Type'])
    
    # Filter rows where a variant exists and count the occurrences of each type
    counts = breed_df_combined[breed_df_combined['Exists']]['Type'].value_counts().reindex(['DEL', 'DUP', 'INV', 'TRA'], fill_value=0)
    
    # Append the counts to the DataFrame
    breed_mutation_counts_corrected = breed_mutation_counts_corrected.append({
        "Breed": breed,
        "DEL": counts.get('DEL', 0),
        "DUP": counts.get('DUP', 0),
        "INV": counts.get('INV', 0),
        "TRA": counts.get('TRA', 0),
        "INV": counts.get('TRA', 0)
    }, ignore_index=True)

breed_mutation_counts_corrected







######################注释
capra_anno_input = pd.read_csv(r'F:\svdata\capra_all_sv_type.vcf.annovar.input',header=None,sep='\t')
capra_anno_input['type'] = capra_281_plink['type']
capra_anno_input['END2'] = capra['END_values']
capra_anno_input['CHR2'] = capra['chr2']
capra_anno_input.to_csv(r'F:\svdata\capra_all_sv_type_0518.vcf.annovar.input',index=None,header=None,sep='\t')


sheep_anno_input = pd.read_csv(r'F:\svdata\sheep-532-have-sv.annovar.input',header=None,sep='\t')
sheep_anno_input['type'] = sheep['type']
sheep_anno_input['END2'] = sheep['end']
sheep_anno_input['CHR2'] = sheep['chr2']
sheep_anno_input.to_csv(r'F:\svdata\sheep_532.sv.vcf.annovar.typeinput',index=None,header=None,sep='\t')

'''
INFO = sheep[sheep['INFO'].str.contains('=TRA')][['#CHROM','POS','INFO']]
#提取TRA信息
# Create a DataFrame
df = INFO

# Function to extract data between given start and end markers
def extract_data(s, start, end):
    pattern = re.compile(r'{}(.*?){}'.format(re.escape(start), re.escape(end)))
    matches = pattern.findall(s)
    return matches[0] if matches else None

# Extract 'CHR2' values
df['CHR2'] = df['INFO'].apply(lambda x: extract_data(x, 'CHR2=', ';'))

# Extract 'END' values, excluding the 'END=' and the trailing ';'
df['END2'] = df['INFO'].apply(lambda x: extract_data(x, 'END=', ';'))
sheep['chr2'] = sheep['#CHROM']
sheep.loc[sheep['type'] == "TRA", 'chr2'] = df['CHR2'].values


'''



###########################################Pi计算
BMPR1B_sheep_low_pi = pd.read_csv(r'F:\svdata\BMPR1B_sheep_low.pi.sites.pi', sep='\t')
BMPR1B_sheep_high_pi = pd.read_csv(r'F:\svdata\BMPR1B_sheep_high.pi.sites.pi', sep='\t')
BMPR1B_goat_low_pi = pd.read_csv(r'F:\svdata\BMPR1B_goat_low.pi.sites.pi', sep='\t')
BMPR1B_goat_high_pi = pd.read_csv(r'F:\svdata\BMPR1B_goat_high.pi.sites.pi', sep='\t')
mfl_pi = pd.read_csv(r'F:\svdata\mfl.pi.sites.pi', sep='\t')
bez_pi = pd.read_csv(r'F:\svdata\bez.pi.sites.pi', sep='\t')

BMPR1B_sheep_pi = pd.merge(BMPR1B_sheep_low_pi, BMPR1B_sheep_high_pi,on = ['CHROM','POS'])
BMPR1B_goat_pi = pd.merge(BMPR1B_goat_low_pi, BMPR1B_goat_high_pi,on = ['CHROM','POS'])

BMPR1B_sheep_mfl = pd.merge(BMPR1B_sheep_pi, mfl_pi, on=['CHROM','POS'])
BMPR1B_goat_bez = pd.merge(BMPR1B_goat_pi, bez_pi, on=['CHROM','POS'])

#提取出BMPR1B区域的所有位点的对应PI值 #34,034,777 30,214,561
BMPR1B_sheep_mfl.query('CHROM==6 and 32030000<POS<36520000').to_csv(r'F:\svdata\BMPR1B_sheep_pi.txt',sep='\t',index=None)
BMPR1B_goat_bez.query('CHROM=="6" and 28130000<POS<32320000').to_csv(r'F:\svdata\BMPR1B_goat_pi.txt',sep='\t',index=None)
#BMPR1B_goat_bez['pls'] = BMPR1B_goat_bez['PI_x']+BMPR1B_goat_bez['PI_y']+BMPR1B_goat_bez['PI']
#BMPR1B_goat_bez_nozero = BMPR1B_goat_bez.query('pls!=0')
#BMPR1B_goat_bez_nozero = BMPR1B_goat_bez_nozero.reset_index()
#BMPR1B_goat_bez_nozero.iloc[5176-25:5176+25].to_csv(r'F:\svdata\BMPR1B_goat_pi.txt',sep='\t',index=None)

native_goat_pi = pd.read_csv(r'F:\svdata\native_goat.pi.sites.pi', sep='\t')
native_sheep_pi = pd.read_csv(r'F:\svdata\native_sheep.pi.sites.pi', sep='\t')

native_goat_merge = pd.merge(native_goat_pi, bez_pi,on = ['CHROM','POS'])
native_sheep_merge = pd.merge(native_sheep_pi, mfl_pi,on = ['CHROM','POS'])
native_sheep_merge['pls'] = native_sheep_merge['PI_x']+native_sheep_merge['PI_y']
native_goat_merge['pls'] = native_goat_merge['PI_x']+native_goat_merge['PI_y']
native_sheep_merge = native_sheep_merge.query('pls!=0')
native_goat_merge = native_goat_merge.query('pls!=0')
native_sheep_merge = native_sheep_merge.reset_index() 
native_goat_merge = native_goat_merge.reset_index() 
#native_sheep_merge.query('POS==218822797') index 14109
#native_goat_merge.query('POS==44898207') index 2278
native_sheep_merge.iloc[14109-25:14109+25].to_csv(r'F:\svdata\BMPR1B_sheep_native_pi.txt',sep='\t',index=None)
native_goat_merge.iloc[2278-25:2278+25].to_csv(r'F:\svdata\BMPR1B_goat_native_pi.txt',sep='\t',index=None)


import pandas as pd

# 读取文件
sheep_hotspots_genes = pd.read_csv(r'F:\svdata\sheep_hotspots_genes.bed', sep='\t', header=None)
capra_hotspots_genes = pd.read_csv(r'F:\svdata\capra_hotspots_genes.bed', sep='\t', header=None)

# 提取基因ID列
sheep_hotspots_genes['gene_id'] = sheep_hotspots_genes[6].str.extract(r'ID=gene-([^;]+)')
capra_hotspots_genes['gene_id'] = capra_hotspots_genes[6].str.extract(r'ID=gene-([^;]+)')

# 找到交集基因
intersection_genes = pd.merge(sheep_hotspots_genes, capra_hotspots_genes, on='gene_id', suffixes=('_sheep', '_capra'))

# 创建一个字典来存储结果
hotspot_dict = {}
for index, row in intersection_genes.iterrows():
    chrom, start, end, gene_id = row[0], row[1], row[2], row['gene_id']
    key = (chrom, start, end)
    if key not in hotspot_dict:
        hotspot_dict[key] = []
    hotspot_dict[key].append(gene_id)

# 将结果转换为DataFrame
result = pd.DataFrame([
    {'chrom': chrom, 'start': start, 'end': end, 'genes': ','.join(genes)}
    for (chrom, start, end), genes in hotspot_dict.items()
])





