import numpy as np
import pandas as pd
import re
capra = pd.read_csv(r'F:\svdata\capra_sample281_flt_nohead.vcf',sep='\t')
sheep = pd.read_csv(r'F:\download\jump_max\sheep-532-have-sv.vcf',sep='\t', header=62)
location = pd.read_csv(r'F:\svdata\location.txt',sep=' ',header=None)
condition = (capra.iloc[:,9:].apply(lambda x: x.str.startswith("./.")).all(axis=1))

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


filt_capra = capra[~condition]
filt_capra

#提取breakpoint中的END数据
pattern = r'END=(\d+);'
# 使用 Pandas 的 str.extract 方法提取符合条件的字符串
filt_capra['END_values'] = filt_capra['INFO'].str.extract(pattern)
#
breakpoint_goat = filt_capra[['#CHROM','POS','END_values']]


capra_281_plink = pd.read_csv(r'F:\svdata\test.vcf',header=560,sep='\t')
capra_281_input = pd.read_csv(r'F:\svdata\capra_sample281.annovar.input',header=None,sep='\t')

#capra_281_plink.insert(2,'end',capra_281_input[2])
capra_281_plink.insert(2,'end',breakpoint_goat.END_values)

SV_type_f = capra_281_plink.iloc[:,3:6]

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



#capra_281_plink.insert(3,'len',capra_281_plink['end']-capra_281_plink['POS'])
#开始进行数目统计
#求长度diff
capra_281_plink_data = capra_281_plink.iloc[:,-282:-1]
#capra_281_plink.iloc[:,10:-3]

#重新将分配给BND的INV找出来
find_INV = capra[capra.applymap(lambda x: pd.notna(x) and ('INV:' in str(x)))]
INV = find_INV.dropna(how='all')
capra_281_plink.iloc[INV.index,-1] = "INV"

#将无法计算长度的TRA长度值赋予NA
capra_281_plink['diff'] = pd.to_numeric(capra_281_plink['end'], errors='coerce', downcast='float') - pd.to_numeric(capra_281_plink['POS'], errors='coerce', downcast='float')
#相减后要+1才是实际值
capra_281_plink['diff'] = capra_281_plink['diff'] + 1
capra_281_plink.iloc[capra_281_plink.query('type=="BND"').index,-1]=None
capra_281_plink_data['type'] = capra_281_plink['type']
capra_281_plink_data = capra_281_plink_data.replace('INV_1','INV')
capra_281_plink = capra_281_plink.replace('INV_1','INV')

#capra_281_plink.to_csv(r'F:\svdata\sheep_532_plink_0110.vcf',sep='\t',index=None)


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

#capra.iloc[capra_281_plink.index,:-1].to_csv(r'F:\svdata\capra_svout_87404_nohead.vcf',index=None,sep='\t')
capra_281_plink = capra_281_plink.reset_index(drop=True)

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
    all_SVtype_data.append([merged_type[[col]].columns[0], non_na_type_counts.DEL, non_na_type_counts.BND, non_na_type_counts.DUP, non_na_type_counts.INV, DEL_len, DUP_len, INV_len])
    
    #print(f"Counts for {col}:")
    #print(non_na_type_counts)
    #print()

SVtype_data = pd.DataFrame(all_SVtype_data)
SVtype_data.columns = ['name', 'DEL', 'BND', 'DUP', 'INV', 'DEL_len', 'DUP_len', 'INV_len']
#将vcf文件的sample中.1，.2这种删掉
SVtype_data.iloc[SVtype_data[SVtype_data.name.str.contains('\.')].index,0] = [i[:-2] for i in SVtype_data[SVtype_data.name.str.contains('\.')].name]

SVtype_data['all_number'] = SVtype_data.DEL + SVtype_data.BND + SVtype_data.DUP + SVtype_data.INV
SVtype_data.columns = ['name', 'DEL', 'BND', 'DUP', 'INV', 'DEL_len', 'DUP_len', 'INV_len', 'all_number']
SVtype_data['all_len'] = SVtype_data.DEL_len + SVtype_data.DUP_len + SVtype_data.INV_len

#根据名称信息将SRR转为对应的物种名
species = pd.read_csv(r'C:\Users\h\Desktop\SVdata\Species.txt',sep='\t')
species = species.iloc[:,[0,5,-1]]
species.columns = ['Species','code','name']
SV_result = pd.merge(species, SVtype_data, on='name')
SV_result.to_excel(r'F:\SV项目\GB submission\Table5.xlsx')
#最终结果
capra_281_plink.to_csv(r'F:\svdata\SV_count_result.csv',index=None)

#统计table S6
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
    j = len(capra_281_plink.query('type==@SV and 100000<diff<=50000f0'))
    k = len(capra_281_plink.query('type==@SV and 500000<diff<=1000000'))
    return [a,b,c,d,e,f,g,h,i,j,k]

capra_281_plink['diff'] = capra_281_plink['diff'].abs()
all_count_length = []
all_count_length.append(count_length('DEL'))
all_count_length.append(count_length('DUP'))
all_count_length.append(count_length('INV'))
pd.DataFrame(all_count_length)
#capra_281_plink[['#CHROM','POS','end','type']].query('type=="DEL"').to_csv(r'F:\svdata\capra_281_DEL.bed',index=None,header=None,sep='\t')
'''
Out[411]: 
      0      1      2     3     4     5     6    7   8   9   10
0  20054  16617  10636  8010  5850  2932  1002  462  77  72  15
1    764   3076   1126   652   527   184    92  111  24  12   1
2     80    224    351   739   676   271   106   86  11  29  11
'''
####Fig 
species = pd.read_excel('F:\svdata\species.xlsx',header=None)
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

##################################table1统计
capra_281_plink = pd.read_csv(r'F:\svdata\SV_count_result.csv',sep=',')
df = capra_281_plink.iloc[:,10:-2]
df.columns = df.columns.map(lambda x: x.replace('.1', '').replace('.2', ''))
goat_breed = pd.read_csv(r'F:\svdata\goat_breed.txt',sep='\t',header=None)
sample_to_breed = goat_breed.set_index(1)[0].to_dict()
df_mapped = df.rename(columns=sample_to_breed)
# Function to check if a variant exists in a group of samples at a given position
# Initializing a DataFrame to store the mutation counts for each breed
breed_mutation_counts_corrected = pd.DataFrame(columns=["Breed", "DEL", "DUP", "INV", "TRA"])


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

name_list = [domestication,Alpine_Ibex,Bezoar,Siberian_Ibex,Iberian_Ibex,Nubian_Ibex,Markhor]

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
        "TRA": counts.get('TRA', 0)
    }, ignore_index=True)

breed_mutation_counts_corrected





