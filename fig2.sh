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
"""
######################对sheep的SV进行分类
find_DEL = sheep[sheep.applymap(lambda x: pd.notna(x) and ('DEL:' in str(x)))].dropna(how="all")
find_INV = sheep[sheep.applymap(lambda x: pd.notna(x) and ('INV:' in str(x)))].dropna(how="all")
find_DUP = sheep[sheep.applymap(lambda x: pd.notna(x) and ('DUP:' in str(x)))].dropna(how="all")
find_TRA = sheep[sheep.applymap(lambda x: pd.notna(x) and ('TRA:' in str(x)))].dropna(how="all")
sheep_532_plink_0110.loc[find_DEL.index,"type"]='DEL'
sheep_532_plink_0110.loc[find_INV.index,"type"]='INV'
sheep_532_plink_0110.loc[find_DUP.index,"type"]='DUP'
sheep_532_plink_0110.loc[find_TRA.index,"type"]='TRA'
sheep_532_plink_0110.
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
def diff_binplot(all_diff, bins):
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
    plt.savefig(r"F:\svdata\svplot\Figure2cS1.pdf")
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
merged_data['count'] = merged_data.DEL + merged_data.INV + merged_data.DUP + merged_data.BND
merged_data = merged_data.sort_values(by=['count'], ascending=False)
# 设置颜色映射
colors = {'DEL': 'blue', 'DUP': 'purple', 'INV': 'green', 'BND': 'orange'}

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