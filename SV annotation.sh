"____________________________________________________________________________"
"注释结果与基因数目和分布统计_______________________________________________________________________"
#Capra
#capra_anno = pd.read_csv(r'F:\svdata\capra_87404.sv.vcf.anno.variant_function',sep='\t',header=None)
capra_anno = pd.read_csv(r'F:\svdata\capra_87404.type.anno.variant_function',sep='\t',header=None)
#capra_anno = capra_anno.iloc[:,3:]
capra_anno.columns = [0,1,2,3,4,5,6,7,8,9,10]
capra_anno_havegene = capra_anno[capra_anno[0] != "intergenic"]
"""
df_exploded[0].value_counts()
Out[196]: 
intronic          30718
exonic             1995
upstream           1994
downstream         1913
UTR3               1490
splicing           1089
ncRNA_intronic      854
UTR5                731
ncRNA_exonic        311
ncRNA_splicing       31
Name: 0, dtype: int64
"""
count_gene_occurrences = capra_anno_havegene[1].str.contains('gene-').sum()
#用正则表达式提取所有基因
pattern = 'gene-(.*?)(?=,|$|\()'
#提取基因
all_gene_location = capra_anno_havegene[1].str.extractall(pattern).astype(str)
all_gene_loc

"____________________________________________________________________________"
"计算有多少个不重复基因，基因被外显子注释、内含子注释、还是都有。_______________________________________________________________________"
#其中一共有12908个基因被SV注释到，其中1725个基因被外显子注释到，454个基因只被外显子注释到。
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
        if gene in common_gene['gene'].values:
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

# Create a new DataFrame from the expanded rows
expanded_df = pd.DataFrame(expanded_rows).reset_index(drop=True)

expanded_df

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

#Sheep
sheep_anno = pd.read_csv(r'F:\svdata\sheep_532.type.anno1.variant_function',sep='\t',header=None)
#sheep_anno = pd.read_csv(r'F:\svdata\sheep-532-have-sv.annovar.variant_function',sep='\t',header=None)
sheep_anno_havegene = sheep_anno[sheep_anno[0] != "intergenic"]
all_gene_location = sheep_anno_havegene[1].str.extractall(pattern).astype(str)
all_gene = all_gene_location.drop_duplicates()
result_str = ','.join(all_gene[0])
gene_df_sheep = pd.DataFrame(result_str.split(','))
gene_df_sheep = gene_df_sheep[~gene_df[0].str.startswith('LOC')]
#计算SV注释到的merge基因
sheep_goat_merge = pd.merge(gene_df_sheep, gene_df, on='gene')

#计算山羊和绵羊SV注释到基因的交集 5948
pd.merge(gene_df, gene_df_sheep)
#common_gene.sort_values(by=[0])