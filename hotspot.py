# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:50:33 2024

@author: h
"""

import pandas as pd

# 假设 VCF 文件名为 'example.vcf'
vcf = pd.read_csv( r'F:\svdata\example.txt', sep='\t', header=5)
#sheep.drop(sheep_532_plink_0110.query('type=="TRA"').index)
vcf = sheep.iloc[:,:8]
vcf_point = sheep
vcf = vcf_point
vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
# 定义一个函数来解析 INFO 字段并提取 SVTYPE, END, CHR2, END2
def parse_info(info_str):
    info_dict = dict(item.split('=') for item in info_str.split(';') if '=' in item)
    return info_dict.get('SVTYPE'), int(info_dict.get('END', -1)), info_dict.get('CHR2'), int(info_dict.get('END', -1))

# 应用这个函数来创建新的列
vcf[['SVTYPE', 'END', 'CHR2', 'END']] = vcf['INFO'].apply(lambda x: pd.Series(parse_info(x)))

# 创建一个空的 DataFrame 来存储断点
breakpoints = pd.DataFrame(columns=['CHROM', 'POS'])

# 对于非 Translocation 类型的 SV
non_tra = vcf[vcf['SVTYPE'] != 'TRA']
breakpoints = breakpoints.append(non_tra[['CHROM', 'POS']], ignore_index=True)
breakpoints = breakpoints.append(non_tra[['CHROM', 'END']].rename(columns={'END': 'POS'}), ignore_index=True)

# 对于 Translocation 类型的 SV
tra = vcf[vcf['SVTYPE'] == 'TRA']
breakpoints = breakpoints.append(tra[['CHROM', 'POS']], ignore_index=True)
breakpoints = breakpoints.append(tra[['CHR2', 'END']].rename(columns={'CHR2': 'CHROM', 'END': 'POS'}), ignore_index=True)

# 输出断点
print(breakpoints)

sheep_brekpoints.to_csv(r'F:\svdata\sheep_brekpoints.txt',sep='\t',header=None)
goat_breakpoints.to_csv(r'F:\svdata\goat_brekpoints.txt',sep='\t',header=None)


window = pd.read_csv(r'F:\svdata\sheep_output_windows.bed',header=None,sep='\t')
#window = pd.read_csv(r'F:\svdata\goat_output_windows.bed',header=None,sep='\t')
window.columns = ['CHROM', 'Start', 'End']
# 初始化用于计数的列表
counts = []
c = 0
# 对每个区间进行循环
for index, row in window.iterrows():
    c += 1
    if c%100 ==0:
        print(c)
    chrom, start, end = row['CHROM'], row['Start'], row['End']
    # 对断点进行筛选，计算在当前区间内的数量
    count = breakpoints[(breakpoints['CHROM'] == chrom) & (breakpoints['POS'] >= start) & (breakpoints['POS'] <= end)].shape[0]
    counts.append(count)

# 将计数结果添加到 window DataFrame
window['Count'] = counts
print(window)
window.sort_values(by='Count',ascending=False).iloc[:518]
hot_window = window.query('Count>85')

def merge_intervals(df):
    # 按染色体排序，以保证连续性
    df = df.sort_values(by=['CHROM', 'Start', 'End'])

    # 初始化合并后的区间列表
    merged_intervals = []

    # 初始化当前合并区间
    current_interval = df.iloc[0]

    for i in range(1, len(df)):
        next_interval = df.iloc[i]

        # 如果当前区间的结束位置与下一个区间的开始位置连续或重叠，并且染色体号相同
        if current_interval['End'] >= next_interval['Start'] - 1 and current_interval['CHROM'] == next_interval['CHROM']:
            # 合并区间并更新断点计数
            current_interval['End'] = max(current_interval['End'], next_interval['End'])
            current_interval['Count'] += next_interval['Count']
        else:
            # 添加当前区间到结果列表并开始一个新的区间
            merged_intervals.append(current_interval)
            current_interval = next_interval

    # 添加最后一个区间
    merged_intervals.append(current_interval)

    return pd.DataFrame(merged_intervals)

# 应用合并函数
merged_hot_windows = merge_intervals(hot_window)
sheep_window = merged_hot_windows
sheep_window.to_csv(r'F:\svdata\sheep_SV_location.bed',sep='\t',header=None)


counts = []
window_result = merged_hot_windows.iloc[:,:3]
# 对每个区间进行循环
for index, row in window_result.iterrows():
    chrom, start, end = row['CHROM'], row['Start'], row['End']
    # 对断点进行筛选，计算在当前区间内的数量
    count = breakpoints[(breakpoints['CHROM'] == chrom) & (breakpoints['POS'] >= start) & (breakpoints['POS'] <= end)].shape[0]
    counts.append(count)

# 将计数结果添加到 window DataFrame
window_result['Count'] = counts
print(window_result)

window_result.to_csv(r'F:\svdata\capra_hotspot_top10_count.bed', index=None, header=None, sep='\t')
window_result[['CHROM','Start','End']].to_csv(r'F:\svdata\goat_SV_location.bed',sep='\t',index=None,header=None)

sheep_window = pd.read_csv(r'F:\svdata\sheep_SV_location.bed',sep='\t',header=None)
goat_window = pd.read_csv(r'F:\svdata\capra_hotspot_top10_count.bed',sep='\t',header=None)

#################################
#telemere
sheep_telemere = pd.read_csv(r'F:\svdata\sheep_telemere_region.txt',sep='\t',header=None)













#########
sheep_hotspot_QTLoverlaps = pd.read_csv(r'F:\svdata\sheep_hotspot_QTLoverlaps.bed',sep='\t',header=None)
sheep_hotspot_QTLoverlaps


sheep_hotspot_sv_anno = pd.read_csv(r'F:\svdata\sheep_hotspot_sv.anno.variant_function',sep='\t',header=None)
sheep_hotspot_sv_anno = sheep_hotspot_sv_anno[[0,1,2,3,4]]

df = sheep_hotspot_sv_anno.iloc[:,:2]
df.columns = [0,1]
df = df[df[0]!='intergenic']

split_genes = df[1].str.split(',').explode()
cleaned_genes = split_genes.apply(clean_gene_name)
unique_genes_count = cleaned_genes.nunique()
df_exploded = pd.DataFrame({0: df[0].repeat(split_genes.groupby(split_genes.index).size()), 1: cleaned_genes})

df_exploded[1].drop_duplicates()
sheep_hotspot_sv_gene = df_exploded[1].drop_duplicates()


goat_hotspot_sv_anno = pd.read_csv(r'F:\svdata\goat_hotspot_sv.anno.variant_function',sep='\t',header=None)
goat_hotspot_sv_anno = goat_hotspot_sv_anno[[0,1,2,3,4]]
df = goat_hotspot_sv_anno.iloc[:,:2]
df.columns = [0,1]
df = df[df[0]!='intergenic']

split_genes = df[1].str.split(',').explode()
cleaned_genes = split_genes.apply(clean_gene_name)
unique_genes_count = cleaned_genes.nunique()
df_exploded = pd.DataFrame({0: df[0].repeat(split_genes.groupby(split_genes.index).size()), 1: cleaned_genes})

df_exploded[1].drop_duplicates()
goat_hotspot_sv_gene = df_exploded[1].drop_duplicates()

sheep_hotspot_sv_gene = pd.DataFrame(sheep_hotspot_sv_gene)
goat_hotspot_sv_gene = pd.DataFrame(goat_hotspot_sv_gene)

sheep_hotspot_sv_gene.columns = ['gene']

goat_hotspot_sv_gene.columns = ['gene']

hotspot_merge_gene = pd.merge(sheep_hotspot_sv_gene, goat_hotspot_sv_gene, on="gene")






sheep_hotspot_sv_anno.query('type!="intergenic"')

sheep_hotspot_region_gene = pd.read_csv(r'F:\svdata\sheep_hotspot_region_gene', sep='\t', header=None)

df = sheep_hotspot_region_gene.iloc[:]


# 提取基因名
df[8] = df[8].str.extract(r'gene-([^;]+?)(?:\(|$)')


# 合并具有相同前三列的行，并将对应的基因提取出来，成为第四列，以逗号分隔
df_grouped = df.groupby([0, 1, 2])[8].apply(lambda x: ','.join(x.unique())).reset_index()

# 打印结果
# 将基因名转换为集合
hotspot_genes_set = set(hotspot_merge_gene['gene'])

# 删除 gene_location 中基因名前缀 'gene-'
gene_location[3] = gene_location[3].str.replace('gene-', '')

# 过滤出 gene_location 中在 hotspot_merge_gene 中的基因
filtered_gene_location = gene_location[gene_location[3].isin(hotspot_genes_set)]

# 打印结果
print(filtered_gene_location)









# 将数据类型转为一致的整数类型
sheep_window = sheep_window.astype({'chr': int, 's': int, 't': int})
filtered_gene_location = filtered_gene_location.astype({'chr': int, 's': int, 't': int})

# 初始化结果列
sheep_window['genes'] = ""

# 创建一个字典来跟踪每个基因是否匹配
gene_matched = {gene: False for gene in filtered_gene_location['gene']}

# 遍历sheep_window和filtered_gene_location并找出交集
for i, row in sheep_window.iterrows():
    chr_window = row['chr']
    start_window = row['s']
    end_window = row['t']
    
    overlapping_genes = filtered_gene_location[
        (filtered_gene_location['chr'] == chr_window) &
        (filtered_gene_location['s'] < end_window) &
        (filtered_gene_location['t'] > start_window)
    ]['gene'].tolist()
    
    # 更新匹配状态
    for gene in overlapping_genes:
        gene_matched[gene] = True
    
    # 打印调试信息
    print(f"Window: chr={chr_window}, start={start_window}, end={end_window}, overlapping_genes={overlapping_genes}")

    sheep_window.at[i, 'genes'] = ','.join(overlapping_genes)

# 检查未匹配的基因
unmatched_genes = [gene for gene, matched in gene_matched.items() if not matched]
print(f"Unmatched genes: {unmatched_genes}")








