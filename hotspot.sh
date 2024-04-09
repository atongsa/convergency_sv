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