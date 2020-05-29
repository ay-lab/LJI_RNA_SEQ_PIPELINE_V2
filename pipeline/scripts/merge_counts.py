import json
import pandas as pd
from functools import reduce

def read_json(file = 'conf_RNA_Seq.json'):
    with open(file) as json_file:
        conf = json.load(json_file)
    return conf

# Load in params
dict_conf = read_json(snakemake.params[0])
sample_list = snakemake.params[1]

# get counts
df_anno = pd.read_csv(dict_conf['config']['annotation_file'],index_col = 0)   
count_tables = [pd.read_csv(f'2.Internal_files/bam_aligned/{sample}/{sample}_ReadsPerGene.out.tab',sep = '\t',skiprows = 4,header = None).set_index(0)[[1]] for sample in sorted(sample_list)]
df_counts = reduce(lambda left,right: pd.merge(left,right,left_index = True, right_index = True), count_tables)
df_counts.columns = sorted(sample_list)
filter_gene_list = df_anno[df_anno["gene_type"].str.contains('tRNA|rRNA')==False].index
df_counts_filter = df_counts.loc[filter_gene_list]
counts = df_counts_filter.values/df_anno['Length'].loc[df_counts_filter.index].values[:,None]
df_TPM = pd.DataFrame((1e6*counts)/(counts.sum(axis = 0)), index = filter_gene_list, columns=df_counts.columns)
df_TPM.to_csv('4.Output/counts/TPM_counts.csv')
df_counts.to_csv('4.Output/counts/raw_counts_unfiltered.csv')
df_counts_filter.to_csv('4.Output/counts/raw_counts.csv')