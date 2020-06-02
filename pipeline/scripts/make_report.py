# general libs
import json
from functools import reduce
import pandas as pd
pd.set_option('display.max_colwidth', -1) # make sure generated links have full string
import numpy as np

# modeling libs
from scipy import stats
from scipy.stats import spearmanr,norm
from scipy.optimize import least_squares

####################### Rules for QC  #######################

# Here is the defination of QC rule logic; for changing thresholds please modify in the configuration file. 

def RNA_QC(sample_tuple,dict_threshold,count_X_low,count_Y_low):
    
    minimal_counts = dict_threshold['minimal_counts']
    
    if minimal_counts == 'fixed':# Reseq based fixed total STAR count 
        if sample_tuple.final_STAR_counts > dict_threshold['final_STAR_counts']:
            A = True
        else:
            A = False 
    elif minimal_counts == 'perc': # Reseq based on % of gene recovery rate if the threshold is higher than the fixed value
        if sample_tuple.final_STAR_counts > max(count_X_low,dict_threshold['final_STAR_counts']):
            A = True
        else:
            A = False 
    if sample_tuple.uniquely_mapped_reads_perc > dict_threshold['uniquely_mapped_reads_perc']:
        B = True
    else:
        B = False
    if sample_tuple.exonic_perc > dict_threshold['exonic_perc']:
        C = True
    else:
        C = False
    if sample_tuple.too_short_reads_perc < dict_threshold['too_short_reads_perc']:
        D = True
    else:
        D = False
    if sample_tuple.t_rRNA_counts_perc < dict_threshold['t_rRNA_counts_perc']:
        E = True
    else:
        E = False
    # 
    if sample_tuple.Total_genes  < max(dict_threshold['Total_genes'][0],min(dict_threshold['Total_genes'][1],count_Y_low)):
        F = False
    else:
        F = True
    if (float(str(sample_tuple.bias_5to3_prim).replace(',','')) < dict_threshold['bias_5to3_prim'][1])&(float(str(sample_tuple.bias_5to3_prim).replace(',','')) > dict_threshold['bias_5to3_prim'][0]): # In rear cases the bias value is small; also not accepted
        G = True
    else:
        G = False
    if (sample_tuple.insert_median > dict_threshold['insert_median'][0]) or (sample_tuple.insert_median == 0):# and (sample_tuple.insert_median < dict_threshold['insert_median'][1]):
        H = True
    else:
        H = False
     
    if A and B and C and D and E and F and G and H:
        return '1.Good','None'
    elif B and C and D and E and F and G and H:
        return '2.Reseq','counts='+str(int(sample_tuple.final_STAR_counts))+'threshold='+str(int(max(count_X_low,dict_threshold['final_STAR_counts'])))
    else:
        return '3.Manual QC','/'.join(name for item,name in zip([A,B,C,D,E,F,G,H],['final_STAR_counts','uniquely_mapped_reads_perc','exonic_perc','too_short_reads_perc','t_rRNA_counts_perc','Total_genes','bias_5to3_prim','insert_median']) if item == False)


####################### basic functions for reading QC files #######################

def Hawk_smash(List):
    '''Flattern lists within a list '''
    return [item for sublist in List for item in sublist]

def platt(coef,X,Y):
    # Modified Platt equation to allow more shape turn at the tipping point
    return coef[0]*(1-np.exp(-(coef[1]*X)/coef[0]))-Y

def fit_platt(list_count_gene):
    # guess init aplha value and set point 1 to zero
    G_alpha = stats.linregress(list_count_gene[0][0:5],list_count_gene[1][0:5]).slope
    G_P = max(list_count_gene[1])
    coef_0 = np.array([G_P,G_alpha], dtype=float)
    # Start non-linear fitting
    res_lsq = least_squares(platt, coef_0, loss='linear',  args=(list_count_gene[0], list_count_gene[1]))
    return res_lsq.x

def read_fastp(sample):
    with open(f'2.Internal_files/fastp_output/{sample}_fastp.json') as json_file:
        data = json.load(json_file)
        try:
            if data['summary']['before_filtering']['read2_mean_length']>0:
                scalar = 2
        except KeyError:
            scalar = 1
        total_reads = int(data['summary']['before_filtering']['total_reads']/scalar)
        filtered_reads = int(data['summary']['after_filtering']['total_reads']/scalar)
        filtered_reads_perc = round(100*filtered_reads/total_reads,2)
        try:
            adaptor_trimm_perc = round(100*data['adapter_cutting']['adapter_trimmed_reads']/data['summary']['before_filtering']['total_reads'],2)
        except KeyError:
            adaptor_trimm_perc = 0
        dup_rate = round(100*data['duplication']['rate'],2)
        link = f'<a href="../3.QC_files/fastp_report/{sample}_fastp.html">seq_QC</a>'
    return [total_reads,filtered_reads,filtered_reads_perc,adaptor_trimm_perc,dup_rate,link]

def read_STAR(sample):
    with open(f'2.Internal_files/bam_aligned/{sample}/{sample}_Log.final.out') as STAR_file:
        for l in STAR_file.readlines():
            if 'Uniquely mapped reads number' in l:
                unique_reads = int(l.split('|')[-1][1:-1])
            if 'Uniquely mapped reads %' in l:
                unique_reads_perc = float(l.split('|')[-1][1:-2])    
            if 'Number of splices: Total' in l:
                spliced_reads = int(l.split('|')[-1][1:-1])
            if 'Number of splices: Annotated (sjdb)' in l:
                anno_spliced_reads = int(l.split('|')[-1][1:-1])
            if 'Number of reads unmapped: too short' in l:
                too_short_reads = int(l.split('|')[-1][1:-1])
            if '% of reads unmapped: too short' in l:
                too_short_reads_perc = float(l.split('|')[-1][1:-2])
    return [unique_reads,unique_reads_perc,spliced_reads,anno_spliced_reads,too_short_reads,too_short_reads_perc]

def read_Qualimap(sample):
    with open(f'3.QC_files/qualimap/{sample}/rnaseq_qc_results.txt') as Qualimap_file:
        for l in Qualimap_file.readlines():
            if 'exonic' in l:
                Exonic_perc = float(l.split('=')[-1].split('(')[1][:-3])
            if 'intronic' in l:
                intronic_perc = float(l.split('=')[-1].split('(')[1][:-3])
            if 'intergenic' in l:
                intergenic_perc = float(l.split('=')[-1].split('(')[1][:-3])
            if "5' bias" in l:
                bias_5_prim = l.split('=')[-1][1:-1]
            if "3' bias" in l and "5'-3' bias" not in l:
                bias_3_prim = l.split('=')[-1][1:-1]
            if "5'-3' bias" in l:
                bias_5to3_prim = l.split('=')[-1][1:-1]
                if bias_5to3_prim == "?":
                    bias_5to3_prim = 100.0
    link = f'<a href="../3.QC_files/qualimap/{sample}/qualimapReport.html">map_QC</a>'
    return [Exonic_perc,intronic_perc,intergenic_perc,bias_5_prim,bias_3_prim,bias_5to3_prim,link]        

def read_bamqc(sample):
    with open(f'3.QC_files/qualimap_bamqc/{sample}/genome_results.txt') as Qualimap_file:
        for l in Qualimap_file.readlines():
            try:
                if 'mean insert size' in l:
                    insert_mean = float(''.join(l.split('=')[-1][:-2].split(',')))
                if 'median insert size' in l:
                    insert_median = float(''.join(l.split('=')[-1][:-1].split(',')))
            except ValueError:
                insert_mean, insert_median = 0, 0
    link = f'<a href="../3.QC_files/qualimap_bamqc/{sample}/qualimapReport.html">bam_QC</a>'
    return [insert_mean,insert_median,link]    

def read_bw(sample,dirin,genome_version):
    relative_dir = "/".join(dirin.split("/")[3:])
    link = f'<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db={genome_version}&position=chr1&hgct_customText=track%20type=bigWig%20name={sample}%20description=%22RNA-Seq%20read%20alignment%22%20visibility=full%20bigDataUrl=https://informaticsdata.liai.org/NGS_analyses/ad_hoc/{relative_dir}/3.QC_files/bed_wiggle/{sample}.bw">bigwig</a>'
    return [link]

# Here will call all previous functions and for each sample generate a list that will be inserted into the dataframe
def read_sample(sample,dirin,genome_version):
    return Hawk_smash([read_fastp(sample),read_STAR(sample),read_Qualimap(sample),read_bamqc(sample),read_bw(sample,dirin,genome_version)])

def read_json(file = 'conf_RNA_Seq.json'):
    with open(file) as json_file:
        conf = json.load(json_file)
    return conf

####################### Make reports #######################

dict_conf = read_json(snakemake.params[0])
sample_list = snakemake.params[1]

dirin = dict_conf['config']['dirin']
genome_version = dict_conf['config']['genome_version']
annotation_file = dict_conf['config']['annotation_file']
gene_recovery_perc = dict_conf['QC_threshold']['gene_recovery_perc']
run_file = dict_conf['config']['run_file']
run_ID = dict_conf['config']['run_ID']

df_anno = pd.read_csv(f'{annotation_file}',index_col = 0)

#get count table and TPM
[df_TPM, df_counts,df_counts_filter] = [pd.read_csv(file,index_col = 0) for file in ['4.Output/counts/TPM_counts.csv','4.Output/counts/raw_counts_unfiltered.csv','4.Output/counts/raw_counts.csv']]

# Make report files
df_QC_report = pd.DataFrame([read_sample(sample,dirin,genome_version) for sample in sample_list])
df_QC_report.index = sample_list
df_QC_report.columns = ['total_reads','filtered_reads','filtered_reads_perc','adaptor_trimm_perc','dup_rate','seq_QC','uniquely_mapped_reads','uniquely_mapped_reads_perc','spliced_reads','anno_spliced_reads','too_short_reads','too_short_reads_perc','exonic_perc','intronic_perc','intergenic_perc','bias_5_prim','bias_3_prim','bias_5to3_prim','map_QC','insert_mean','insert_median','bam_QC','bigwig']

df_QC_report['STAR_counts'] =  df_counts[df_QC_report.index].sum(axis = 0)
df_QC_report['STAR_counts_perc'] =  round(100*df_QC_report['STAR_counts']/df_QC_report['uniquely_mapped_reads'],2)
df_QC_report['final_STAR_counts'] = df_counts_filter[df_QC_report.index].sum(axis = 0)
df_QC_report['t_rRNA_counts'] = df_counts[df_QC_report.index].sum(axis = 0) - df_counts_filter[df_QC_report.index].sum(axis = 0)
df_QC_report['t_rRNA_counts_perc'] = round(100*df_QC_report['t_rRNA_counts']/df_QC_report['STAR_counts'],2)

df_genetype = df_anno[['gene_type_4']].merge(df_counts_filter,left_index = True, right_index = True).groupby('gene_type_4').sum().T
df_QC_report['protein_coding_perc'] = (df_genetype['protein_coding']/df_counts_filter.sum(axis = 0)).apply(lambda x: round(100*x,2))
df_QC_report['pseudogene_perc'] = (df_genetype['pseudogene']/df_counts_filter.sum(axis = 0)).apply(lambda x: round(100*x,2))
df_QC_report['long-noncoding_perc'] = (df_genetype['long-noncoding']/df_counts_filter.sum(axis = 0)).apply(lambda x: round(100*x,2))
df_QC_report['short-noncoding_perc'] = (df_genetype['short-noncoding']/df_counts_filter.sum(axis = 0)).apply(lambda x: round(100*x,2))
df_QC_report['Total_genes'] = (df_counts>10).sum(axis = 0)

# Calculate minimal total STAR counts that yields desired percentage of gene recovery rate using Platt func for fitting, then take samples passing the initial filter and get minimal genes per sample based on normal distribution - those < 0.95 confidence interval will be marked 

para = fit_platt([df_QC_report['final_STAR_counts'].values,df_QC_report['Total_genes'].values])
fit_X = np.arange(0,max(df_QC_report['final_STAR_counts'].values),max(df_QC_report['final_STAR_counts'].values)/100)
Y_max = max(platt(para,fit_X,0))
Y_threshold = Y_max*gene_recovery_perc
X_threshold = -np.log(1-gene_recovery_perc)*Y_max/para[1]

data_init_filter = df_QC_report[df_QC_report['final_STAR_counts']>X_threshold]['Total_genes']
para_norm = norm.fit(data_init_filter)
Y_lower = norm.ppf(0.05,para_norm[0],para_norm[1])

# Calculate outliers based on spearman correlation and any sample with mean corr significantly different from the mean at 95% confidence interval will be marked "outlier"

# df_corr = pd.read_csv('4.Output/counts/TPM_counts.csv',index_col = 0).corr(method = 'spearman').fillna(0)    
df_corr = df_TPM.corr(method = 'spearman').fillna(0)
corr_norm = norm.fit(df_corr.reindex(df_QC_report.index).mean(axis = 1))
df_QC_report['Outlier'] = ['Yes' if x < min(norm.ppf(0.05,corr_norm[0],corr_norm[1]),0.6) else 'No' for x in df_corr.reindex(df_QC_report.index).mean(axis = 1)] # if mean spearman corr for a sample < 0.6 or less than 0.95 confidence that is lower than 0.6, then mark as outlier 

df_QC_report['recommendation'] = [RNA_QC(row,dict_conf['QC_threshold'],X_threshold,Y_lower)[0] for row in df_QC_report.itertuples()]
df_QC_report['Note'] = [RNA_QC(row,dict_conf['QC_threshold'],X_threshold,Y_lower)[1] for row in df_QC_report.itertuples()]
df_QC_report[run_ID] = pd.read_csv(run_file,index_col = 0).loc[df_QC_report.index][run_ID]

df_QC_report = df_QC_report[['total_reads','filtered_reads','filtered_reads_perc','adaptor_trimm_perc','dup_rate','uniquely_mapped_reads','uniquely_mapped_reads_perc','spliced_reads','anno_spliced_reads','too_short_reads','too_short_reads_perc','exonic_perc','intronic_perc','intergenic_perc','bias_5_prim','bias_3_prim','bias_5to3_prim','STAR_counts','STAR_counts_perc','t_rRNA_counts','t_rRNA_counts_perc','protein_coding_perc','pseudogene_perc','long-noncoding_perc','short-noncoding_perc','final_STAR_counts','insert_mean','insert_median','Total_genes','seq_QC','map_QC','bam_QC','bigwig','recommendation',run_ID,'Outlier','Note']]

df_QC_report = df_QC_report.sort_index()
df_QC_report[['total_reads','filtered_reads','filtered_reads_perc','adaptor_trimm_perc','dup_rate','uniquely_mapped_reads','uniquely_mapped_reads_perc','spliced_reads','anno_spliced_reads','too_short_reads','too_short_reads_perc','exonic_perc','intronic_perc','intergenic_perc','bias_5_prim','bias_3_prim','bias_5to3_prim','STAR_counts','STAR_counts_perc','t_rRNA_counts','t_rRNA_counts_perc','protein_coding_perc','pseudogene_perc','long-noncoding_perc','short-noncoding_perc','final_STAR_counts','insert_mean','insert_median','Total_genes','recommendation',run_ID,'Outlier','Note']].to_csv('4.Output/QC_report.csv')

df_html = df_QC_report[['seq_QC','map_QC','bam_QC','bigwig','recommendation',run_ID,'Outlier','Note']]
df_html.index.name = f'<a href="./check_QC_PCA.html">PCA_plot</a>\t\t<a href="./QC_report.csv">Download_table</a>\t\t<a href="4.Output/QC_plots">QC_plots</a>'
df_html.to_html(f'4.Output/QC_report.html',escape=False,notebook = True)      