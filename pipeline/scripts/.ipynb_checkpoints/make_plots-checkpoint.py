# general libs
import json
from functools import reduce
import pandas as pd
import numpy as np

# modeling libs
from scipy import stats
from scipy.stats import spearmanr,norm
from scipy.optimize import least_squares

# plotting libs
import matplotlib.pyplot as plt
import seaborn as sns

####################### basic functions for reading QC files #######################

def read_json(file = 'conf_RNA_Seq.json'):
    with open(file) as json_file:
        conf = json.load(json_file)
    return conf

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

def read_bamcoverage(sample): 
    df = pd.read_csv(f'3.QC_files/qualimap/{sample}/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt',sep = '\t')
    return df[df.columns[1]]*100/sum(df[df.columns[1]])

####################### plot functions #######################

def QC_plot(vis_parameters,dict_threshold):
    
    try:
        df_QC_report = pd.read_csv('4.Output/QC_report.csv').set_index('Unnamed: 0')

        QC_category = list(set(df_QC_report['recommendation']))
        
        colors = sns.color_palette("Set1")


        # Making dot plots for QC parameters and save to png
        
        if (len(QC_category) > 1):
            if ('1.Good' in QC_category): 
                category_length = list(df_QC_report.groupby('recommendation').count().mean(axis = 1)[1:].values)
                category_length = [sum(category_length)/5] + category_length
            else:
                category_length = list(df_QC_report.groupby('recommendation').count().mean(axis = 1).values)
                
                
            f,ax = plt.subplots(1,len(QC_category),figsize = (20,3),gridspec_kw={'width_ratios': category_length},sharey = True)
            for i,category in enumerate(sorted(QC_category)):
                for c,parameter in enumerate(vis_parameters):
                    ax[i].plot(df_QC_report[df_QC_report['recommendation'] == category][parameter],'o',label = parameter,color = colors[c])
                    ax[i].hlines(y = dict_threshold[parameter],xmin = -0.5,xmax = len(df_QC_report[df_QC_report['recommendation'] == category])-0.5,linestyles = 'dashed',colors = colors[c],label = dict_threshold[parameter])

                xlabels = df_QC_report[df_QC_report['recommendation'] == category].index
                ax[i].set_title(category)
            
            
                if ('1.Good' in QC_category):     
                    if i%2 == 0:   
                        ax[i].set_title(category)
                    if i%2 == 1:   
                        ax[i].set_title(category+'\n')
                    if i != 0:
                        ax[i].set_xticklabels(xlabels,rotation = 90)
                    else:
                        ax[0].set_xticks('');
                        ax[0].set_xticklabels('');   
                else:
                    ax[i].set_xticklabels(xlabels,rotation = 90)
                    if i%2 == 0:   
                        ax[i].set_title(category)
                    if i%2 == 1:   
                        ax[i].set_title(category+'\n')
                        
            plt.legend()
            plt.subplots_adjust(wspace=0.0, hspace=0)
            plt.legend(bbox_to_anchor=(1., 1.))
            plt.savefig(f'4.Output/QC_plots/{vis_parameters}.png',format='png', bbox_inches='tight')
        
        elif len(QC_category) == 1:
            f,ax = plt.subplots(figsize = (20,3))
            for c,parameter in enumerate(vis_parameters):
                ax.plot(df_QC_report[parameter],'o',label = parameter,color = colors[c])
                ax.hlines(y = dict_threshold[parameter],xmin = -0.5,xmax = len(df_QC_report)-0.5,linestyles = 'dashed',colors = colors[c],label = dict_threshold[parameter])
            xlabels = df_QC_report.index
            ax.set_title(QC_category[0])
            ax.set_xticklabels(xlabels,rotation = 90)
            plt.legend()
            plt.subplots_adjust(wspace=0.0, hspace=0)
            plt.legend(bbox_to_anchor=(1., 1.))
            plt.savefig(f'4.Output/QC_plots/{vis_parameters}.png',format='png', bbox_inches='tight')
            
    except FileNotFoundError:
        print('Please generate QC report')
                                   
def bam_plot(sample_list):
    try:
        df_QC_report = pd.read_csv('4.Output/QC_report.csv',index_col=0)
        
        # Making dot plots for bam coverage QC and save to png
        QC_category = list(set(df_QC_report['recommendation']))
        df_bamcoverage = pd.DataFrame(columns = sample_list,index = range(100))
        for col in df_bamcoverage.columns:
            df_bamcoverage[col] = read_bamcoverage(col) 
        for category in QC_category:
            sample_short_list = list(df_QC_report[df_QC_report['recommendation'] == category].index)
            df_bamcoverage[sample_short_list].plot(figsize = (8,5))
            plt.title(category)
            plt.legend(bbox_to_anchor=(1., 1.),ncol=max(1,int(len(sample_short_list)/12)))
            plt.savefig(f'4.Output/QC_plots/bam_coverage_{category}.png',format='png', bbox_inches='tight')
            
    except FileNotFoundError:
        print('Please generate QC report')
        
def RNA_QC_spearman():
    '''Using spearman correlation to detect outliers'''
    try:
        # calculate corr for all TPM counts
        df_corr = pd.read_csv('4.Output/counts/TPM_counts.csv',index_col = 0).corr(method = 'spearman')
        sns.clustermap(df_corr,vmin = 0, vmax = 1,figsize = (int(len(df_corr)/3),int(len(df_corr)/3)))
        plt.savefig(f'4.Output/QC_plots/Spearman_correlation.png',format='png', bbox_inches='tight')
        
    except FileNotFoundError:
        print('Please generate TPM table')
        
        
def soft_threshold(dict_conf):
    '''Detects minimal soft threshold for total STAR counts'''
    try:
        gene_recovery_perc = dict_conf['QC_threshold']['gene_recovery_perc']
        df_QC_report = pd.read_csv('4.Output/QC_report.csv',index_col = 0)
        para = fit_platt([df_QC_report['final_STAR_counts'].values,df_QC_report['Total_genes'].values])
        fit_X = np.arange(0,max(df_QC_report['final_STAR_counts'].values),max(df_QC_report['final_STAR_counts'].values)/100)
        Y_max = max(platt(para,fit_X,0))
        Y_threshold = Y_max*gene_recovery_perc # Y percentage threshold based on input
        X_threshold = -np.log(1-gene_recovery_perc)*Y_max/para[1] # Calculated X threshold based on Y threshold value
        
        data_init_filter = df_QC_report[df_QC_report['final_STAR_counts']>X_threshold]['Total_genes']
        para_norm = norm.fit(data_init_filter)
        Y_lower = norm.ppf(0.05,para_norm[0],para_norm[1]) # New Y lower limit based on normal distribution
        
        plt.subplots(figsize = (12,6))
        plt.plot(df_QC_report['final_STAR_counts'].values,df_QC_report['Total_genes'].values,'o')
        X_fit = np.arange(0,max(df_QC_report['final_STAR_counts'].values),max(df_QC_report['final_STAR_counts'].values)/100)
        plt.plot(X_fit,platt(para,X_fit,0),'-')
        plt.axvline(x = max(X_threshold,dict_conf['QC_threshold']['final_STAR_counts']),linestyle = '--',color = 'r')
        plt.axhline(y = Y_lower,linestyle = '--',color = 'r')
        perc_threshold ="{:.1%}".format(gene_recovery_perc)
        
        if X_threshold > int(dict_conf['QC_threshold']['final_STAR_counts']):
            text = f'{int(X_threshold):,}'
            plt.text(X_threshold*1.1, Y_max/2, f'<- {gene_recovery_perc} gene recovery threshold = {text}')
        else:
            X_threshold = int(dict_conf['QC_threshold']['final_STAR_counts'])
            text = f'{int(X_threshold):,}'
            plt.text(int(X_threshold*1.1), Y_max/2, f'<- minimal set threshold = {text}')
            
        plt.text(X_threshold*1.5, Y_max/1.3, f'^ minimal gene threshold of {int(Y_lower)} at 0.95 confidence interval',color = 'r')
        plt.xlabel('STAR counts')
        plt.ylabel('Number of genes')
        plt.savefig(f'4.Output/QC_plots/STAR_minimal_counts_soft_threshold.png',format='png', bbox_inches='tight')
    except FileNotFoundError:
        print('Please generate QC report')

####################### Make plots #######################

dict_conf = read_json(snakemake.params[0])
sample_list = snakemake.params[1]

bam_plot(sample_list)
vis_parameters = ['final_STAR_counts']
QC_plot(vis_parameters,dict_conf['QC_threshold'])

vis_parameters = ['uniquely_mapped_reads_perc','exonic_perc']
QC_plot(vis_parameters,dict_conf['QC_threshold'])

vis_parameters = ['too_short_reads_perc','t_rRNA_counts_perc']
QC_plot(vis_parameters,dict_conf['QC_threshold'])

vis_parameters = ['Total_genes']
QC_plot(vis_parameters,dict_conf['QC_threshold'])

vis_parameters = ['bias_5to3_prim']
QC_plot(vis_parameters,dict_conf['QC_threshold'])

vis_parameters = ['insert_median']
QC_plot(vis_parameters,dict_conf['QC_threshold'])

RNA_QC_spearman()

soft_threshold(dict_conf)