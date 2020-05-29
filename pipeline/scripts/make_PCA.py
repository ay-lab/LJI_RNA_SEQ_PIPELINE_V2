
# general libs
import json
import pandas as pd
import numpy as np

# modeling libs
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

# plotting libs
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go

def read_json(file = 'conf_RNA_Seq.json'):
    with open(file) as json_file:
        conf = json.load(json_file)
    return conf

def Plot_3D(df_meta_input = None,n_comps = 10, dim = 3):
    try:
        df_QC_report = pd.read_csv('4.Output/QC_report.csv',index_col=0)
        
        # Import metadata; if not provided only QC will be used
        try:
            df_meta_plot = pd.read_csv(df_meta_input,index_col=0)
            df_meta_plot = df_meta_plot.select_dtypes('object').merge(df_QC_report[['recommendation','Outlier']],left_index = True,right_index = True)
            df_meta_plot = df_meta_plot.astype(str)
            
        except FileNotFoundError:
            print('No metadata file found')
            df_meta_plot = df_QC_report[['recommendation','Outlier']]
            
        df_plot = pd.read_csv('4.Output/counts/TPM_counts.csv',index_col = 0)[df_meta_plot.index]
        
        # Compute PCA results for plotiing
        data_array = df_plot.apply(lambda x: np.log2(x+1)).values.T
        pca = PCA(n_components = min(n_comps,len(df_QC_report)-1)) # PC number can not exceed sample size -1
        PCA_result = pca.fit_transform(data_array)
        Ratio = pca.explained_variance_ratio_
        mu = PCA_result.mean(axis=0)
        PCA_result = PCA_result - mu
        X = [xx[0] for xx in PCA_result]
        Y = [xx[1] for xx in PCA_result]
        Z = [xx[2] for xx in PCA_result]

        # generate PCA - QC correlation matrix heatmap here
        df_PCA = pd.DataFrame(PCA_result,columns=[f'PC_{i+1}' for i in range(PCA_result.shape[1])],index = df_QC_report.index)
        corr_matrix = np.zeros([int(df_QC_report.shape[1]-3),int(df_PCA.shape[1])])
        p_matrix = np.zeros([int(df_QC_report.shape[1]-3),int(df_PCA.shape[1])])
        for i,col in enumerate(df_QC_report.columns[:-3]): 
            for j,PC in enumerate(df_PCA.columns):
                corr_matrix[i,j] = spearmanr(df_QC_report[col].values,df_PCA[PC].values).correlation
                p_matrix[i,j] = spearmanr(df_QC_report[col].values,df_PCA[PC].values).pvalue
        _,ax = plt.subplots(figsize = (25,6))
        sns.heatmap(pd.DataFrame(p_matrix,index = df_QC_report.columns[:-3], columns= df_PCA.columns).apply(lambda x: -np.log10(x)).T,cmap = 'GnBu')
        for j in range(len(pd.DataFrame(corr_matrix).T.columns)):
            for i in range(len(pd.DataFrame(corr_matrix).T.index)):
                text = ax.text(j+0.5, i+0.5, pd.DataFrame(corr_matrix).T.applymap(lambda x:'%.2f'%x).values[i, j],
                               ha="center", va="center", color="red")
        plt.savefig(f'4.Output/QC_plots/PCA_QC_correlation.png',format='png', bbox_inches='tight')
        
        # Back to plot 3D funciton
        df_PCA = pd.DataFrame({'sample':df_meta_plot.index.values,'X':X,'Y':Y,'Z':Z}).set_index('sample')
        group_factors = df_meta_plot.columns

        scatter_data = []
        Scatter_data_length = []
        if dim == 3:
            for select_factor in group_factors:
                trace = []
                for N in sorted(set(df_meta_plot[select_factor])):
                    trace.append(go.Scatter3d(
                    x = df_PCA[df_meta_plot[select_factor] == N]['X'].values,
                    y = df_PCA[df_meta_plot[select_factor] == N]['Y'].values,
                    z = df_PCA[df_meta_plot[select_factor] == N]['Z'].values,
                    mode='markers',
                    name = N,
                    text=df_meta_plot[df_meta_plot[select_factor] == N].index.values,
                    hoverinfo='text',
                    marker=dict(
                        size=4,
                        line=dict(
                            width=0.5
                        ),
                        opacity=0.8
                    )
                ))

                scatter_data+=trace
                Scatter_data_length.append(len(trace))

        if dim == 2:
            for select_factor in group_factors:
                trace = []
                for N in sorted(set(df_meta_plot[select_factor])):
                    trace.append(go.Scatter(
                    x = df_PCA[df_meta_plot[select_factor] == N]['X'].values,
                    y = df_PCA[df_meta_plot[select_factor] == N]['Y'].values,
                    mode='markers',
                    name = N,
                    text=df_meta_plot[df_meta_plot[select_factor] == N].index.values,
                    hoverinfo='text',
                    marker=dict(
                        size=4,
                        line=dict(
                            width=0.5
                        ),
                        opacity=0.8
                    )
                ))

                scatter_data+=trace
                Scatter_data_length.append(len(trace))

        Vis_table = []
        for i in np.arange(len(Scatter_data_length)):
            Vis_list = []
            for j,elements in enumerate(Scatter_data_length):
                if j == i:
                    Vis_list.append([True]*elements)
                else:
                    Vis_list.append([False]*elements)
            Vis_table.append([item for sublist in Vis_list for item in sublist])
        dict_Vis = dict(zip(group_factors,Vis_table))

        updatemenus = [dict(active=0,
                 buttons=list([   
                    dict(label = Group_factor,
                         method = 'update',
                         args = [ {'visible': dict_Vis[Group_factor]+[True]},{'title': str(Group_factor)}]) for Group_factor in group_factors
                                     ])
                     )
            ]

        layout = go.Layout(
                legend=dict(x=-0.15, y=0),
                scene = dict(
                xaxis=dict(
                title='PC1  '+'%.2f'%Ratio[0],
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                    )
                ),

                yaxis=dict(
                title='PC2  '+'%.2f'%Ratio[1],
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                    )
                ),

                zaxis=dict(
                title='PC3  '+'%.2f'%Ratio[2],
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                    )
                ),),

                updatemenus=updatemenus
                )
        fig = go.Figure(data=scatter_data,layout = layout)
        py.plot(fig,filename='4.Output/check_QC_PCA.html', auto_open=False)
        
    except FileNotFoundError:
        print('Please generate QC report')

dict_conf = read_json(snakemake.params[0])        
meta_input = dict_conf['config']['metadata_dir']
Plot_3D(meta_input,n_comps = 10) 