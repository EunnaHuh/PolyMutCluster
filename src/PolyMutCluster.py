'''
@author: Eunna Huh
'''
import os, sys
import pandas as pd
import seaborn as sns
import os,math,sys, numpy as np
import itertools
from scipy import stats
import multiprocessing as mp
from multiprocessing import current_process
import warnings
import time
from tqdm import tqdm 
from datetime import date
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans
from scipy.cluster import hierarchy
from scipy.spatial import distance
import matplotlib.pyplot as plt
import argparse
today=date.fromtimestamp(time.time())

def parseArguments():
    parser = argparse.ArgumentParser(description='PolyMutCluster')
    parser.add_argument('-In', help='Directory of assay readouts file in csv format', type=str, nargs='+')
    parser.add_argument('-Err', help='Type of error "SEM":standard error of the mean or "SD":standard deviation', type=str, nargs='+')
    parser.add_argument('-Rep', help='Directory of experiment replication file in csv format when type of error is SEM. Default is 4', type=str, nargs='?')
    parser.add_argument('-Norm', help='ND: Normalized Difference, MinMax: normalized each experiment based on min and max', type=str, nargs='+',default='MinMax')
    
    parser.add_argument('-p', help='The number of CPU to be used. Default is to assign all available CPUs', type=int, nargs='?',default=mp.cpu_count())
    parser.add_argument('-itr', help='The nummber of error-propagated matrix generated. Default is 500', type=int, nargs='?',default=500)
    parser.add_argument('-k', help='The maximum number of K. Default is 3', type=int, nargs='?',default=3)
    parser.add_argument('-l', help='scipy.cluster.hierarchy.linkage. Options: single, complete, average, weighted, centroid, ward, median. Default is average' , type=str, nargs='?',default='average')
    parser.add_argument('-pdist', help='A measurement for a pairwise correlation of cluster frequency matrix. Options: braycurtis, canberra, chebyshev, cityblock, correlation, cosine, dice, euclidean, hamming, jaccard, jensenshannon, kulsinski, mahalanobis, matching, minkowski, rogerstanimoto, russellrao, seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule, pearson, kendall, spearman. Default is euclidean', 
                        type=str, nargs='?',default='euclidean')
    parser.add_argument('-log', help='If you want to change scale of measurement from log to linear in the process of normalization, type the names of experiments', type=str, nargs='+')
    parser.add_argument('-S', help='If you want to save every normalized matrices, please type "ON"', type=str, nargs='?',default='OFF')
    args=parser.parse_args()
    args=vars(args)
    return (args) 

def read_csv(inf=None,err=None,rep=None):
    df=pd.read_csv(inf,index_col=0)
    Mut_ls= df.index.tolist()
    M_mat=df[df.columns[[i*2 for i in range(int(df.shape[1]/2))]]].copy()
    Exp_ls=M_mat.columns.tolist()
    std_mat=df[df.columns[[(i*2)+1 for i in range(int(df.shape[1]/2))]]].copy()
    if err =='SEM':
        if rep == None:
            rep_df=pd.DataFrame(np.full((len(df.index),int(df.shape[1]/2)),2),columns=df.columns[[(i*2)+1 for i in range(int(df.shape[1]/2))]],index=df.index)
        else: rep_df=pd.read_csv(rep,index_col=0)**0.5
        std_mat=std_mat.mul(rep_df)
    return (Mut_ls,Exp_ls,M_mat.values,std_mat.values)
    
def err_propagating():
    err_arr=np.zeros(M_mat.shape)
    for i,j in list(itertools.product(list(range(M_mat.shape[0])),list(range(M_mat.shape[1])))):
        if np.isnan(M_mat[i][j]):
            err_arr[i][j]=np.nan
        else:
            LocalProcRandom = np.random.RandomState()
            err_arr[i][j]= LocalProcRandom.normal(M_mat[i][j],std_mat[i][j],1)
    err_mat=pd.DataFrame(err_arr,columns=Exp_ls,index=Mut_ls)
    if log_col != None:
        Exp_log=[s for s in Exp_ls if any (a in s for a in log_col)] 
        err_mat=err_mat.apply(lambda x : 10**x if x.name in Exp_log else x)[Exp_ls]
    
    if Norm_type =='ND':
        Norm_mat=err_mat.T.apply(lambda x: ((x-x['WT'])/(x+x['WT']))+1,axis=1).T.copy()
    if Norm_type =='MinMax':
        Norm_mat=err_mat.copy()
        Norm_mat=(Norm_mat-Norm_mat.min())/(Norm_mat.max()-Norm_mat.min())
    Norm_mat[Norm_mat<0]=0 
    Norm_mat=Norm_mat.fillna(0).round(decimals=5)
    return(Norm_mat)

def Run_clust(iteration=None,num_processers=None):
    print ('Start Clustering ...')
    if num_processers==None:
        num_processers=mp.cpu_count()
    iter_p_processor=math.ceil(iteration/num_processers)
    Freq_1D=mp.Array('d', len(Mut_ls) * len(Mut_ls))
    processers_ls=[]
    for num_p in range(num_processers):
        p=mp.Process(target=_par_processer,args=(num_p,Freq_1D, iter_p_processor))
        p.start()
        processers_ls.append(p)
    for p in processers_ls:
        p.join()
    
    Freq_2D=np.reshape(Freq_1D,(-1,len(Mut_ls)))
    Final_DF=pd.DataFrame(Freq_2D/(iter_p_processor*num_processers),columns=Mut_ls,index=Mut_ls)
    Final_DF.to_csv('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/ClusterFreq_Matrices/AvgClustFrequency.csv'%((today,K,args['itr'])))
    return (Final_DF)

def _par_processer(num_p,Freq_1D, iter_p_processor): 
    current = current_process()
    with tqdm(desc=current.name, total=iter_p_processor, position=current._identity[0] - 1) as progress:
        for n_iter in np.arange(iter_p_processor):
            Norm_mat=err_propagating()
            K_mat=np.zeros((len(Mut_ls),len(Mut_ls)))
            for K_clust in range(2,K+1):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore') 
                    model = NMF(n_components=K, init='random', random_state=0,max_iter=200)
                    W=model.fit_transform(Norm_mat.values)
                    clusterdata=KMeans(n_clusters=K,max_iter=100).fit_predict(W)
                for c1 in range (0,len(clusterdata)):
                    for c2 in range (0,len(clusterdata)):
                        if clusterdata[c1]==clusterdata[c2]:
                            K_mat[c1][c2]=K_mat[c1][c2]+1
            K_avg_mat=K_mat/(K-1)
            K_avg_DF=pd.DataFrame(K_avg_mat,columns=Mut_ls,index=Mut_ls)
            if args['S'] == 'ON':
                Norm_mat.to_csv('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/Normalized_Matrices/Processor-%s_iter-%s_Matrix.csv'%(today,K,args['itr'],num_p+1,n_iter))
                K_avg_DF.to_csv('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/ClusterFreq_Matrices/Processor-%s_iter-%s_ClustFrequency.csv'%(today,K,args['itr'],num_p+1,n_iter))
            Freq_1D[:]+=K_avg_mat.flatten()
            progress.update(1)

def ClusterMap(DF=None):
    sns.set(font='monospace',font_scale=1.2)
    if pdist == 'pearson' or pdist == 'kendall' or pdist == 'spearman' :
        row_linkage = hierarchy.linkage(distance.pdist(np.asarray(DF.corr(method=pdist))), method=linkage) 
        ax=sns.clustermap(DF.corr(method=pdist), col_cluster=True, row_cluster=True, cmap='viridis', linewidth=.75,cbar_kws={'ticks':[-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0]},figsize=(16,16))
    else:
        row_linkage = hierarchy.linkage(distance.pdist(DF,metric=pdist), method=linkage)
        ax=sns.clustermap(DF, col_cluster=True, row_cluster=True, cmap='viridis', linewidth=.75,cbar_kws={'ticks':[0,0.25,0.5,0.75,1.0]},figsize=(16,16))
    plt.savefig('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/Graph/'%(today,K,args['itr'])+'AvgClustFrequency_heatmap.pdf',format='pdf')
    plt.close(ax.fig)
        
    sns.set_style("white")
    plt.figure(figsize=(18,10))
    hierarchy.dendrogram(row_linkage,labels=DF.index,leaf_font_size=12)
    plt.savefig('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/Graph/'%(today,K,args['itr'])+'AvgClustFrequency_dendrogram.png',format='png')
    plt.close()
    
if __name__ == '__main__':
    args=parseArguments()
    K=args['k']
    linkage=args['l']
    pdist=args['pdist']
    Norm_type=args['Norm'][0]
    log_col=args['log']
    
    try: os.mkdir('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/'%(today,K,args['itr']))
    except:''
    try: os.mkdir('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/Graph/'%(today,K,args['itr']))
    except:''
    try: os.mkdir('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/ClusterFreq_Matrices/'%(today,K,args['itr']))
    except:''
    if args['S'] == 'ON':
        try: os.mkdir('/'.join(args['In'][0].split('/')[:-1])+'/%s_Output_K%s_itr%s/Normalized_Matrices/'%(today,K,args['itr']))
        except:''
         
    Mut_ls,Exp_ls,M_mat,std_mat=read_csv(inf=args['In'][0],err=args['Err'][0],rep=args['Rep'])
    Clust_DF=Run_clust(iteration=args['itr'],num_processers=args['p'])
    ClusterMap(DF=Clust_DF)
    
    
     
    