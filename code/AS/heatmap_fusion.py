import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# from code.gene_list import longRNA
from sklearn.preprocessing import MinMaxScaler, StandardScaler
def heatmap_plot(data, datalabel,type, selectfeature,gene,savename,figsize=(13, 13),scaler=True, col_cluster=True, row_cluster=True,yticklabels=False):

    data = data.fillna(0)
    data = data.loc[selectfeature,:]
    gene_name=gene
    if scaler:
        # X = np.log2(data + 0.001).T
        X = data.T
        X = MinMaxScaler().fit_transform(X)
        # X[X>1]=1
        # X[X<-1]=-1
        X_new = pd.DataFrame(X.T)
    else:
        X_new = pd.DataFrame(data)
    X_new.columns = [i for i in range(len(X_new.T))]
    X_new.index=gene_name
    labelnum = len(set(datalabel))
    network_pal = sns.husl_palette(labelnum, s=.45)
    # lut = dict(zip(map(str, datalabel.unique()), network_pal))
    lut = {'NC':'g','CRC_0':'b','CRC_1':'r'}
    col_colors = datalabel.map(lut)

    # plt.switch_backend('agg')
    import sys
    sys.setrecursionlimit(10000)
    sns.set_context("notebook", font_scale=2, rc={"lines.linewidth": 2.5})
    g = sns.clustermap(X_new,
                       col_colors=col_colors,
                       # cmap='flare',
                       cmap='vlag',
                       col_cluster=col_cluster,
                       row_cluster=row_cluster,
                       xticklabels=False,
                       yticklabels=yticklabels,
                       # linewidths=0.05,
                       # center=0,
                       # linewidths=0,
                       figsize=figsize
                       )
    plt.legend(datalabel)
    for label in datalabel.unique():
        g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                                label=label, linewidth=10)
    g.ax_col_dendrogram.legend(loc="upper center", ncol=labelnum)

    # if not os.path.exists(savepath): os.makedirs(savepath)
    # plt.savefig(savename)
    plt.tight_layout()
    # plt.show()
    plt.savefig(savename)
    plt.close()


data_exo_2=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/pico_data.csv',sep='\t')


label = data_exo_2['label']
data_exo_2=data_exo_2[AS]
data_exo=data_exo_2
# data_exo[data_exo>20]=20
data_exo=data_exo.T
AS_=list(pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/feature_0.7.txt',header=None)[0])

fusion_name=[i.split('|')[0]+'|'+i.split('|')[2]+'|'+i.split('|')[5]+'|'+i.split('|')[6] for i in AS]
savename='/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/heatmap_pico.png'
heatmap_plot(data_exo, label, 'CRC', AS_,fusion_name,savename, figsize=(15, 10), scaler=True, col_cluster=False,
             row_cluster=False, yticklabels=True)


