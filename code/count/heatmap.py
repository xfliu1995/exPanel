import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def heatmap_plot(data, datalabel,type, selectfeature,savename,figsize=(13, 13),scaler=True, col_cluster=True, row_cluster=True,yticklabels=False):
    from sklearn.preprocessing import MinMaxScaler
    data = data.fillna(0)
    data = data.loc[selectfeature,:]
    gene_name=[i.split('|')[1]+'|'+i.split('|')[2] for i in selectfeature]
    if scaler:
        # X = np.log2(data + 0.001).T
        X = data.T
        X = MinMaxScaler().fit_transform(X)
        X_new = pd.DataFrame(X.T)
    else:
        X_new = pd.DataFrame(data)
    X_new.columns = [i for i in range(len(X_new.T))]
    X_new.index=gene_name
    labelnum = len(set(datalabel))
    network_pal = sns.husl_palette(labelnum, s=.45)
    # lut = dict(zip(map(str, datalabel.unique()), network_pal))
    lut = {'Healthy':'g',type:'r'}
    col_colors = datalabel.map(lut)

    # plt.switch_backend('agg')
    import sys
    sys.setrecursionlimit(10000)
    sns.set_context("notebook", font_scale=2, rc={"lines.linewidth": 2.5})
    g = sns.clustermap(X_new,
                       col_colors=col_colors,
                       cmap='vlag',
                       col_cluster=col_cluster,
                       row_cluster=row_cluster,
                       xticklabels=False,
                       yticklabels=yticklabels,
                       # linewidths=0.05,
                       linewidths=0,
                       figsize=figsize)
    plt.legend(datalabel)
    for label in datalabel.unique():
        g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                                label=label, linewidth=10)
    g.ax_col_dendrogram.legend(loc="upper center", ncol=labelnum)

    # if not os.path.exists(savepath): os.makedirs(savepath)
    plt.savefig(savename)
    plt.close()

######筛选基因##################
result = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/pico/result/merge_result_TMM.txt',sep='\t')
result_new = result.loc[((result['NC_TPM_mean']>2)|(result['HCC_TPM_mean']>2))&(result['FDR']<0.05)&
                        ((result['logFC']<-1)|(result['logFC']>1))&((result['ratio_NC']>0.1)|(result['ratio_cancer']>0.1))
                        &(result['AUC']>0.8)&((result['gini_NC']<0.3)&(result['gini_HCC']<0.3)),:]

# result_new = result.loc[((result['NC_TPM_mean']>2)|(result['HCC_TPM_mean']>2))&(result['FDR']<0.05)&
#                         ((result['logFC']<-1)|(result['logFC']>1))&((result['ratio_NC']>0.1)|(result['ratio_cancer']>0.1))
#                         &(result['AUC']>0.9),:]
result_new =result_new.sort_values('AUC',ascending=False)
gene = list(result_new['feature'])
print(len(gene))
gene_new=[
# 'ENSG00000153147.5|mRNA|SMARCA5|ENSG00000153147.5|ENSG00000153147.5|0|7923',
# 'ENSG00000255717.6|lncRNA|SNHG1|ENSG00000255717.6|ENSG00000255717.6|0|2768',
'ENSG00000169045.17|mRNA|HNRNPH1|ENSG00000169045.17|ENSG00000169045.17|0|5564',
'ENSG00000047597.5|mRNA|XK|ENSG00000047597.5|ENSG00000047597.5|0|5206',
'ENSG00000222041.10|lncRNA|CYTOR|ENSG00000222041.10|ENSG00000222041.10|0|1096',
# 'ENSG00000227373.5|lncRNA|AL121983.2|ENSG00000227373.5|ENSG00000227373.5|0|688',
'ENSG00000265185.5|snoRNA|SNORD3B-1|ENSG00000265185.5|ENSG00000265185.5|0|758',
# 'ENSG00000164362.18|mRNA|TERT|ENSG00000164362.18|ENSG00000164362.18|0|4018',
'ENSG00000214049.6|lncRNA|UCA1|ENSG00000214049.6|ENSG00000214049.6|0|2299',
'ENSG00000251562.7|lncRNA|MALAT1|ENSG00000251562.7|ENSG00000251562.7|0|8708',
'ENSG00000224078.13|lncRNA|SNHG14|ENSG00000224078.13|ENSG00000224078.13|0|13074',
'G053567|tucpRNA|G053567|G053567|G053567|0|27939',
'ENSG00000251164.1|lncRNA|HULC|ENSG00000251164.1|ENSG00000251164.1|0|556',
'ENSG00000196620.9|mRNA|UGT2B15|ENSG00000196620.9|ENSG00000196620.9|0|2077',
'ENSG00000109072.13|mRNA|VTN|ENSG00000109072.13|ENSG00000109072.13|0|2149',
'ENSG00000130649.9|mRNA|CYP2E1|ENSG00000130649.9|ENSG00000130649.9|0|7918'
]
#######heatmap##################
data=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/pico/count_matrix/exp_TMM.txt',sep='\t')
label = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/pico/data/label_data.txt',sep='\t')
NC = list(label.loc[(label['label']=='NC'),'sample_id'])
HCC = list(label.loc[(label['label']=='HCC'),'sample_id'])
label_ = pd.Series(['Healthy' for i in range(len(NC))] + ['HCC' for i in range(len(HCC))])
data = data.set_index('feature')
data = data[NC+HCC]
# savename='/BioII/lulab_b/liuxiaofan/project/expanel/pico/result/heatmap_0.8.png'
# heatmap_plot(data, label_,'HCC', gene,savename,figsize=(13, 10),scaler=True, col_cluster=False, row_cluster=True,yticklabels=True)
savename='/BioII/lulab_b/liuxiaofan/project/expanel/pico/result/heatmap_select.png'
heatmap_plot(data, label_,'HCC', gene_new,savename,figsize=(15, 10),scaler=True, col_cluster=False, row_cluster=False,yticklabels=True)

result_2=result.loc[result.feature.isin(gene_new),:]
result_2=result_2.set_index('feature')
result_2 =result_2.loc[gene_new,:]
result_2.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/pico/result/result_selectfeature.txt',sep='\t',index=True)






