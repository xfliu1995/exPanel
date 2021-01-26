import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

NC_TPM_mean=2
HCC_TPM_mean=2
FDR=0.05
logFC=1
ratio_NC_1=0
ratio_cancer_1=0
ratio_NC_2=0.5
ratio_cancer_2=0.5
AUC=0.8
gini_HCC=1
gini_NC=1

longRNA_list=['tucpRNA', 'lncRNA', 'Y_RNA', 'snRNA', 'pseudogene', 'Mt_tRNA', 'srpRNA', 'snoRNA', 'mRNA']
cirRNA_list=['circRNA']

result_1 = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/exoRbase/result/merge_result_TMM.txt',sep='\t')
result_1['gene_ID']=result_1['feature'].map(lambda x:x.split('|')[0])
result_1['gene_type']=result_1['feature'].map(lambda x:x.split('|')[1])
result_1_=result_1.loc[result_1.gene_type.isin(longRNA_list+cirRNA_list),:]
result_1_new = result_1_.loc[((result_1_['NC_TPM_mean']>NC_TPM_mean)|(result_1_['HCC_TPM_mean']>HCC_TPM_mean))
                             &(result_1_['FDR']<FDR)
                             &((result_1_['logFC']<-logFC)|(result_1_['logFC']>logFC))
                             &((result_1_['ratio_NC']>ratio_NC_1)&(result_1_['ratio_cancer']>ratio_cancer_1))
                             &((result_1_['ratio_NC']>ratio_NC_2)|(result_1_['ratio_cancer']>ratio_cancer_2))
                             &(result_1_['AUC']>AUC)
                             &(result_1_['gini_HCC']<gini_HCC)&(result_1_['gini_NC']<gini_NC),:]
result_1_new =result_1_new.sort_values('AUC',ascending=False)
gene_1 = list(result_1_new['gene_ID'])




result_3 = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/merge_result_TMM.txt',sep='\t')
result_3['gene_ID']=result_3['feature'].map(lambda x:x.split('|')[0])
result_3['gene_type']=result_3['feature'].map(lambda x:x.split('|')[1])
result_3_=result_3.loc[result_3.gene_type.isin(longRNA_list+cirRNA_list),:]
result_3_new = result_3_.loc[((result_3_['NC_TPM_mean']>NC_TPM_mean)|(result_3_['HCC_TPM_mean']>HCC_TPM_mean))
                             &(result_3_['FDR']<FDR)
                             &((result_3_['logFC']<-logFC)|(result_3_['logFC']>logFC))
                             &((result_3_['ratio_NC']>ratio_NC_1)&(result_3_['ratio_cancer']>ratio_cancer_1))
                             &((result_3_['ratio_NC']>ratio_NC_2)|(result_3_['ratio_cancer']>ratio_cancer_2))
                             &(result_3_['AUC']>AUC)
                             &(result_3_['gini_HCC']<gini_HCC)&(result_3_['gini_NC']<gini_NC),:]
result_3_new =result_3_new.sort_values('AUC',ascending=False)
gene_3 = list(result_3_new['gene_ID'])


print('exoRbase')
x1_1=set(gene_1)

print('pico')
x3_1=set(gene_3)


print('exoRbase&pico')
y1_1=set(gene_1)&set(gene_3)


result=[len(x1_1),len(x3_1),len(y1_1)]
print(result)

list(result_3.loc[result_3.gene_ID.isin(list(y1_1)),'feature'].map(lambda x:x.split('|')[1]+'|'+x.split('|')[2]))
list(result_3.loc[result_3.gene_ID.isin(list(y1_1)),'gene_ID'])
