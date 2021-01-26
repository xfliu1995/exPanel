import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

P=0.05
FDR=1
logFC=1
ratio_NC_1=0
ratio_cancer_1=0
ratio_NC_2=0
ratio_cancer_2=0
AUC=0.7
gini_CRC=1
gini_NC=1

result_1 = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/merge_result.txt',sep='\t')
result_1_new = result_1.loc[(result_1['P']<P)&(result_1['FDR']<FDR)
                             &((result_1['logFC']<-logFC)|(result_1['logFC']>logFC))
                             &((result_1['ratio_NC']>ratio_NC_1)&(result_1['ratio_cancer']>ratio_cancer_1))
                             &((result_1['ratio_NC']>ratio_NC_2)|(result_1['ratio_cancer']>ratio_cancer_2))
                             &(result_1['AUC']>AUC)
                             &(result_1['gini_HCC']<gini_CRC)&(result_1['gini_NC']<gini_NC),:]
result_1_new =result_1_new.sort_values('AUC',ascending=False)
feature_1 = list(result_1_new['feature'])

result_2 = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/exo_AS/result/merge_result.txt',sep='\t')
result_2_new = result_2.loc[(result_2['P']<P)&(result_2['FDR']<FDR)
                             &((result_2['logFC']<-logFC)|(result_2['logFC']>logFC))
                             &((result_2['ratio_NC']>ratio_NC_1)&(result_2['ratio_cancer']>ratio_cancer_1))
                             &((result_2['ratio_NC']>ratio_NC_2)|(result_2['ratio_cancer']>ratio_cancer_2))
                             &(result_2['AUC']>AUC)
                             &(result_2['gini_HCC']<gini_CRC)&(result_2['gini_NC']<gini_NC),:]
result_2_new =result_2_new.sort_values('AUC',ascending=False)
feature_2 = list(result_2_new['feature'])

print('pico')
x1_1=set(feature_1)
print(len(x1_1))
print('exo')
x1_2=set(feature_2)
print(len(x1_2))
print('pico&exo')
x1_3=set(feature_1)&set(feature_2)
print(len(x1_3))

pd.DataFrame(feature_1).to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/feature_0.7.txt',sep='\t',index=False,header=False)
