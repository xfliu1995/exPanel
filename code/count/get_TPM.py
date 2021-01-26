import numpy as np
import pandas as pd
import argparse
import sys
print("Load data ...")
df = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data.txt',sep="\t")
length_data=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/HCC/code/gene_length.csv',sep='\t')

print("Done .")
print("Calculate TPM ...")
gene = list(df['feature'].map(lambda x:x.split('|')[0]))
length_data =length_data.set_index('gene')
length=length_data.loc[gene,:]
df=df.set_index('feature')
lengthScaledDf = pd.DataFrame((df.values/length.values.reshape((-1,1))),index=df.index,columns=df.columns)
data_1 = (1000000*lengthScaledDf.div(lengthScaledDf.sum(axis=0))).round(4)
label = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt',sep='\t')

NC = label.loc[(label['label']=='NC'),'sample_id']
CRC = label.loc[(label['label']=='CRC'),'sample_id']

data_NC=data_1[NC]
data_HCC=data_1[CRC]

data_NC_mean = data_NC.mean(axis=1)
data_NC_mean = pd.DataFrame(data_NC_mean)
data_NC_mean=data_NC_mean.reset_index()
data_NC_mean.columns=['feature','NC_TPM_mean']
data_HCC_mean = data_HCC.mean(axis=1)
data_HCC_mean = pd.DataFrame(data_HCC_mean)
data_HCC_mean=data_HCC_mean.reset_index()
data_HCC_mean.columns=['feature','HCC_TPM_mean']

data=pd.merge(data_NC_mean,data_HCC_mean,on='feature',how='inner')
data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/TPM_mean.txt',sep='\t',index=False)

