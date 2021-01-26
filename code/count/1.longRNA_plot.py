import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/exoRbase/count_matrix/exp_data.txt',sep="\t")
length_data=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/HCC/code/gene_length.csv',sep='\t')

gene = list(df['feature'].map(lambda x:x.split('|')[0]))
length_data =length_data.set_index('gene')
length=length_data.loc[gene,:]
df=df.set_index('feature')
lengthScaledDf = pd.DataFrame((df.values/length.values.reshape((-1,1))),index=df.index,columns=df.columns)
data_exo = (1000000*lengthScaledDf.div(lengthScaledDf.sum(axis=0))).round(4)

label_exo = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/exoRbase/data/label_data.txt',sep='\t')
sample_id_exo = label_exo.loc[(label_exo['label']=='Healthy')|(label_exo['label']=='CRC'),'sample_id']
NC_exo = label_exo.loc[(label_exo['label']=='Healthy'),'sample_id']
HCC_exo = label_exo.loc[(label_exo['label']=='CRC'),'sample_id']

df = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data.txt',sep="\t")
length_data=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/HCC/code/gene_length.csv',sep='\t')
gene = list(df['feature'].map(lambda x:x.split('|')[0]))
length_data =length_data.set_index('gene')
length=length_data.loc[gene,:]
df=df.set_index('feature')
lengthScaledDf = pd.DataFrame((df.values/length.values.reshape((-1,1))),index=df.index,columns=df.columns)
data_pico = (1000000*lengthScaledDf.div(lengthScaledDf.sum(axis=0))).round(4)

label_pico = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data_stage.txt',sep='\t')
sample_id_pico = label_pico.loc[(label_pico['label']=='NC')|(label_pico['label']=='CRC'),'sample_id']
NC_pico = label_pico.loc[(label_pico['label']=='NC'),'sample_id']
HCC_pico_0 = label_pico.loc[(label_pico['label']=='CRC')&(label_pico['Stage']=='0'),'sample_id']
HCC_pico_1 = label_pico.loc[(label_pico['label']=='CRC')&(label_pico['Stage']=='1'),'sample_id']

data_exo = data_exo.reset_index()
data_exo['gene']=data_exo['feature'].map(lambda x:x.split('|')[0])
data_exo=data_exo.set_index('gene')
exo_data = data_exo.loc[CRC_long,list(NC_exo)+list(HCC_exo)].T
exo_data['label']=['NC' for i in range(len(NC_exo))]+['CRC' for i in range(len(HCC_exo))]

data_pico = data_pico.reset_index()
data_pico['gene']=data_pico['feature'].map(lambda x:x.split('|')[0])
data_pico=data_pico.set_index('gene')
pico_data = data_pico.loc[CRC_long,list(NC_pico)+list(HCC_pico_0)+list(HCC_pico_1)].T
pico_data['label']=['NC' for i in range(len(NC_pico))]+['CRC_0' for i in range(len(HCC_pico_0))]+['CRC_1' for i in range(len(HCC_pico_1))]

exo_data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/result_long/exo_data.csv',sep='\t')
pico_data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/result_long/pico_data.csv',sep='\t')


# sns.set_theme(style="ticks", palette="pastel")
# ax = sns.boxplot(x="label", y="gene", data=data_,order=["NC", "HCC"])
# plt.savefig('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/code/test.png')
# plt.close()