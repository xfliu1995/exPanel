import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/data/matrix.txt',sep="\t")

label_pico = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_fusion/data/label_data_stage.txt',sep='\t')
sample_id_pico = label_pico.loc[(label_pico['label']=='NC')|(label_pico['label']=='CRC'),'sample_id']
NC_pico = label_pico.loc[(label_pico['label']=='NC'),'sample_id']
CRC_pico_0 = label_pico.loc[(label_pico['label']=='CRC')&(label_pico['Stage']=='0'),'sample_id']
CRC_pico_1 = label_pico.loc[(label_pico['label']=='CRC')&(label_pico['Stage']=='1'),'sample_id']

AS_=list(pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/feature_0.7.txt',header=None)[0])
data_pico=df.set_index('feature')
pico_data = data_pico.loc[AS,list(NC_pico)+list(CRC_pico_0)+list(CRC_pico_1)].T
pico_data['label']=['NC' for i in range(len(NC_pico))]+['CRC_0' for i in range(len(CRC_pico_0))]+['CRC_1' for i in range(len(CRC_pico_1))]

pico_data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/pico_data.csv',sep='\t')


# sns.set_theme(style="ticks", palette="pastel")
# ax = sns.boxplot(x="label", y="gene", data=data_,order=["NC", "CRC"])
# plt.savefig('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/code/test.png')
# plt.close()