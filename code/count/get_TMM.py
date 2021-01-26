import pandas as pd
import numpy as np

data_1=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/pico/count_matrix/exp_TMM.txt',sep='\t')
label = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/pico/data/label_data.txt',sep='\t')
data_1 =data_1.set_index('feature')
NC = label.loc[(label['label']=='NC'),'sample_id']
HCC = label.loc[(label['label']=='HCC'),'sample_id']

data_NC=data_1[NC]
data_HCC=data_1[HCC]

data_NC_mean = data_NC.mean(axis=1)
data_NC_mean = pd.DataFrame(data_NC_mean)
data_NC_mean=data_NC_mean.reset_index()
data_NC_mean.columns=['feature','NC_TMM_mean']
data_HCC_mean = data_HCC.mean(axis=1)
data_HCC_mean = pd.DataFrame(data_HCC_mean)
data_HCC_mean=data_HCC_mean.reset_index()
data_HCC_mean.columns=['feature','HCC_TMM_mean']

data=pd.merge(data_NC_mean,data_HCC_mean,on='feature',how='inner')
data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/pico/count_matrix/TMM_mean.txt',sep='\t',index=False)
