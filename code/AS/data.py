import pandas as pd
import numpy as np
def get_ratio(data,nc,cancer,n=0):
    feature = list(data.index)
    data_NC = data.loc[feature,nc]
    data_cancer = data.loc[feature, cancer]
    ratio = pd.DataFrame(columns={'feature', 'ratio_all', 'ratio_NC', 'ratio_cancer'})
    ratio_all = []
    ratio_NC = []
    ratio_cancer = []
    j = 0
    # n = 0
    for i in feature:
        if j%1000==0:
            print(j)
        j=j+1
        X_all = np.array(data.loc[i,:])
        X_all_ = X_all[X_all>n]
        ratio_all.append(len(X_all_)/len(X_all))
        X_NC = np.array(data_NC.loc[i,:])
        # print(X_NC)
        X_NC_ = X_NC[X_NC>n]
        ratio_NC.append(len(X_NC_)/len(X_NC))
        X_CRC = np.array(data_cancer.loc[i,:])
        X_CRC_ = X_CRC[X_CRC>n]
        ratio_cancer.append(len(X_CRC_)/len(X_CRC))

    ratio['feature']=feature
    ratio['ratio_all']=ratio_all
    ratio['ratio_NC']=ratio_NC
    ratio['ratio_cancer']=ratio_cancer
    return ratio

data_1=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/pico_feature_select/row_data_1020/AS/splicing.txt',sep='\t')
data_1=data_1.fillna(0)
data_1=data_1.rename(columns={'Unnamed: 0':'feature'})
label = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_fusion/data/sample_classes.txt',sep='\t')
sample_id = label.loc[(label['label']=='NC')|(label['label']=='CRC'),'sample_id']
NC = label.loc[(label['label']=='NC'),'sample_id']
CRC = label.loc[(label['label']=='CRC'),'sample_id']
label_new = label.loc[(label['label']=='NC')|(label['label']=='CRC'),:]
data_1 = data_1[['feature']+list(sample_id)]
data_1 = data_1.set_index('feature')
ratio=get_ratio(data_1,NC,CRC,n=0)
ratio.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/ratio.txt',sep='\t',index=False)

data_1=data_1.reset_index()
data_1.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/data/matrix.txt',sep='\t',index=False)
