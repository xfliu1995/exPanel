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

data_1=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/all-sample-counts.txt',sep='\t')
data_1=data_1.fillna(0)
label = pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/sample_classes.txt',sep='\t')
sample_id = label.loc[(label['label']=='NC')|(label['label']=='CRC'),'sample_id']
NC = label.loc[(label['label']=='NC'),'sample_id']
CRC = label.loc[(label['label']=='CRC'),'sample_id']
label_new = label.loc[(label['label']=='NC')|(label['label']=='CRC'),:]
data_1 = data_1[['feature']+list(sample_id)]
data_1 = data_1.set_index('feature')
ratio=get_ratio(data_1,NC,CRC,n=0)
ratio.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/ratio.txt',sep='\t',index=False)
p=0.1
ratio_filiter = ratio.loc[(ratio['ratio_NC']>p)|(ratio['ratio_cancer']>p),:]
data = data_1.loc[list(ratio_filiter['feature']),:]

data_1.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data.txt',sep='\t',index=True)
data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data_filter.txt',sep='\t',index=True)
label_new.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt',sep='\t',index=False)


data=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM.txt',sep='\t')
data=data.reset_index()
data=data.rename(columns={'index':'feature'})
data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM.txt',sep='\t',index=False)
