import pandas as pd

diff_exp=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/edgeR.txt',sep='\t')
diff_exp=diff_exp.reset_index()
diff_exp=diff_exp.rename(columns={'index':'feature'})

ratio=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/ratio.txt',sep='\t')
mean_data=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/TPM_mean.txt',sep='\t')

result=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/result.txt',sep='\t')

result_TMM=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/result_TMM.txt',sep='\t')

# result_TMM_RUV=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/result_TMM_RUVg.txt',sep='\t')


data = pd.merge(pd.merge(pd.merge(diff_exp[['feature', 'logFC','PValue', 'FDR']],mean_data[['feature', 'NC_TPM_mean', 'HCC_TPM_mean']],on='feature',how='left')
                ,ratio[['feature','ratio_NC','ratio_cancer']],on='feature',how='left')
                ,result[['feature','gini_NC', 'gini_HCC', 'AUC','AUC_train_LDA', 'AUC_test_LDA', 'AUC_train_LR','AUC_test_LR']],on='feature',how='left')
data_TMM = pd.merge(pd.merge(pd.merge(pd.merge(diff_exp[['feature', 'logFC','PValue', 'FDR']],mean_data[['feature', 'NC_TPM_mean', 'HCC_TPM_mean']],on='feature',how='left')
                ,ratio[['feature','ratio_NC','ratio_cancer']],on='feature',how='left')
                ,result[['feature', 'gini_NC', 'gini_HCC']],on='feature',how='left')
                ,result_TMM[['feature', 'AUC','AUC_train_LDA', 'AUC_test_LDA', 'AUC_train_LR','AUC_test_LR']],on='feature',how='left')
# data_TMM_RUVg = pd.merge(pd.merge(pd.merge(pd.merge(diff_exp[['feature', 'logFC','PValue', 'FDR']],mean_data[['feature', 'NC_TPM_mean', 'HCC_TPM_mean']],on='feature',how='left')
#                 ,ratio[['feature','ratio_NC','ratio_cancer']],on='feature',how='left')
#                 ,result[['feature', 'gini_NC', 'gini_HCC']],on='feature',how='left')
#                 ,result_TMM_RUV[['feature','AUC','AUC_train_LDA', 'AUC_test_LDA', 'AUC_train_LR','AUC_test_LR']],on='feature',how='left')

data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/merge_result.txt',sep='\t',index=False)
data_TMM.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/merge_result_TMM.txt',sep='\t',index=False)
# data_TMM_RUVg.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/merge_result_TMM_RUVg.txt',sep='\t',index=False)

