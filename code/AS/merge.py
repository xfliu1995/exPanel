import pandas as pd

diff_exp=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/diff_exp.txt',sep='\t')
diff_exp=diff_exp.rename(columns={'ID':'feature'})

ratio=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/ratio.txt',sep='\t')

result=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/result.txt',sep='\t')



data = pd.merge(pd.merge(diff_exp[['feature', 'logFC','P', 'FDR']],ratio[['feature','ratio_NC','ratio_cancer']],on='feature',how='left')
                ,result[['feature','gini_NC', 'gini_HCC', 'AUC','AUC_train_LDA', 'AUC_test_LDA', 'AUC_train_LR','AUC_test_LR']],on='feature',how='left')

data.to_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/merge_result.txt',sep='\t',index=False)

