import pandas as pd
import numpy as np
from scipy import stats
from rpy2 import robjects
import argparse
import warnings
import os
import seaborn as sns
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore") #不显示警告

def BH(p_vulae):
    p_vulae_str = [str(x) for x in p_vulae.tolist()]
    p_vulae_r = 'c('+','.join(p_vulae_str)+')'
    r_script = '''
    p.adjust(%s, "BH")
    '''
    p_BH = np.array(robjects.r(r_script % (p_vulae_r)))

    return  p_BH

def wilcoxontest_feature(X,healthy_number):
    m1 = np.array(X[:healthy_number])
    m2 = np.array(X[healthy_number:])

    if (len(set(m1)) == 1) & (len(set(m2)) == 1):
        p = np.nan
    # elif ((len(set(m1)) == 1) | (len(set(m2)) == 1))&(len(set(m1))+len(set(m2))>2):
    #     t, p = stats.ranksums(m1, m2)
    else:
        t, p = stats.mannwhitneyu(m1, m2)
    p_vulae = p
    return p_vulae

def logFC_test(X,healthy_number):
    m1 = np.array(X[:healthy_number])
    m2 = np.array(X[healthy_number:])
    if m1.mean() == 0: m1_mean = m1.mean() + 0.01
    else: m1_mean = m1.mean()
    if m2.mean() == 0: m2_mean = m2.mean() + 0.01
    else: m2_mean = m2.mean()

    logFC = np.log2(m2_mean/m1_mean)

    return logFC


def main(args):

    data=pd.read_csv(args.matrix,sep='\t')
    label = pd.read_csv(args.sample_classes,sep='\t')
    sample_id = label.loc[(label['label']==args.negative_class)|(label['label']==args.positive_class),'sample_id']
    NC = list(label.loc[(label['label']==args.negative_class),'sample_id'])
    HCC = list(label.loc[(label['label']==args.positive_class),'sample_id'])
    data=data.set_index(data.columns[0])
    data=data[NC+HCC]
    healthy_number=len(NC)

    wilcoxon = pd.DataFrame(columns=['ID', 'P', 'FDR', 'logFC'])
    p_wilcoxon = data.apply(lambda x: wilcoxontest_feature(x, healthy_number), axis=1)
    logFC = data.apply(lambda x: logFC_test(x, healthy_number), axis=1)
    wilcoxon['ID'] = data.index
    wilcoxon['P'] = p_wilcoxon.values
    wilcoxon['logFC'] = logFC.values
    wilcoxon = wilcoxon.loc[pd.notnull(wilcoxon['P']), :]
    p_wilcoxon_bh = BH(np.array(wilcoxon['P']))
    wilcoxon['FDR'] = p_wilcoxon_bh
    wilcoxon['FDR'] = wilcoxon['FDR'].astype('float')
    wilcoxon.to_csv(args.output_dir, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Machine learning module')
    parser.add_argument('--matrix', type=str, required=True,
        help='input feature matrix (rows are samples and columns are features',dest='matrix')
    parser.add_argument('--sample-classes', type=str, required=True,
        help='input file containing sample classes with 2 columns: sample_id, sample_class',dest='sample_classes')
    parser.add_argument('--positive-class', type=str,
        help='comma-separated list of sample classes to use as positive class',dest='positive_class')
    parser.add_argument('--negative-class', type=str,
        help='comma-separates list of sample classes to use as negative class',dest='negative_class')
    parser.add_argument('--output-dir', '-o', type=str, metavar='DIR',
        required=True, help='output directory',dest='output_dir')
    args = parser.parse_args()
    main(args)