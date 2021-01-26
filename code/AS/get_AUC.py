import pandas as pd
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics
from sklearn.metrics import roc_curve,roc_auc_score,auc,precision_recall_curve,average_precision_score,accuracy_score
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVC, SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model  import LogisticRegression
from sklearn.utils.class_weight import compute_sample_weight
from numba import jit
import argparse, sys, os, errno

# @jit(nopython=False)
def gini_index(x, ratio_nan=0.9):
    ########
    x_nonull = x[np.isnan(x) == False]
    ratio = len(x_nonull) / len(x)
    if ratio <= ratio_nan:
        gini = np.nan
    else:
        sorted = np.sort(x_nonull)
        height, area = 0, 0
        for i in range(0, len(sorted)):
            height += sorted[i]
            area += height - sorted[i] / 2.
        fair_area = height * len(sorted) / 2.
        if fair_area == 0:
            gini = np.nan
        else:
            gini = (fair_area - area) / fair_area
    return gini

def get_gini(x, y,ratio_nan=0.9):
    x_nc = x[y==0]
    x_hcc=x[y==1]
    gini_nc=gini_index(x_nc, ratio_nan=0.9)
    gini_hcc = gini_index(x_hcc, ratio_nan=0.9)
    return gini_nc,gini_hcc
# @jit(nopython=False)
def getprob(x, y):
    data = x.copy()
    data = pd.DataFrame(data)
    data['label'] = y
    datamean = data.groupby('label').apply(lambda x: x.mean())
    if (x.max() - x.min())==0:
        prob = np.zeros([len(x),1])
    else:
        if datamean.loc[datamean['label'] == 1, :].iloc[0, 0] > \
                datamean.loc[datamean['label'] == 0, :].iloc[0, 0]:
            prob = (x - x.min()) / (x.max() - x.min())
        else:
            prob = 1 - (x - x.min()) / (x.max() - x.min())
    return prob
# @jit(nopython=False)
def plot_roc(prob, y):
    fpr, tpr, _ = roc_curve(y, prob)
    precision, recall, _ = precision_recall_curve(y, prob)  # recall: Identical to sensitivity,TPR; precision: PPV
    roc_auc = metrics.auc(fpr, tpr)
    return roc_auc


# @jit(nopython=False)
def auc(x,y):
    if len(set(list(x)))==1:
        auc=np.nan
    else:
        prob=getprob(x, y)
        auc=plot_roc(prob, y)
    return auc

# @jit(nopython=False)
def model(x_train,x_test,y_train,y_test):
    if len(set(list(x_train[:,0])))==1:
        auc_train_1, auc_test_1, auc_train_2, auc_test_2=np.nan,np.nan,np.nan,np.nan
    else:
        clf_1 = LinearDiscriminantAnalysis()
        clf_2 = LogisticRegression(penalty='l2', solver='liblinear', C=1)
        clf_1.fit(x_train,y_train)
        clf_2.fit(x_train,y_train)
        predict_train_1 = clf_1.predict_proba(x_train)
        predict_test_1 = clf_1.predict_proba(x_test)
        fpr, tpr, thresholds = metrics.roc_curve(y_train, predict_train_1[:, 1], pos_label=1)
        auc_train_1 = metrics.auc(fpr, tpr)
        fpr, tpr, thresholds = metrics.roc_curve(y_test, predict_test_1[:, 1], pos_label=1)
        auc_test_1 = metrics.auc(fpr, tpr)
        predict_train_2 = clf_2.predict_proba(x_train)
        predict_test_2 = clf_2.predict_proba(x_test)
        fpr, tpr, thresholds = metrics.roc_curve(y_train, predict_train_2[:, 1], pos_label=1)
        auc_train_2 = metrics.auc(fpr, tpr)
        fpr, tpr, thresholds = metrics.roc_curve(y_test, predict_test_2[:, 1], pos_label=1)
        auc_test_2 = metrics.auc(fpr, tpr)

    return auc_train_1,auc_test_1,auc_train_2,auc_test_2


def main(args):

    data=pd.read_csv(args.matrix,sep='\t')
    label = pd.read_csv(args.sample_classes,sep='\t')

    feature = data['feature']
    data = data.set_index('feature')
    sample= list(data.columns)
    label=label.set_index('sample_id')
    label=label.loc[sample]

    X=np.array(data.T).reshape(len(label),len(feature))
    Y=np.array(label['label'].replace(args.negative_class,0).replace(args.positive_class,1)).reshape(len(label),)
    if args.cv==1:
        result=pd.DataFrame(columns=['feature','AUC_train_LDA','AUC_test_LDA','AUC_train_LR','AUC_test_LR','gini_NC','gini_HCC','AUC'])
        result['feature']=feature
        for j in range(len(feature)):
            result.loc[j, 'gini_NC'], result.loc[j, 'gini_HCC'] = get_gini(X[:, j], Y)
            result.loc[j, 'AUC'] = auc(X[:,j],Y)
            if j%1000==0:
                print(j)
            auc_train_LDA=[]
            auc_test_LDA=[]
            auc_train_LR=[]
            auc_test_LR=[]
            for n in range(0, 1):
                skf = StratifiedKFold(n_splits=5, random_state=n, shuffle=True)
                for train, test in skf.split(X, Y):
                    x_train = X[train,j].reshape(len(train),1)
                    y_train = Y[train]
                    x_test = X[test, j].reshape(len(test),1)
                    y_test = Y[test]
                    auc_train_LDA_, auc_test_LDA_, auc_train_LR_, auc_test_LR_ = model(x_train,x_test,y_train,y_test)
                    auc_train_LDA.append(auc_train_LDA_)
                    auc_test_LDA.append(auc_test_LDA_)
                    auc_train_LR.append(auc_train_LR_)
                    auc_test_LR.append(auc_test_LR_)
            auc_train_LDA=np.array(auc_train_LDA)
            auc_test_LDA=np.array(auc_test_LDA)
            auc_train_LR=np.array(auc_train_LR)
            auc_test_LR=np.array(auc_test_LR)
            result.loc[j,'AUC_train_LDA']= np.mean(auc_train_LDA[np.isnan(auc_train_LDA) == False])
            result.loc[j,'AUC_test_LDA'] = np.mean(auc_test_LDA[np.isnan(auc_test_LDA) == False])
            result.loc[j,'AUC_train_LR']= np.mean(auc_train_LR[np.isnan(auc_train_LR) == False])
            result.loc[j,'AUC_test_LR'] = np.mean(auc_test_LR[np.isnan(auc_test_LR) == False])
        result.to_csv(args.output_dir,sep='\t',index=False)
    else:
        result=pd.DataFrame(columns=['feature','gini_NC','gini_HCC','AUC'])
        result['feature']=feature
        for j in range(len(feature)):
            if j%1000==0:
                print(j)
            result.loc[j, 'gini_NC'],result.loc[j, 'gini_HCC'] = get_gini(X[:,j],Y)
            result.loc[j, 'AUC'] = auc(X[:,j],Y)
        result.to_csv(args.output_dir, sep='\t', index=False)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Machine learning module')
    parser.add_argument('--matrix', type=str, required=True,
        help='input feature matrix (rows are samples and columns are features',dest='matrix')
    parser.add_argument('--sample-classes', type=str, required=True,
        help='input file containing sample classes with 2 columns: sample_id, sample_class',dest='sample_classes')
    parser.add_argument('--cv', type=int, required=True,
        help='cv',dest='cv')
    parser.add_argument('--positive-class', type=str,
        help='comma-separated list of sample classes to use as positive class',dest='positive_class')
    parser.add_argument('--negative-class', type=str,
        help='comma-separates list of sample classes to use as negative class',dest='negative_class')
    parser.add_argument('--output-dir', '-o', type=str, metavar='DIR',
        required=True, help='output directory',dest='output_dir')
    args = parser.parse_args()
    main(args)
