import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data=pd.read_csv('D:/pythonfile/biology/expanle/exo_data.csv',sep='\t')
# sns.despine(offset=10, trim=True)
ax = sns.boxplot(x="label", y="ENSG00000104760.16", data=data,order=["NC", "HCC"],palette=["m", "g"])
plt.ylabel('FGL1')
# plt.ylim([])
plt.show()
plt.close()

# sns.despine(offset=10, trim=True)
ax = sns.boxplot(x="label", y="ENSG00000104760.16", data=data,order=["NC", "HCC"],palette=["m", "g"])
plt.ylabel('FGL1')
# plt.ylim([])
plt.show()
plt.close()