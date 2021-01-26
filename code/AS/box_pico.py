import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data=pd.read_csv('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/pico_data.csv',sep='\t')


plt.figure(figsize=(4,5))
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2.5})
ax = sns.boxplot(x="label", y="SE|ENSG00000087191.12|PSMC5|chr17|+|63827828|63827923|63827445|63827514|63828137|63828209", data=data,order=["NC", "CRC_0","CRC_1"],palette=["#808080", "#239B56","#2980B9"])
plt.ylabel('SE|PSMC5|63827828|63827923')
# plt.ylim([0,400])
plt.tight_layout()
sns.despine(offset=5, trim=True)
plt.savefig('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/SE|PSMC5|63827828|63827923.png')
plt.close()

plt.figure(figsize=(4,5))
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2.5})
ax = sns.boxplot(x="label", y="SE|ENSG00000087191.12|PSMC5|chr17|+|63827673|63827923|63827459|63827514|63828137|63828209", data=data,order=["NC", "CRC_0","CRC_1"],palette=["#808080", "#239B56","#2980B9"])
plt.ylabel('SE|PSMC5|63827673|63827923')
# plt.ylim([0,400])
plt.tight_layout()
sns.despine(offset=5, trim=True)
plt.savefig('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/SE|PSMC5|63827673|63827923.png')
plt.close()

plt.figure(figsize=(4,5))
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2.5})
ax = sns.boxplot(x="label", y="SE|ENSG00000087191.12|PSMC5|chr17|+|63827850|63827923|63827459|63827514|63828137|63828209", data=data,order=["NC", "CRC_0","CRC_1"],palette=["#808080", "#239B56","#2980B9"])
plt.ylabel('SE|PSMC5|63827850|63827923')
# plt.ylim([0,400])
plt.tight_layout()
sns.despine(offset=5, trim=True)
plt.savefig('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/SE|PSMC5|63827850|63827923.png')
plt.close()

plt.figure(figsize=(4,5))
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2.5})
ax = sns.boxplot(x="label", y="MXE|ENSG00000164924.17|YWHAZ|chr8|-|100924915|100925039|100948595|100948900|100924134|100924298|100953268|100953341", data=data,order=["NC", "CRC_0","CRC_1"],palette=["#808080", "#239B56","#2980B9"])
plt.ylabel('MXE|YWHAZ|100924915|100925039')
# plt.ylim([0,400])
plt.tight_layout()
sns.despine(offset=5, trim=True)
plt.savefig('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/MXE|YWHAZ|100924915|100925039.png')
plt.close()

plt.figure(figsize=(4,5))
sns.set_context("talk", font_scale=1, rc={"lines.linewidth": 2.5})
ax = sns.boxplot(x="label", y="SE|ENSG00000082074.15|FYB1|chr5|-|39212677|39213017|39202773|39202987|39219442|39219563", data=data,order=["NC", "CRC_0","CRC_1"],palette=["#808080", "#239B56","#2980B9"])
plt.ylabel('SE|FYB1|39212677|39213017')
# plt.ylim([0,400])
plt.tight_layout()
sns.despine(offset=5, trim=True)
plt.savefig('/BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/result_AS/SE|FYB1|39212677|39213017.png')
plt.close()

