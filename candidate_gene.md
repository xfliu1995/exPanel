## Candidate gene

### 1）筛选流程

这里说明如何进行Candidate gene筛选。我们筛选Candidate gene原则上满足：分类能力好，覆盖度高，异质性低。因此我们围绕这三个方面构造特征。

在Candidate gene筛选过程，我们基本遵从以下过程：

* 高通量测序数据预处理：
    * Mapping得到count matrix或者splicing matrix等
    * count matrix归一化处理
 
 * 筛选指标计算
     * 分类指标
     * 覆盖度指标
     * 异质性指标
 
 * 初筛后画图，主观筛选
     * heatmap，boxplot
     * IGV

### 2）指标计算说明
#### 2.1）  覆盖度指标

ratio_NC：NC样本非0值的比例
ratio_cancer：HCC样本非0值的比例
NC_TMM_mean：NC样本的平均TMM值
HCC_TMM_mean：HCC样本的平均TMM值

脚本：/BioII/lulab_b/liuxiaofan/project/expanel/HCC/exoRbase/code/data.py

#### 2.2） 分类指标和异质性指标

**差异表达计算**

* count：diff-exp.R

>/BioII/lulab_b/chenyinghui/software/conda2/bin/Rscript diff-exp.R  /BioII/lulab_b/liuxiaofan/project/expanel/exoRbase/count_matrix/exp_data.txt  /BioII/lulab_b/liuxiaofan/project/expanel/exoRbase/data/label_data.txt  NULL HCC  Healthy NULL 1 /BioII/lulab_b/liuxiaofan/project/expanel/exoRbase/result/DESeq.txt /BioII/lulab_b/liuxiaofan/project/expanel/exoRbase/result/edgeR.txt




* AS and fusion : diff-exp.py 

>/BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 diff-exp.py --matrix /BioII/lulab_b/liuxiaofan/project/expanel/HCC_mutation/exo_AS/data/matrix.txt  --sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/HCC_mutation/exo_AS/data/sample_classes.txt  --positive-class HCC  --negative-class Healthy --output-dir /BioII/lulab_b/liuxiaofan/project/expanel/HCC_mutation/exo_AS/result/diff_exp.txt



**单个基因分类AUC以及gini index计算**

 * 单基因AUC，单基因每个类别gini index：get_AUC.py 
 
> /BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 get_AUC.py --matrix /BioII/lulab_b/liuxiaofan/project/expanel/HCC_mutation/exo_AS/data/matrix.txt --sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/HCC_mutation/exo_AS/data/sample_classes.txt  --cv 1 --positive-class HCC  --negative-class Healthy --output-dir /BioII/lulab_b/liuxiaofan/project/expanel/HCC_mutation/exo_AS/result/result.txt

>cv取1则计算交叉验证结果（LR,LDA两个分类器），取0仅计算整个AUC（不训练分类器，直接用分布计算）



### 3）画图说明

box plot脚本：/BioII/lulab_b/liuxiaofan/project/expanel/HCC_mutation/result_fusion/box_pico.py

heatmap脚本：/BioII/lulab_b/liuxiaofan/project/expanel/HCC_mutation/result_fusion/heatmap_fusion.py


主要是根据画图调整box坐标轴，并在画heatmap时将异常值进行处理。


### 4）文件说明

在网页中，我们上次了count和AS两种情况下的计算脚本。

