/BioII/lulab_b/chenyinghui/software/conda2/bin/Rscript normalization.R  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data.txt  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  \
NULL \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM.txt  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM_RUVg.txt

/BioII/lulab_b/chenyinghui/software/conda2/bin/Rscript normalization.R  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data_filter.txt  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  \
NULL \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM_filter.txt  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM_RUVg_filter.txt

/BioII/lulab_b/chenyinghui/software/conda2/bin/Rscript diff-exp.R  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data.txt  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  \
NULL \
CRC  \
NC \
NULL \
1 \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/DESeq.txt \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/edgeR.txt


/BioII/lulab_b/chenyinghui/software/conda2/bin/Rscript diff-exp.R  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data_filter.txt  \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  \
NULL \
CRC  \
NC \
NULL \
1 \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/DESeq_filter.txt \
/BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/edgeR_filter.txt

/BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 get_AUC.py --matrix /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_data.txt  --sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  --cv 1 --positive-class CRC --negative-class NC --output-dir /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/result.txt

/BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 get_AUC.py --matrix /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM.txt  --sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  --cv 1 --positive-class CRC --negative-class NC --output-dir /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/result_TMM.txt

/BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 get_AUC.py --matrix /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM_filter.txt  --sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  --cv 1 --positive-class CRC --negative-class NC --output-dir /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/result_TMM_filter.txt

/BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 get_AUC.py --matrix /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM_RUVg.txt  --sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  --cv 1 --positive-class CRC --negative-class NC --output-dir /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/result_TMM_RUVg.txt

/BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 get_AUC.py --matrix /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/count_matrix/exp_TMM_RUVg_filter.txt  --sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/data/label_data.txt  --cv 1 --positive-class CRC --negative-class NC --output-dir /BioII/lulab_b/liuxiaofan/project/expanel/CRC/pico/result/result_TMM_RUVg_filter.txt


