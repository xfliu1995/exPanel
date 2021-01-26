/BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 diff-exp.py \
--matrix /BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/data/matrix.txt  \
--sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/data/sample_classes.txt  \
--positive-class CRC  \
--negative-class NC \
--output-dir /BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/diff_exp.txt



/BioII/lulab_b/liuxiaofan/software/conda3/bin/python3 get_AUC.py \
--matrix /BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/data/matrix.txt  \
--sample-classes /BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/data/sample_classes.txt  \
--cv 1 \
--positive-class CRC  \
--negative-class NC \
--output-dir /BioII/lulab_b/liuxiaofan/project/expanel/CRC_mutation/pico_AS/result/result.txt



