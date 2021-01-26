#! /usr/bin/env Rscript

library(org.At.tair.db)
args<-commandArgs(TRUE)
matrix <-args[1]
classes <-args[2]
samples <- args[3]
positive_class <-args[4]
negative_class <-args[5]
batch <-args[6]
batch_index <-args[7]
output_file_DESeq <-args[8]
output_file_edgeR <-args[9]

message('read count matrix: ', matrix)
mat <- read.table(matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')
message('read class information: ', classes)
class_info <- read.table(classes, row.names=1, check.names=FALSE, header = TRUE, sep='\t', as.is=TRUE)
class_info <- class_info[colnames(mat),]
names(class_info) <- colnames(mat)

if(samples=='NULL'){
  class_info <- class_info
}else{
  message('read sample ids: ', samples)
  samples <- read.table(samples, check.names=FALSE, header=FALSE, as.is=TRUE)[,1]
  mat <- mat[, samples]
  class_info <- class_info[samples]
}
# get positive and negative class
for(cls in strsplit(positive_class, ',')[[1]]){
  class_info[class_info == cls] <- 'positive'
}
positive_samples <- names(class_info)[class_info == 'positive']
if(length(positive_samples) == 0){
  stop('No positive samples found')
}
message('Number of positive samples: ', length(positive_samples))
negative_samples <- NULL
for(cls in strsplit(negative_class, ',')[[1]]){
  class_info[class_info == cls] <- 'negative'
}
negative_samples <- names(class_info)[class_info == 'negative']
if(length(negative_samples) == 0){
  stop('No negative samples found')
}
message('Number of negative samples: ', length(negative_samples))

samples <- c(positive_samples, negative_samples)

group <- class_info[samples]
mat <- as.matrix(mat[,samples])
#class_info <- as.matrix(class_info)
#colnames(class_info) <- 'label'
mat <- as.matrix(mat)

# read batch information
if(batch=='NULL'){
  batch <- NULL
}else{
  message('read batch information from: ', batch)
  batch <- read.table(batch, check.names=FALSE, header=TRUE, as.is=TRUE, row.names=1, sep='\t')
  if((batch_index < 1) || (batch_index > ncol(batch))){
    stop('Batch index out of bound')
  }
  batch <- batch[samples, batch_index]
}


# Required columns for a differential expression file: baseMean, log2FoldChange, pvalue, padj

message('DESeq2')
suppressPackageStartupMessages(library(DESeq2))
if(!is.null(batch)){
  # include batch into regression
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = data.frame(group=group, batch=batch),
                                design = ~group + batch)
} else {
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = data.frame(group=group),
                                design = ~group)
}
dds <- DESeq(dds)
res <- results(dds, contrast=c('group', 'positive', 'negative'))
#res <- res[order(res$padj)]
write.table(as.data.frame(res), output_file_DESeq, sep='\t', quote=FALSE, row.names=TRUE)


message('edger_glmlrt')
suppressPackageStartupMessages(library(edgeR))
y <- DGEList(counts=mat, samples=samples, group=group)
y <- calcNormFactors(y, method='TMM')
if(!is.null(batch)){
  # regress out batch information as an additive term
  design <- model.matrix(~group + batch)
} else {
  design <- model.matrix(~group)
}
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
test <- glmLRT(fit, coef=2)
res <- topTags(test, n=nrow(mat), sort.by='none')
mapped_names <- colnames(res)
for(i in 1:ncol(res)){
  if(colnames(res)[i] == 'logFC'){
    mapped_names[i] <- 'log2FoldChange'
  }else if(colnames(res)[i] == 'PValue'){
    mapped_names[i] <- 'pvalue'
  }else if(colnames(res)[i] == 'FDR') {
    mapped_names[i] <- 'padj'
  }else{
    mapped_names[i] <- colnames(res)[i]
  }
}
#colnames(res) <- mapped_names

# write results to file
message('Write results to output file: ', output_file_edgeR)
write.table(res, output_file_edgeR, sep='\t', quote=FALSE, row.names=TRUE)



