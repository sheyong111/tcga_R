rm(list = ls())#一键清空
#安装并加载R包
if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("rtracklayer")) BiocManager::install("rtracklayer")
if(!require("tidyr")) BiocManager::install("tidyr")
if(!require("dplyr")) BiocManager::install("dplyr")
if(!require("DESeq2")) BiocManager::install("DESeq2")
if(!require("limma")) BiocManager::install("limma")
if(!require("pheatmap")) BiocManager::install("pheatmap")
if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
if(!require("clusterProfiler")) BiocManager::install("clusterProfiler")

library(limma)
library(edgeR)
library(dplyr)
a <- read.delim2("combined_RNA_Expr.txt", row.names = 1)
dim(a)
a[0:4,1:4]










