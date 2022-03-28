library(DESeq2)
library(limma)
library(edgeR)
library(dplyr)
(.packages())#查看已经加载的r包




setwd("H:/毕业设计/毕设生信/TCGA-BRCA")
counts <- read.csv(file = "result/combined_RNA_Expr.csv")
dim(counts)
counts[1:4,1:4]
colnames(counts)[1] <- 'gene_id'#将基因列名设为gene_id，也可在csv文件中手动修改。

#####gene id 转换
gtf <- rtracklayer::import('Homo_sapiens.GRCh38.105.chr.gtf.gz')
gtf <- as.data.frame(gtf)#转化为数据框
#查看文件，保存文件为Rdata，将来方便我们直接打开
dim(gtf)#gtf文件有300多万行，27列
#save(gtf,file = "Homo_sapiens.GRCh38.105基因组注释文件.Rda")
a <- dplyr::filter(gtf,type=="gene",gene_biotype=="protein_coding")#筛选编码蛋白的基因
dim(a)
b <- dplyr::select(a,c(gene_name,gene_id,gene_biotype))#只选择gene_name，gene_id和gene_biotype这三列
b[1:4,]
c <- dplyr::inner_join(b, counts, by = "gene_id")#用counts和b文件中共有的gene_id列还合并文件
c[1:5,1:5]
d=select(c,-gene_id,-gene_biotype)#去掉第二列第三列
mRNAdata=distinct(d,gene_name,.keep_all = T)#基因名去重
mRNAdata[1:4,1:4]

##暂时不用这个
library(org.Hs.eg.db)#基因注释包
mRNAdata1 <- bitr(counts$gene_id,fromType="ENSEMBL",toType="SYMBOL", OrgDb = "org.Hs.eg.db")


#注意gene_name有时候会出现NA。需要去除后在分析。
sum(is.na(mRNAdata$gene_name))
which(is.na(mRNAdata$gene_name))
dim(mRNAdata)
mRNAdata <- na.omit(mRNAdata)#去除NA
dim(mRNAdata)
rownames(mRNAdata)<- mRNAdata[,1]#将第一列作为行名
mRNAdata<-mRNAdata[,-1]#去掉第一列
save(mRNAdata,file = "result/mRNAdata.Rda")
write.csv(mRNAdata,"result/mRNAdata.csv",quote = F,row.names = T)



rm(list = ls())
exprSet <- load(file = "result/mRNAdata.Rda")
exprSet <- mRNAdata
exprSet <- read.table(file = "result/mRNAdata.csv", header = T, row.names = 1, sep = "," )
# exprSet <- read.csv("result2/mRNAdata.csv")
# dim(exprSet)
# exprSet[1:4,1:4]
# row.names(exprSet) <- exprSet[, 1]#第一列设置为行名,与下面一行放在一块
# exprSet <- exprSet[,-1]#去掉第一列多余的
# exprSet[1:4,1:4]
##########表达矩阵筛选
exprSet <- ceiling(exprSet)#?ceiling（）取整，不小于变量本身

#round()以四舍五入形式保存指定小数位数
#floor()不大于变量本身最大整数
#table(is.na(exprSet))
#pmax(),使得矩阵最小值是0？
save <- exprSet
exprSet[exprSet<0] <- 0
table(exprSet<0)
table(is.na(exprSet))
exprSet<-na.omit(exprSet)


##将样品分组
#group_list <- ifelse(as.numeric(str_sub(colnames(exprSet),14,15)<10,"tumor","normal"))
#TCGA可以根据第14和第15位判断是癌组织还是癌旁组织。01表示癌症组织，11表示正常组织
colnames(exprSet)
#提取列名
substr(colnames(exprSet),14,15)
#提取列名字符串的第14和15位置字符
as.numeric(substr(colnames(exprSet),14,15))
#现在提取的字符变成数字 '11'和 11不等同
ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal')
#判断第14和15号位的数值<10 即等于01时候，返回tumor 否则返回normal
factor(ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal'))
#question 把它变成因子型，因子型是有顺序的。这里需要变换成因子型
group_list=factor(ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,'tumor','normal'))
#group_list记录了表达矩阵中患者ID的对应的组织类型，和表达矩阵的顺序一样。
table(group_list)



#######基因差异表达分析DESeq2#################
library(DESeq2)
#构建一个病例号和肿瘤分类的对应关系
colData <- data.frame(row.names = colnames(exprSet),group_list= group_list)
#colData:是患者ID号和组织类型的对应关系

#构建DESeq()函数要求的表达式
dds <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData,design = ~ group_list)
#countData = exprSet,指出DESeq的表达矩阵
#colData = colData,知名每个表达矩阵的分类，比如实验组&对照组，正常组织&癌症组织
#question
#design= ~group_list，因子型，指出不同组的区别，是有顺序的。
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~group_list)
dds <- DESeq(dds)#进行差异基因分析
resultsNames(dds)
#用group_list来做引导文件，用tumor来比较normal组织
res <-  results(dds, contrast=c("group_list","tumor","normal"))
resOrdered <- res[order(res$padj),]#把res差异分析文件通过padj来排序
head(resOrdered)
resOrdered=as.data.frame(resOrdered)#把resOrdered变成数据框，
DEG <- na.omit(resOrdered)


#DEG里添加一列名为change列，标记基因上调下调:
#sd()为数字函数，意为标准差。abs() 也是数字函数，意为绝对值。mean函数是求算数平均值。
#Mean+2*sd可以反应95%以上观测值，比较可信。
#with把所有操作都限制在数据框上(with()的括号内外，信息是完全隔开的，避免程序定位错误的情况)调用数据框DEG里的logFC。
logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)))
logFC_cutoff
#as.factor将函数参数转换为factor类型。
#差异基因筛选，P值小于0.05。pvalue即概率，一般以p<0.05为有统计学差异，对两组差别有显著差异。
#此代码目的对DEG的每一行判定为UP/DOWN/NOT，判定结果添加为change列，并转换为因子格式，得到DOWN、UP、NOT的三种因子。  P.Value<0.05且l logFC l >logFC_cutoff且logFC>logFC_cutoff,输出“UP”    P.Value<0.05且llogFC l >logFC_cutoff且logFC=logFC_cutoff,输出”DOWN”   P.Value>/=0.05或l logFC l =logFC_cutoff，输出“NOT”。
DEG$change = as.factor(
  ifelse(DEG$padj < 0.01 & abs(DEG$log2FoldChange) > logFC_cutoff,
         ifelse(DEG$log2FoldChange >= logFC_cutoff ,'UP','DOWN'),'NOT')
)
head(DEG)
table(DEG$change)
DESeq2_DEG <- DEG
write.csv(DEG, file="result/DESeq_DEG.csv", quote = F, row.names = T)



#######基因差异表达分析edgeR包##################
library(edgeR)
dge <- DGEList(counts = exprSet,group = group_list)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group_list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group_list)
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast=c(-1,1))
DEG2=topTags(fit2, n=nrow(exprSet))
DEG2=as.data.frame(DEG2)
logFC_cutoff2 <- with(DEG2,mean(abs(logFC)) + 2*sd(abs(logFC)))
DEG2$change = as.factor(
  ifelse(DEG2$PValue < 0.05 & abs(DEG2$logFC) > logFC_cutoff2,
         ifelse(DEG2$logFC > logFC_cutoff2 ,'UP','DOWN'),'NOT'))
head(DEG2)
table(DEG2$change)
edgeR_DEG <- DEG2
write.csv(edgeR_DEG, file="result/edgeR_DEG.csv", quote = F, row.names = T)



#######基因差异表达分析limma包###############
library(limma)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design)<- colnames(exprSet)
dge <- DGEList(counts = exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
v <- voom(dge, design, normalize="quantile")
fit <- lmFit(v, design)
constrasts <- paste(rev(levels(group_list)), collapse = "-")
cont.matrix <- makeContrasts(contrasts = constrasts, levels = design)
fit3 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit3)
DEG3 <- topTable(fit3, coef = constrasts, n = Inf)
DEG3 <- na.omit(DEG3)
logFC_cutoff3 <- with(DEG3,mean(abs(logFC)) + 2*sd(abs(logFC)))
logFC_cutoff3
DEG3$change <- as.factor(
  ifelse(DEG3$P.Value < 0.05 & abs(DEG3$logFC) > logFC_cutoff3,
         ifelse(DEG3$logFC > logFC_cutoff3 ,'UP','DOWN'),'NOT')
  )
head(DEG3,n = 10)
table(DEG3$change)
limma_voom_DEG <- DEG3
write.csv(limma_voom_DEG, file="result/limma_voom_DEG.csv", quote = F, row.names = T)

save(DESeq2_DEG,edgeR_DEG,limma_voom_DEG,group_list,file = "result/DEG.Rdata")



###下面的这些也可以封装为函数然后调用见go_kegg.R
###################热图########################################
cg1 <- rownames(DESeq2_DEG)[DESeq2_DEG$change != "NOT"]
cg2 <- rownames(edgeR_DEG)[edgeR_DEG$change != "NOT"]
cg3 <- rownames(limma_voom_DEG)[limma_voom_DEG$change != "NOT" ]

library(pheatmap)
library(RColorBrewer)
#定义热图颜色
color <- colorRampPalette(c('#436eee','white','#EE0000'))(100)
#第一个热图——DESeq2矩阵热图
mat <- exprSet[cg1,]
n <- t(scale(t(mat)))
n[n>1] <- 1
n[n < -1] <- -1
ac <- data.frame(group_list)
rownames(ac) <- colnames(mat)
ht1 <- pheatmap(n,show_colnames = F, show_rownames = F,
                cluster_rows = F, cluster_cols = T,
                annotation_col = ac, color = color)
#画edgeR的热图
mat2 <- exprSet[cg2,]
n2 <- t(scale(t(mat2)))
n2[n2 > 1] <- 1
n2[n2 < -1] <- -1
ht2 <- pheatmap(n2, show_rownames = F, show_colnames = F,
                cluster_rows = F, cluster_cols = T,
                annotation_row = ac, color = color)
#画limma的热图
mat3 <- exprSet[cg3,]
n3 <- t(scale(t(mat3)))
n3[n3 > 1] <- 1
n3[n3 < -1] <- -1
ht3 <- pheatmap(n3, show_rownames = F, show_colnames = F,
                cluster_rows = F, cluster_cols = T,
                annotation_col = ac, color = color)




###############火山图##########################

library(EnhancedVolcano)
library(airway)

EnhancedVolcano(DESeq2_DEG,
                lab = rownames(DESeq2_DEG),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-8,8),
                title = 'DESeq火山图',
                pCutoff = 10e-17,
                FCcutoff = 2.5,
                )


#############三大R包差异基因对比##########

########GO富集和KEGG富集
diff <- read.table(file = "result/DESeq2备份.csv",sep = ",")
diff[1:4,1:4]
colnames(diff[1]) <- 'gene_id'





