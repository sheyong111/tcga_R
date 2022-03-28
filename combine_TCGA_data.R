setwd("H:/毕业设计/毕设生信/TCGA-BRCA/")
#这个函数有三个参数，
#@metadata是从TCGA数据下载的sample sheet
#@path是保存样本表达文件的路径
#@data.type是要合并的数据的类型，这里支持RNAseq和miRNAs两种
merge_TCGA <- function(metadata, path, data.type){
  #通过合并path,还有sample sheet前两列得到每一个文件的完整路径
  filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                         fsep = .Platform$file.sep)
  #判断需要合并的是什么数类型，如果是RNAseq执行下面的代码
  if (data.type=='RNAseq') {
    message ('############### Merging RNAseq data ################\n',
             '### This step may take a few minutes ###\n')
    #通过lapply循环去读每一个样本的表达，然后通过cbind合并成矩阵
    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) 
      read.table(gzfile(fl))$V2))
    #获取第一个文件的第一列作为矩阵的行名
    rownames(rnaMatrix) <- read.table(gzfile(filenames[1]))$V1
    #去掉Ensembl ID后面的.和数字，eg.ENSG00000000003.13
    rownames(rnaMatrix) <- sapply(strsplit(rownames(rnaMatrix), '.', fixed=TRUE), '[',1)
    
    #将sample sheet的sample id这一列作为表达矩阵的列名
    colnames(rnaMatrix) <- metadata$sample_id
    #过滤掉不是基因的行
    index=grepl("^ENSG",rownames(rnaMatrix))
    rnaMatrix <- rnaMatrix[index,]
    
    #统计样本数和基因数
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    
    #输出样本数和基因数
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    #返回最后的基因表达矩阵
    return (rnaMatrix)
    
   }else if (data.type=='miRNAs') { #如果需要合并的是miRNA的数据，执行下面代码
      message ('############### Merging miRNAs data ###############\n')
      #利用lapply来读取每个样本miRNA的表达数据，这里需要用到filtermir这个函数
      #filtermir主要提取mature mir的counts数
      mirMatrix <- lapply(filenames, function(fl) filtermir(fl))
      #mirbase是目前人的说有miRNA成熟提的ID号，eg.MIMAT0027618
      mirs <- mirbase$V1
      #cbind合并成矩阵合并成矩阵，这里注意并不是每个样本中都表达所有的miRNA
      #如果不存在某个miRNA,在这一步表达量会用NA表示
      mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                           function(expr) expr[mirs]))
      
      #设置表达矩阵的行名为miRNA的名字，mirbase的第二列就是所有人的miRNA的名字
      rownames(mirMatrix) <- mirbase$V2
      #设置表达矩阵的列名为sample sheet的sample id这一列
      colnames(mirMatrix) <- metadata$sample_id
      
      #将表达量为NA的地方转换成0
      mirMatrix[is.na(mirMatrix)] <- 0
      
      #统计样本数和miRNA数目
      nSamples = ncol(mirMatrix)
      nGenes = nrow(mirMatrix)
      
      #输出样本数和miRNA的数目
      message (paste('Number of samples: ', nSamples, '\n', sep=''),
               paste('Number of miRNAs: ', nGenes, '\n', sep=''))
      #返回miRNA的表达矩阵
      return (mirMatrix)
    }else{  #如果data.type不是上面提到的两种，就报错，停止执行
      stop('data type error!')
    }
}

#定义filtermir函数
#@fl参数为所有样本miRNA表达文件的链接
filtermir <- function(fl) {
  #readtable读取文件内容
  expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
  #找到miRNA_region这一列以mature开头的行，eg.mature,MIMAT0000062
  expr <- expr[startsWith(expr$miRNA_region, "mature"),]
  #将同一个成熟体的所有counts累加起来，作为这个成熟体counts
  expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
  
  #将miRNA_region这一列（mature,MIMAT0000062）按头号分开，取第二列，得到成熟体ID号
  mirs <- sapply(strsplit(expr$Group.1, ',', fixed=TRUE),'[',2)
  
  #去掉第一列，只保留miRNA成熟体的counts数
  expr <- expr[,-1]
  #加上对应的成熟体ID作为counts的名字
  names(expr) <- mirs
  #返回过滤之后的miRNA成熟体的counts
  return(expr)
}

#读取RNA的sample sheet
  RNA_metadata=read.table("gdc_sample_sheet.2022-03-14 (1).tsv",header=T,sep="\t")
#替换.为下划线，转换成小写
names(RNA_metadata)=tolower(gsub("\\.","_",names(RNA_metadata)))
#调用merge_TCGA函数合并RNA的表达矩阵
rnaCounts=merge_TCGA(metadata=RNA_metadata, path="Raw_data", data.type="RNAseq")
#保存RNA表达矩阵
#write.table(file="result2/combined_RNA_Expr.txt",rnaCounts,sep="\t",quote=F)
write.csv(file="result/combined_RNA_Expr.csv", rnaCounts, quote = F, row.names = T)


#加载mirbase.rds文件，里面保存了人的所有miRNA的成熟体ID和miRNA名字
load("mirbase.rds")
#读取miRNA的sample sheet
mir_metadata=read.table("miRNAs_sample_sheet.tsv",header=T,sep="\t")
#替换.为下划线，转换成小写
names(mir_metadata)=tolower(gsub("\\.","_",names(mir_metadata)))
#调用merge_TCGA函数合并miRNA的表达矩阵
mirCounts=merge_TCGA(metadata=mir_metadata, path="miRNAs_data", data.type="miRNAs")
#保存miRNA表达矩阵
write.table(file="combined_miRNA_Expr.txt",mirCounts,sep="\t",quote=F)


