#作图函数是在网上找的，在3-plotfunction.R里
rm(list = ls())
load("result/mRNAdata.Rda")
expr <- mRNAdata
load("result/DEG.Rdata")
source("R_script/3-plotfunction.R")
logFC_cutoff <- 1
expr[1:4,1:4]


dat = log(expr+1)
pca.plot = draw_pca(dat,group_list)

cg1 = rownames(DESeq2_DEG)[DESeq2_DEG$change !="NOT"]
cg2 = rownames(edgeR_DEG)[edgeR_DEG$change !="NOT"]
cg3 = rownames(limma_voom_DEG)[limma_voom_DEG$change !="NOT"]

h1 = draw_heatmap(expr[cg1,],group_list)
h2 = draw_heatmap(expr[cg2,],group_list)
h3 = draw_heatmap(expr[cg3,],group_list)
v1 = draw_volcano(test = DESeq2_DEG[,c(2,5,7)],pkg = 1)
v2 = draw_volcano(test = edgeR_DEG[,c(1,4,6)],pkg = 2)
v3 = draw_volcano(test = limma_voom_DEG[,c(1,4,7)],pkg = 3)

library(patchwork)
(h1 + h2 + h3) / (v1 + v2 + v3) +plot_layout(guides = 'collect')



# 三大R包差异基因交集
UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}

up = intersect(intersect(UP(DESeq2_DEG),UP(edgeR_DEG)),UP(limma_voom_DEG))
down = intersect(intersect(DOWN(DESeq2_DEG),DOWN(edgeR_DEG)),DOWN(limma_voom_DEG))

hp = draw_heatmap(expr[c(up,down),],group_list)

#上调、下调基因分别画维恩图

up.plot <- venn(UP(DESeq2_DEG),UP(edgeR_DEG),UP(limma_voom_DEG),
                "UPgene"
)
down.plot <- venn(DOWN(DESeq2_DEG),DOWN(edgeR_DEG),DOWN(limma_voom_DEG),
                  "DOWNgene"
)
#维恩图拼图，

library(cowplot)
library(ggplotify)
up.plot = as.ggplot(as_grob(up.plot))
down.plot = as.ggplot(as_grob(down.plot))
library(patchwork)
#up.plot + down.plot
pca.plot + hp+up.plot +down.plot
ggsave("result/deg.png",height = 10,width = 10)




rm(list = ls())
diff <- read.csv(file = "result/DESeq_DEG.csv")
diff[1:4,1:4]
colnames(diff)[1] <- 'geneSymbol'
#1、加载包
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包
#2、基因id转换，kegg和go富集用的ID类型是ENTREZID）
gene.df <- bitr(diff$geneSymbol,fromType="SYMBOL",toType="ENTREZID", OrgDb = "org.Hs.eg.db")
gene <- gene.df$ENTREZID
#3、GO富集
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P值可以取0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)



#4有时候我们富集到的通路太多了，不方便绘制图片，可以选择部分绘制，
display_number = c(22, 22, 22)#这三个数字分别代表选取的BP、CC、MF的通路条数，这个自己设置就行了
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}

##开始绘制GO柱状图
###横着的柱状图
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

###竖着的柱状图 
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  theme_bw() + 
  xlab("GO term") + 
  ylab("Num of Genes") + 
  labs(title = "The Most Enriched GO Terms")+ 
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))#angle是坐标轴字体倾斜的角度，可以自己设置


#####################KEGG富集和可视化###########################
#1、KEGG富集
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)

#2、可视化
###柱状图
hh <- as.data.frame(kk)#自己记得保存结果哈！
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "red",high ="blue" )+#颜色自己可以换
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

###气泡图
hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()
ggsave(file="result/kegg气泡图.png" , ggplot(hh,aes(y=order,x=Count))+
         geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
         scale_color_gradient(low="green",high = "red")+
         labs(color=expression(p.adjust,size="Count"), 
              x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
         theme_bw())

