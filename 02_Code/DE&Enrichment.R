# 0. Package&Function ----------------------------------------------------------
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)
library(openxlsx)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(devtools)
remotes::install_github("wenbostar/metaX")  # github
BiocManager::install("wenbostar/metaX")
source("./02_Code/QC_PCA.R")
source("./02_Code/QC_boxplot.R")
source("./02_Code/QC_heatmap.R")
source("./02_Code/run_DE.R")
source("./02_Code/run_enrichment_analysis.R")


# 1. Data input ----------------------------------------------------------------
## 1.1 Group input ----
# 导入分组信息
data_group <- read_excel("./01_Data/MOLM13_group.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
data_group$group<- gsub("_WT","",data_group$group)
data_group[grep("WT",data_group$id),2] <- "MOLM13_WT"
#table(data_group$group)

# 配色设置
# 配色设置
value_colour <- c("MOLM13_WT" = "#E64B35FF",
                  "MOLM13_2W" = "#4DBBD5FA",
                  "MOLM13_6W" = "#F2A200")
rownames(data_group) <- data_group$id

## 1.2 DIA matrix input ----
fill_norm <- read.csv("./01_Data/MOLM13_pg_matrix_fill_norm.csv")
colnames(fill_norm)
rownames(fill_norm) <- fill_norm$X
fill_norm <- fill_norm[,-1]
#MV4_11_fill_norm <- fill_norm[,grep("MV4_11",colnames(fill_norm))]
#write.csv(MV4_11_fill_norm,file = "./01_Data/MV4_11_pg_matrix_fill_norm.csv")

# remove abnormal sample
data_group <- data_group[-grep("4W",data_group$id),]
fill_norm <- fill_norm[,-grep("4W",colnames(fill_norm))]

# Protein anno input
data_anno <- read.xlsx("./01_Data/data_anno.xlsx")
data_anno <- as.data.frame(data_anno)
rownames(data_anno) <- data_anno$Protein.Group

# 2. Set output category -------------------------------------------------------------
#dir.create("./03_Result/GO&KEGG/MOLM13")
dir <- "./03_Result/QC/MV4_11/"

# 3. QC ------------------------------------------------------------------------
## 3.1 Boxplot -----------------------------------------------------------------
pdf(file = paste0(dir,"QC_boxplot_normalization.pdf"),
    width = 6,
    height = 4)
QC_boxplot(fill_norm,data_group = data_group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()
## 3.2 Heatmap -----------------------------------------------------------------
pdf(file = paste0(dir,"QC_heatmap_normalization.pdf"),
    width = 6,
    height = 6)
QC_heatmap(fill_norm,data_group = data_group,
           value_colour = value_colour)
dev.off()

## 3.3 PCA ---------------------------------------------------------------------
pdf(file = paste0(dir,"QC_pca_normalization.pdf"),
    width = 7,
    height = 7)
QC_PCA(data = fill_norm,
       data_group = data_group,
       value_colour = value_colour)
dev.off()

# 4. DE ------------------------------------------------------------------------
# 注意，data和data_anno的行名应一致
data_anno <- data_anno[rownames(data_anno)%in%rownames(oci_fill_norm),]
data_anno <- read_xlsx("./01_Data/data_anno.xlsx")
colnames(data_anno)
data_anno <- as.data.frame(data_anno)
rownames(data_anno) <- data_anno$Protein.Group
## 4.1 Set output catagory ----
dir <- "./03_Result/DE/OCI_AML2/OCI_VEN_vs_OCI_WT/"

## 4.2 Set group ----
table(data_group$group)
data_group[grep("6W",data_group$id),2] <- "MV4_11_VEN"

# 根据分组选择要进行差异分析的组别
# group 1为实验组
# group 2为对照组
group_1 <- "MOLM13_2W"
group_2 <- "MOLM13_WT"
result_merge <- run_DE(data = fill_norm,
                       data_group = data_group,
                       log2 = T,
                       data_anno = data_anno,
                       group_1 = group_1,group_2 = group_2,
                       dir = "./03_Result/DE/MOLM13/")

## 4.3 Volcano plot ------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

#dif_logFC_up <- subset(result_merge,result_merge$logFC > 0)
#dif_logFC_up <- dif_logFC_up[order(dif_logFC_up$logFC,decreasing = T),]
#dif_logFC_down <- subset(result_merge,result_merge$logFC < 0)
#dif_logFC_down <- dif_logFC_down[order(dif_logFC_down$logFC),]

# y <- c(rownames(dif_logFC_up)[1:10],rownames(dif_logFC_down)[1:10])
# y <- na.omit(y)


## **DE_res input（Done） ----
result_merge <- read.csv("./03_Result/DE/OCI_AML2/OCI_VEN_vs_OCI_WT/result_DE.csv")
result_merge <- result_merge[,-1]

if(T){
res_data <- result_merge
data <- res_data[res_data[,"P.Value"] <= 1,]

y <- result_merge$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
data$gene <- gene

# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3
logfc <- 0.585                                 # 1.5倍变化
pvalue <- 0.05

data$sig[data$P.Value >= pvalue | abs(data$logFC) < logfc] <- "Not"
data$sig[data$P.Value < pvalue & data$logFC >= logfc] <- "Up"
data$sig[data$P.Value < pvalue & data$logFC <= -logfc] <- "Down"

# 火山图 
p <- ggplot(data = data, aes(x =logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.7, size = 3.5, 
             aes(color = sig)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("#4DBBD5", "grey80", "#E64B35"))+
  geom_vline(xintercept = c(-logfc, logfc), lty = 4, 
             col = "black", lwd = 0.8,alpha=0.4) +
  geom_hline(yintercept = -log10(pvalue), lty = 4, 
             col = "black", lwd = 0.8,alpha=0.4) +
  labs(title = paste0(group_1,"-",group_2),) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中

ggsave(filename = paste0(dir,"/volc.pdf"),
       plot = p, device = "pdf", 
       width = 6, height = 5)
}

# 统计上下调基因数量
table(data$sig)

## 4.4 Annotated volcano plot ----
library(readxl)
list <- read_xlsx("01_rawdata/A list of positive control proteins.xlsx")
y <- list$`UniProt accession`
anno_gene_symbol <- unlist(lapply(y,function(y) strsplit(as.character(y),"_")[[1]][1]))
Gene <- as.data.frame(anno_gene_symbol)
library(readr)
rownames(result_merge) <- result_merge$Protein.Names

volc +
  geom_text_repel(data = data[data$Protein.Group %in% Gene[,1],],
                  aes(label=gene),
                  size = 3,
                  color = "black",
                  max.overlaps = 50) +
  geom_point(data = data[data$Protein.Group %in% Gene[,1],],
             alpha = 0.9,color = "black")
ggsave(plot = last_plot(),filename = paste0(dir,"volc.pdf"),
       height = 5,
       width = 6)

# 5. KEGG GO -------------------------------------------------------------------
## 5.1 Set output catagory----
#OCI_AML2/MV4_11/MOLM13
# 指定文件夹路径
dir.create("./03_Result/GO&KEGG/MOLM13/VEN_vs_WT/")
dir <- "./03_Result/GO&KEGG/OCI_AML2/VEN_vs_WT/"

## 5.2 DE_res input ----
DE_result <- read.csv('./')

## 5.3 set P.Value ----
GeneSymbol <- subset(data, P.Value < 0.05)

## 5.4 set cutoff值 ----
cutoff <- 0.585                  # 1.5倍变化

# 转换基因名 
#y <- result_merge$Genes
#gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
#data$gene <- gene

if(T){
## 5.5 down genes ----
down_genes <- subset(GeneSymbol, logFC < -cutoff)

# 设置数据库 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa'         # KEGG是hsa数据库

# gene ID转换 
gene <- bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

### 5.5.1 GO ----
# GO富集分析
kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                keyType = "ENTREZID", # 设定读取的gene ID类型
                ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)# 设定q值阈值

### 5.5.2 KEGG ----
# KEGG富集分析
kk <- enrichKEGG(gene = gene$ENTREZID,
                 keyType = "kegg",
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

## GO、KEGG结果整合 
result <- list(enrichGO = kkd, enrichKEGG = kk)

# 结果标记为下调 
result_down <- result
kkd_down <- result_down$enrichGO
kk_down <- result_down$enrichKEGG

### 5.5.3 down res_output ----
# 导出下调enrichGO 
write.csv(kkd_down@result, file = paste0(dir, "/GO_down.csv"), quote = F, row.names = F)

# dotplot
pdf(file = paste0(dir, "/GO_down.pdf"), width = 6, height = 7)
p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
print(p1)
dev.off()

# 导出下调enrichKEGG
write.csv(kk_down@result, file = paste0(dir, "/KEGG_down.csv"), quote = F, row.names = F)

# dotplot
pdf(file = paste0(dir, "/KEGG_down.pdf"), width = 6, height = 5)
p2 <- dotplot(kk_down,showCategory = 10)
print(p2)
dev.off()
}

if(T){
## 5.6 up genes ----
up_genes <- subset(GeneSymbol, logFC > cutoff)

### 5.6.1 GO-up ----
# gene ID转换 
gene <- bitr(up_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

# GO富集分析
kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                keyType = "ENTREZID", # 设定读取的gene ID类型
                ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2, # 设定q值阈值
                readable = T)
### 5.6.2 KEGG-up ----

# KEGG富集分析
kk <- enrichKEGG(gene = gene$ENTREZID,
                 keyType = "kegg",
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

# GO、KEGG结果整合
result <- list(enrichGO = kkd, enrichKEGG = kk)

# 结果标记为上调
result_up <- result
kkd_up <- result_up$enrichGO
kk_up <- result_up$enrichKEGG

### 5.6.3 up res_output ----

# 导出上调enrichGO
write.csv(kkd_up@result, file = paste0(dir, "/GO_up.csv"), quote = F, row.names = F)

# dotplot
pdf(file = paste0(dir, "/GO_up.pdf"), width = 6, height = 7)
p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p3)
dev.off()

# 导出上调enrichKEGG
write.csv(kk_up@result, file = paste0(dir, "/KEGG_up.csv"), quote = F, row.names = F)

# dotplot
pdf(file = paste0(dir, "/KEGG_up.pdf"), width = 6, height = 5)
p4 <- dotplot(kk_up,showCategory = 10)
print(p4)
dev.off()
}

kegg <- kk_up@result[kk_up@result$pvalue<0.05,]
p2 <- ggplot(kegg,aes(x=GeneRatio,y=Description))+
      geom_point(aes(size=Count,color= -log10(pvalue)))+
      theme_bw()+labs(y="",x="GeneRatio")+ 
      scale_color_gradient(low="blue",high="red")+
  #scale_size_continuous(range = c(3, 12))+  # 调整气泡的大小范围
      theme(axis.text.y = element_text(angle = 0, hjust = 1))  # 调整Y轴标签角度

ggsave(plot = p2,filename = "./03_Result/GO&KEGG/MOLM13/6W_vs_WT/downkegg.pdf")

## 5.7 统计上下调的通路的数量 ----
# 下调GO 
table(kkd_down@result$p.adjust<0.05)
# 下调KEGG 
table(kk_down@result$p.adjust<0.05)
# 上调GO 
table(kkd_up@result$p.adjust<0.05)
# 上调KEGG 
table(kk_up@result$p.adjust<0.05)

# 6.heatmap------
library(reshape2)
library(ComplexHeatmap)
library(openxlsx)
library(gtable)
library(grid)
## df structure
# gene    module
# PCDH9   tan
# ALX1    tan
# CGREF1  darkorange
# 
df <- data.frame(gene = colnames(cleanedData), module = mergedColors)
df$gene <- gsub(".*\\|", "", df$gene)
kegg_list <- list()
i <- 0
all_paths <- c()
for (module in unique(df$module)) {
  i <- i + 1
  genes <- df[df$module == module,]$gene
  trans_id <- bitr(genes, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
  kegg <- enrichKEGG(gene = trans_id$ENTREZID, organism = "hsa", 
                     keyType = "kegg", 
                     qvalueCutoff = 0.05)
  result <- kegg@result
  result$module <- module
  kegg_list[[i]] <- result
  significant_paths <- result[result$pvalue < 0.05, "Description"]
  all_paths <- c(all_paths, significant_paths)
  #all_paths <- c(all_paths, kegg@result$Description[c(1,2)])
  #all_paths <- c(all_paths, kegg@result$Description[1])
}

paths <- c("PI3K-Akt signaling pathway", "MAPK signaling pathway")
all_paths <- unique(c(all_paths, paths))
kegg_all <- do.call(rbind, kegg_list)
kegg_all <- kegg_all[kegg_all$Description %in% all_paths,]
plot_data <- dcast(kegg_all[,c("Description", "pvalue", "module")], Description ~ module, value.var = "pvalue", fill = 1)

write.xlsx(plot_data, file = "03_result/WGCNA/KEGGplotdata.xlsx",sep=",", quote=F)
plot_data <- read_excel("03_result/WGCNA/KEGGplotdataselected.xlsx")
plot_data<- as.data.frame(plot_data)
rownames(plot_data) <- plot_data$Description
plot_data$Description <- NULL
plot_data <- -log(plot_data)
plot_data[plot_data>5] <- 5
#colnames(plot_data) <- paste0("Cluster", seq(ncol(plot_data)))
plot_data <- as.matrix(plot_data)
pdf(file =  "03_result/WGCNA/10.pathwayheatmap.pdf", width = 10, height = 10)
p<- pheatmap(plot_data,name =" -log(pvalue)",
             cluster_cols = T,
             cluster_rows = F,
             color = colorRampPalette(c("white", "blue"))(100),
             
