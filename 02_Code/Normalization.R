# 0. Package&Function ----
library(readxl)
library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(openxlsx)
library(mixOmics)
library(metaX)
source("./02_Code/QC_PCA.R") 
source("./02_Code/QC_boxplot.R")
source("./02_Code/QC_heatmap.R")
source("./02_Code/run_DE.R")
source("./02_Code/run_enrichment_analysis.R")

# 1. Data input ----
## 1.1 Group input ----
# 导入分组信息
data_group <- read_excel("./01_Data/sam_qc_infor_neg.xlsx")
data_group$id <- gsub("neg_","", data_group$id)            # 去除样本id多余信息
data_group$id <- gsub("cas_","", data_group$id)
data_group <- as.data.frame(data_group)

# 提取分组样本
data_group_OCI_VEN <- data_group[c(9:17),]
data_group_OCI_VEN[grep("WT",data_group_OCI_VEN$id),2] <- "WT"
data_group_molm13_VEN[-grep("WT",data_group_molm13_VEN$id),2] <- "6W"

# 配色设置 
value_colour <- c(#"6W" = "#E64B35FF",
                  "2W" = "#F2A200",
                  "WT" = "#4DBBD5FA")
                  #"OCI" = "#F2A200",
                  #"QC" = "#8D8D8D"
rownames(data_group) <- data_group$id



## 1.2 meta matrix input ----
data_input <- read.xlsx("./01_Data/meta_intensity_class_neg.xlsx")
data_input <- as.data.frame(data_input)
rownames(data_input) <- data_input$Name
data_input <- data_input[,-grep("4w",colnames(data_input))]   # 删除4w样本
colnames(data_input) <- gsub("neg_","", colnames(data_input)) # 去除样本id多余信息
colnames(data_input) <- gsub("cas_","", colnames(data_input))

# 保留注释
data_anno <- data_input[,1:13]
data_anno <- as.data.frame(data_anno)
rownames(data_anno) <- data_anno$Name
data_input <- data_input[,-1:-13]
write.xlsx(data_anno,file = "./01_Data/data_anno_pos.xlsx")


# 2. Normalization -----------------------------------------------------------
# 设置输出目录
dir.create("./03_Result/QC/neg/OCI/6W_WT")
dir <- "./03_Result/QC/neg/OCI/2W_WT/"

## 2.1 Intensity normalization ----
log_data <- log2(data_input)

# 计算校正前各样本的intensity median
column_medians <- apply(log_data, 2, median, na.rm = TRUE)
column_medians

# 目标中位数
target_median <- max(column_medians)
data_after <- log_data
for (col in names(data_after)) {
  if (is.numeric(data_after[[col]])) {
    # 防止除以零的情况
    if (column_medians[col] != 0) {
      data_after[[col]] <- data_after[[col]] / column_medians[col] * target_median
    }
  }
}

# 验证校正后各样本intensity median是否一致
column_medians_2 <- apply(data_after, 2, median, na.rm = TRUE)
column_medians_2

# 返回log2之前的数据
data_input_norm <- 2 ^ data_after
write.xlsx(data_input_norm,file = "./01_Data/meta_intensity_norm_class_pos.xlsx")

# 3. QC --------------------------------------------------------------------------
# QC函数data_group项需要满足包含id,group列！！！！！！！！！！！！！！！！！！！
colnames(data_group)
data_group <- data_group[,-(3:4)]
## 3.1 Boxplot -----------------------------------------------------------------
# 函数进行了log2
pdf(file = paste0(dir,"QC_boxplot_normalization.pdf"),
    width = 6,
    height = 4)
QC_boxplot(data_input,data_group = data_group_OCI_2w,
           value_colour = value_colour,
           title = "normalized data")
dev.off()

pdf(file = paste0(dir,"QC_boxplot_normalization.pdf"),
    width = 6,
    height = 4)
QC_boxplot(data_input_norm,data_group = data_group,
           value_colour = value_colour,
           title = "normalized data")
dev.off()

## 3.2 Heatmap -----------------------------------------------------------------
# 函数中有log2（x+1）
pdf(file = paste0(dir,"QC_heatmap_normalization.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_input,data_group = data_group_OCI_2w,
           value_colour = value_colour)
dev.off()

pdf(file = paste0(dir,"QC_heatmap_normalization.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_input_norm,data_group = data_group,
           value_colour = value_colour)
dev.off()

## 3.3 PCA ---------------------------------------------------------------------
pdf(file = paste0(dir,"QC_pca_normalization.pdf"),
    width = 7,
    height = 7)
QC_PCA(data = log2(data_input+1),
       data_group = data_group_OCI_2w,
       value_colour = value_colour)
dev.off()

pdf(file = paste0(dir,"QC_pca_normalization.pdf"),
    width = 7,
    height = 7)
QC_PCA(data = data_input_norm,
       data_group = data_group,
       value_colour = value_colour)
dev.off()

## 3.4 RSD ---------------------------------------------------------------------
# QC样本中有超过70%的潜在特征峰（化合物）的相对标准偏差（RSD）不超过30%,说明检测体系的稳定性良好
if(T){
# 提取 QC 样本
data_qc <- data_input_norm[,grep("QC",colnames(data_input_norm))]

# 计算 QC 样本的相对标准偏差 (RSD)
rsd <- apply(data_qc, 1, function(x) sd(x) / mean(x) * 100)
qc_metrics <- data.frame(Metabolite = rownames(data_qc), RSD = rsd)

# 判断稳定性是否良好：超过 70% 的潜在特征峰 RSD <= 30%
valid_ratio <- mean(rsd <= 30)

# 输出稳定性评估结果
print(ifelse(valid_ratio > 0.7, "检测体系稳定性良好", "检测体系稳定性较差"))
}

# 对所有代谢物进行Shapiro-Wilk检验
results_1 <- apply(data_input[,1:8], 1, function(x) shapiro.test(x)$p.value)
table(results_1<0.05)

# 4. DE ----------------------------------------------------------------------
# 注意，data和data_anno的行名应一致
# 根据分组选择要进行差异分析的组别
source("./02_Code/run_DE.R")
table(data_group$group)
targeted_group <- data_group_p53

# VEN-WT 组别
targeted_group$group <- paste0(targeted_group$group, "_VEN")
targeted_group[grep("WT",targeted_group$id),2] <- gsub("_VEN","_WT",
                                                       targeted_group[grep("WT",targeted_group$id),2])
table(targeted_group$group)

# 2W/6W-WT 组别
targeted_group[grep("6w",targeted_group$id),2] <- gsub("_VEN","_6W",
                                                       targeted_group[grep("6w",targeted_group$id),2])
table(targeted_group$group)


#dir.create("./03_Result/DE/OCI_AML2")
## 4.1 Set group ----
group_1 <- "TP53"          # 实验组
group_2 <- "WT"          # 对照组
dir_DE <- "./03_Result/DE/pos/TP53_vs_WT/"

# 若选择wilcoxon检验，检查是否有平局值 
anyDuplicated(data_input_norm)    # 结果大于0代表有

result_merge <- run_DE(data = data_input,
                       data_group = targeted_group,
                       data_anno = data_anno,
                       group_1 = group_1,group_2 = group_2,
                       log2 = TRUE,
                       logfc_threshold = 0.263,        # log2fc值,1.2倍fc
                       pvalue_threshold = 0.05, 
                       qvalue_threshold = NULL,
                       test_method = "t-test",     # "t-test" or "wilcoxon"
                       paired = FALSE ,            # 是否配对检验，TRUE or FALSE 必须为逻辑值
                       dir = "03_result/DE/pos/") # 每次需要更改
# 统计上下调Meta个数
table(result_merge$change)

## 4.2 PLS-DA ----
# 载入metaX package
# 定量数据第一列列名必须为“name”
# 分组数据QC样本分组必须设置为"NA"
para <- new("metaXpara")
pfile <- "./01_Data/meta_intensity_class_neg.txt"
sfile <- "./01_Data/sam_qc_infor_neg.txt"
rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
sampleListFile(para) <- sfile
para <- reSetPeaksData(para)
#para <- missingValueImpute(para)
para <- metaX::normalize(para,method="pqn")
plsdaPara <- new("plsDAPara")
# 选择最优主成分
best_comp <- selectBestComponent(
  para = para,               # 输入的metaXpara对象
  np = 10,                     # 最大主成分数为10
  sample = c("MOLM13_6W", 
             "MOLM13_WT"),       # 选择的样本组
  scale = "uv",           # 数据标准化方法：pareto
  valueID = "value",          # 代谢物定量数据列名称
  k = 7                       # 7折交叉验证
)
plsdaPara@nperm <- 200

plsda.res <- runPLSDA(para = para,plsdaPara = plsdaPara,
                      sample = c("MOLM13_6W","MOLM13_WT"),valueID="value")
cat("R2Y: ",plsda.res$plsda$res$R2,"\n")
## R2Y: 0.9840716

cat("Q2Y: ",plsda.res$plsda$res$Q2,"\n")
## Q2Y: 0.9806666

## permutation test

cat("P-value R2Y: ",plsda.res$pvalue$R2,"\n")
## P-value R2Y: 0

cat("P-value Q2Y: ",plsda.res$pvalue$Q2,"\n")
## P-value Q2Y: 0

# PLS-DA Loading Plot 

## 4.3 Volc Plot ----
# change列因子化
result_merge$change <- factor(
  result_merge$change,
  levels = c("up", "down", "stable"))  # 强制按此顺序排列

# 计算对称的横坐标范围（基于logFC绝对值的最大值）
max_abs_logfc <- max(abs(result_merge$logFC), na.rm = TRUE)

# 扩展5%的余量，避免点紧贴坐标轴边缘
x_limit <- max_abs_logfc * 1.05

logfc_threshold <- 0.263
pvalue_threshold <- 0.05

p <- ggplot(data = result_merge, aes(x =logFC, y = -log10(pvalue))) +
  geom_point(alpha = 0.5, aes(color = change, size = VIP )) +
  ylab("-log10(Pvalue)")+
  # 按因子顺序指定颜色
  scale_color_manual(
    name = "Change",                              # 图例标题
    values = c(
      "up" = "#B30000",      # 红
      "stable" = "grey",     # 灰
      "down" = "#003366"),    # 蓝,
    labels = c(                                  # 显示上下调的个数
      "up" = paste0("Up ：", sum(result_merge$change == "up")),
      "stable" = paste0("Stable ：", sum(result_merge$change == "stable")),
      "down" = paste0("Down ：", sum(result_merge$change == "down"))))+
  scale_size_continuous(name = "VIP",        # 为 size 变量添加图例标题
                        range = c(1, 6))+    # 设置点的大小范围，可以根据需要调整
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), lty = 4, 
             col = "black", lwd = 0.8,alpha=0.4) +
  geom_hline(yintercept = -log10(pvalue_threshold), lty = 4, 
             col = "black", lwd = 0.8,alpha=0.4) +
  labs(title = paste0(group_1,"-",group_2)) +
  xlim(-x_limit, x_limit)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        aspect.ratio = 1.2)  # 设置纵横比，调整为更高

ggsave(filename = paste0(dir_DE,"/volc.pdf"),
       plot = p, device = "pdf", 
       width = 6, height = 5)
}

# 5. KEGG GO -----------------------------------------------------------------
GeneSymbol <- subset(result_merge,adj.P.Val< 0.05)
# GeneSymbol$Genes <- GeneSymbol$GN
y <- GeneSymbol$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))

# 样本为人，OrgDb为Hs
# 样本为小鼠，OrgDb为Mm
run_enrichment_analysis(data = GeneSymbol,
                        OrgDb = "Hs",
                        dir = paste0("03_result/DE/",
                                     group_1,"_vs_",group_2,"/"))

