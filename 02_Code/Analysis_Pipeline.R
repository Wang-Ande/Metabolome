# 0. Package&Function ----
pak::pak("xia-lab/MetaboAnalystR")
library(readxl)
library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(openxlsx)
library(pmp)
source("./02_Code/QC_PCA.R") 
source("./02_Code/QC_boxplot.R")
source("./02_Code/QC_heatmap.R")
source("./02_Code/run_DE.R")
# source("./02_Code/run_enrichment_analysis.R")

# 1. Data input ----
## 1.1 Group input ----
# 导入分组信息
data_group <- read.xlsx("./01_Data/01.MetQuant/sam_infor_combined.xlsx")
data_group <- as.data.frame(data_group)

## 1.2 meta matrix input ----
data_input <- read.csv("./01_Data/01.MetQuant/meta_intensity_combined.csv",row.names = 1)
data_input <- as.data.frame(data_input)
colnames(data_input) <- gsub("neg_","", colnames(data_input))   # 去除样本id多余信息
colnames(data_input) <- gsub("cas_","", colnames(data_input))

# 保留注释
data_anno <- data_input[,1:12]
data_anno <- as.data.frame(data_anno)
write.xlsx(data_anno,file = "./01_Data/01.MetQuant/data_anno_pos.xlsx")
data_anno <- read.xlsx("./01_Data/01.MetQuant/meta_anno_combined.xlsx",rowNames = TRUE)

# 2. Normalization -----------------------------------------------------------
## 2.1 Intensity normalization ----
# method 1
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
write.csv(data_input_norm,file = "./01_Data/01_MetQuant/meta_intensity_median_combined.csv")

# method 2
# 使用 PQN 进行归一化，指定 QC 样本为参考样本
class <- colnames(data_input)
colnames(data_input)[30] <- "QC"
normalized_data <- pmp::pqn_normalisation(df = data_input,
                                          classes = class,
                                          qc_label = "QC",
                                          ref_method = "mean")
# methods 3
# 计算 QC 样本的中位数
qc_index <- c(28,29,30)
qc_medians <- apply(data_input[, qc_index], 1, median, na.rm = TRUE)

# 进行 QC-based Median 矫正
normalized_data <- sweep(data_input, 1, qc_medians, "/")

# 3. QC --------------------------------------------------------------------------
# 设置输出目录
dir.create("./03_Result/QC/combined/OCI_M2")
dir_qc <- "./03_Result/1.QC/combined/OCI_M2/"

# QC函数data_group项需要满足包含id,group列！！！！！！！！！！！！！！！！！！！
data_qc <- data_filter
colnames(data_group)
group_qc <- data_group[,c(1,3)]
colnames(group_qc)[2] <- "group"
group_qc <- group_qc[grep("OCI",group_qc$id),]

# 配色设置 
value_colour <- c("High" = "#E64B35FF",
                  "Low" = "#F2A200",
                  "Con" = "#4DBBD5FA")
#"QC" = "#8D8D8D"
rownames(data_group) <- data_group$id

## 3.1 Boxplot -----------------------------------------------------------------
# 函数进行了log2
pdf(file = paste0(dir_qc,"QC_boxplot_filter.pdf"),
    width = 6,
    height = 4)
QC_boxplot(data_qc,data_group = group_qc,
           value_colour = value_colour,
           title = "normalized data")
dev.off()

## 3.2 Heatmap -----------------------------------------------------------------
# 函数中有log2（x+1）
pdf(file = paste0(dir_qc,"QC_heatmap_filter.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_qc,data_group = group_qc,
           value_colour = value_colour)
dev.off()

## 3.3 PCA ---------------------------------------------------------------------
pdf(file = paste0(dir_qc,"QC_pca_filter.pdf"),
    width = 7,
    height = 7)
QC_PCA(data = log2(data_qc+1),
       data_group = group_qc,
       value_colour = value_colour)
dev.off()

## 3.4 RSD ---------------------------------------------------------------------
# QC样本中有超过70%的潜在特征峰（化合物）的相对标准偏差（RSD）不超过30%,说明检测体系的稳定性良好
if(T){
# 提取 QC 样本
qc_sam <- data_input[,grep("QC",colnames(data_input))]

# 计算 QC 样本的相对标准偏差 (RSD)
rsd <- apply(qc_sam, 1, function(x) sd(x) / mean(x) * 100)
qc_metrics <- data.frame(Metabolite = rownames(qc_sam), RSD = rsd)

# 判断稳定性是否良好：超过 70% 的潜在特征峰 RSD <= 30%
valid_ratio <- mean(rsd <= 30)

# 输出稳定性评估结果
print(ifelse(valid_ratio > 0.7, "检测体系稳定性良好", "检测体系稳定性较差"))
}

# 去除QC样本中RSD超过30%的样本
Good_metabo <- qc_metrics[qc_metrics$RSD <= 30,]
Bad_metaba <- data_anno[qc_metrics$RSD > 30,]
View(Bad_metaba$Name)
data_filter <- data_input[rownames(Good_metabo),]

## 3.5 Shapiro-Wilk ----
# 对所有代谢物进行Shapiro-Wilk检验(判断是否满足正态性)
results_1 <- apply(data_input[,1:8], 1, function(x) shapiro.test(x)$p.value)
table(results_1<0.05)

# 4. DE ----------------------------------------------------------------------
# 注意，data和data_anno的行名应一致
# 根据分组选择要进行差异分析的组别
source("./02_Code/run_DE.R")
table(data_group$group)
targeted_group <- data_group[grep("OCI_M2",data_group$id),]
targeted_group <- targeted_group[,c(1,3)]
# colnames(targeted_group)[2] <- "group"

## 4.1 Set group ---------------------------------------------------------------
group_1 <- "High"        # treatment
group_2 <- "Con"        # control

# 若选择wilcoxon检验，检查是否有平局值 
anyDuplicated(data_input)    # 结果大于0代表有

## 4.1 LogFC & P-value ---------------------------------------------------------
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
                       dir = "03_result/2.DE/combined/OCI_M2/") # 每次需要更改
# 统计上下调Meta个数
table(result_merge$change)

## 4.2 PLS-DA ------------------------------------------------------------------
# Input Normalization Data
# Data format:samples in rows and variables in columns.
# X: gene expression matrix 
# Y: factor indicating sample class membership

### 4.2.1 load data ----
X <- read.csv("./03_Result/2.DE/combined/OCI_M2/High_vs_Con/DE_results.csv",row.names = 1)
X <- X[,grep("_WT|_6W",colnames(X))]
X <- t(X)
X <- log2(X)

Y <- targeted_group[grep("Con|High",targeted_group$group),]
Y <- factor(Y$group,levels = c("Con","High"))

# check
dim(X); length(Y)
summary(Y)

### 4.2.2 method 1 mixOmics ----------------------------------------------------
library(mixOmics)
# Initial exploration with PCA 
pca <- pca(X, ncomp = 3, scale = TRUE)
plotIndiv(pca, group = data_group$cell.line, ind.names = FALSE,
          legend = TRUE, 
          title = 'PCA')
ggsave(paste0(dir_DE,"PCA.pdf"))

plsda <- mixOmics::plsda(X,Y, ncomp = 2)
vip <- mixOmics::vip(plsda)

# 交叉验证评估（稳定性能指标）
set.seed(123)  
perf.plsda <- mixOmics::perf(plsda, validation = "Mfold",folds = 7,# n<3,时，k=2n
                             nrepeat = 50, progressBar = FALSE) # TRUE显示进度
plot(perf.plsda, sd = TRUE, legend.position = 'horizontal')

final.plsda <- plsda(X,Y, ncomp = 2)

plotIndiv(final.plsda, ind.names = FALSE, legend=TRUE,
          comp=c(1,2), ellipse = TRUE, 
          title = 'PLS-DA',
          X.label = 'comp 1', Y.label = 'comp 2')

### 4.2.3 method 2 ropls -------------------------------------------------------
library(ropls)
# PLS-DA 分析
set.seed(123)
plsda_model <- opls(X, Y, predI = 2, crossvalI = 5, permI = 200, scaleC = "pareto")
# scaleC = c("pareto","center","standard","none")
# PLS-DA
# 5 samples x 1120 variables and 1 response
# standard scaling of predictors and response(s)
# R2X(cum) R2Y(cum) Q2(cum)   RMSEE pre ort  pR2Y   pQ2
# Total    0.645        1   0.966 0.00685   2   0 0.355 0.135

# 提取 VIP 值
vip_scores <- plsda_model@vipVn
table(vip_scores > 1)
summary(vip_scores)

## 4.3 Res output ---------------------------------------------------------------
dir_DE <- "./03_Result/2.DE/combined/MOLM13/Low_vs_Con/"
result_merge <- read.csv("./03_Result/2.DE/combined/OCI_M2/High_vs_Con/DE_results.csv",row.names = 1)
result_merge$VIP <- vip_scores 

# 将VIP<1 的change列改为 stable
table(result_merge$change)
result_merge <- result_merge %>%
  mutate(change = ifelse(VIP < 1, "stable", change))
table(result_merge$change)
write.csv(result_merge, file = paste0(dir_DE,"DE_results.csv"))
save(plsda_model, file = paste0(dir_DE,"plsda_model.rds"))

## 4.4 Volc Plot ---------------------------------------------------------------
# data input 
result_merge <- read.csv("./03_Result/2.DE/combined/MOLM13/Low_vs_Con/DE_results.csv",row.names = 1)

group_1 <- "Low"        # treatment
group_2 <- "Con"        # control
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
  geom_point(aes(color = change, size = VIP),  
             alpha = 0.6, shape = 16) +
  ylab("-log10(Pvalue)")+
  # 按因子顺序指定颜色
  scale_color_manual(
    name = "Change",                              # 图例标题
    values = c(
      "up" = "#B30000",      # 红
      "stable" = "grey",     # 灰
      "down" = "#003366"),    # 蓝,
    labels = c(                                  # 显示上下调的个数
      "up" = paste0("Up:", sum(result_merge$change == "up")),
      "stable" = paste0("Stable:", sum(result_merge$change == "stable")),
      "down" = paste0("Down:", sum(result_merge$change == "down"))))+
  scale_size_continuous(name = "VIP",        # 为 size 变量添加图例标题
                        range = c(1, 6))+    # 设置点的大小范围，可以根据需要调整
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), lty = 4, 
             col = "black", lwd = 0.8,alpha=0.4) +
  geom_hline(yintercept = -log10(pvalue_threshold), lty = 4, 
             col = "black", lwd = 0.8,alpha=0.4) +
  labs(title = paste0(group_1,"-",group_2)) +
  xlim(-x_limit, x_limit)+
  theme_bw() +
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        aspect.ratio = 1.2,                      # 设置纵横比，调整为更高
        axis.title = element_text(size = 12, face = "plain"),  # 轴标题加粗，字号14
        # 图例标题 
        legend.title = element_text(
          family = "serif",        # 字体类型（如serif衬线体）
          size = 12,               # 字体大小
          face = "bold",           # 加粗
          color = "#333333",       # 深灰色
          margin = margin(b = 5)   # 标题与标签的间距
        ),
        
        # 图例标签 
        legend.text = element_text(
          family = "sans",         # 无衬线字体（如Arial）
          size = 11,               
          color = "#666666",       # 中灰色
          margin = margin(r = 10)  # 标签右侧间距
        ),
        
        # 整体布局 
        legend.position = "right",
        legend.spacing.y = unit(0.3, "cm"),  # 图例项垂直间距
        legend.key.size = unit(0.8, "cm"),   # 图例符号大小
        legend.background = element_rect(fill = "white", color = NA))  # 背景优化  
print(p)
ggsave(filename = paste0(dir_DE,"/volc.pdf"),
       plot = p, device = cairo_pdf, 
       width = 6, height = 5)

# 5. KEGG GO -----------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)  # 人类数据库，其他物种可更换
library(KEGGREST)
library(pathview)
library(MetaboAnalystR)

# MetaboAnalyst操作建议在在线网站上操作，本地也需要连接服务器，且不稳定
# 本地通路富集用Clusterprofiler包，需要提供背景代谢物集

# load data
DEM_result <- read.csv("./03_Result/2.DE/combined/OCI_M2/High_vs_Con/DE_results.csv",row.names = 1)
BGM <- DEM_result[,c("Name","change")]
DEM <- DEM_result[DEM_result$change != "stable",c("Name","change")]
write.xlsx(DEM, file = "./03_Result/2.DE/combined/OCI_M2/High_vs_Con/DEM_name.xlsx")
write.xlsx(BGM,file = "./03_Result/2.DE/combined/BGM_name.xlsx")

# MetaboAnalyst 的ID转换 or 在线网站 https://www.metaboanalyst.ca  
# mSet<-InitDataObjects("conc", "msetora", FALSE)
# cmpd.vec<-DEM_result$Name
# mSet<-Setup.MapData(mSet, cmpd.vec)
# mSet<-CrossReferencing(mSetObj = mSet, q.type = "name", hmdb = T, pubchem = T, kegg = T)
 mSet<-CreateMappingResultTable(mSet)
# mSet<-SetMetabolomeFilter(mSet, F)

## BgM list ---- 
BGM_map <- mSet$dataSet$map.table
BGM <- as.data.frame(BGM)

## DEM list ----
DEM_map <- read.csv("./03_Result/2.DE/combined/OCI_M2/High_vs_Con/name_map.csv")
pathway <- pathway[,c("PathwayName","Metabolite")]


## Enrichment ----
# 可以使用在线网站 https://www.metaboanalyst.ca 
# 富集分析
kegg_res <- clusterProfiler::enricher(gene = DEM_map$KEGG,
                                      TERM2GENE = pathway,
                                      minGSSize = 1,
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)
# 结果导出
write.csv(as.data.frame(x@result) %>% select(-1,-2),
          file = paste0(uploadfile2,"/MP_KEGG_enrichment_result.csv"))
