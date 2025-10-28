# temporary code

## 4.2 PLS-DA ----
if(T){
  # 确保数据矩阵没有零方差特征（重要！否则会导致plsda报错）
  samples <- targeted_group$id[targeted_group$group %in% c(group_1, group_2)]
  target_data <- data_input[, colnames(data_input) %in% samples, drop=FALSE]
  target_data <- log2(target_data)
  
  row_vars <- apply(target_data, 1, var, na.rm = TRUE)
  if(any(row_vars == 0)){
    zero_var_features <- sum(row_vars == 0)
    warning(paste("Removing", zero_var_features, "features with zero variance"))
    target_data <- target_data[row_vars > 0, ]
  }
  
  # 模型评估和验证
  # 模型构建与交叉验证
  plsda_model <- plsda(
    X = t(target_data), 
    Y = factor(targeted_group$group[targeted_group$id %in% samples]), 
    ncomp = 2,                    
    scale = TRUE,                 # 强制标准化数据
    near.zero.var = TRUE          # 自动处理接近零方差特征
  )
  
  # VIP
  vip_scores <- vip(plsda_model)
  
  # 确保行名匹配（重要！）
  if(!all(rownames(target_data) %in% rownames(vip_scores))){
    stop("Feature names mismatch between target_data and VIP scores")
  }
  
  # 合并VIP值（默认使用第一主成分）
  result_merge$VIP <- vip_scores[rownames(result_merge), 1] 
  write.csv(result_merge,file = paste0(dir_DE, "DE_results.csv"))
  
  # PLS-DA Score Plot 
  plsda_scores <- as.data.frame(plsda_model$variates$X)
  plsda_scores$Group <- factor(targeted_group$group[targeted_group$id %in% samples])
  #plsda_scores$Sample <- data_group$id[data_group$id %in% samples]
  
  # 提取主成分解释方差比例（新增）
  explained_var <- plsda_model$prop_expl_var$X * 100
  
  p_plsda <- ggplot(plsda_scores, aes(x = comp1, y = comp2, color = Group)) +
    geom_point(size = 4, alpha = 0.8, shape = 17) +  # 改为三角形符号
    stat_ellipse(level = 0.95, linewidth = 0.8) +    # 加粗椭圆线
    geom_text_repel(aes(label = targeted_group$id[targeted_group$id %in% samples]), 
                    size = 3, max.overlaps = Inf) +  # 使用 geom_text_repel 避免标签重叠
    labs(
      x = paste0("Component 1 (", round(explained_var[1],1), "%)"),
      y = paste0("Component 2 (", round(explained_var[2],1), "%)"),
      title = paste("PLS-DA Score Plot:", group_1, "vs", group_2)
    ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +  # 自定义颜色
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      aspect.ratio = 1
    )
  
  ggsave(filename = file.path(dir_DE,"PLSDA_score_plot.pdf"), 
         plot = p_plsda, device = "pdf", width = 6, height = 5, dpi = 300)
}
  # PLS-DA Loading Plot 
  
# Pathway analysis ----
library(ggplot2)
# example
df <- data.frame(
  pathway = c("Glycolysis", "TCA cycle", "Glutathione metabolism"),
  Impact = c(0.56, 0.82, 0.13),
  pvalue = c(0.001, 0.045, 0.2),
  Hits = c(6, 3, 1) 
)
df$logP <- -log10(data$pvalue)


pathway_res <- read.csv("./03_Result/4.Pathway analysis/OCI_AML2/High_vs_Con/OCI_Pathway_analysis_up_DE_KEGG/pathway_results.csv")
colnames(pathway_res)[1] <- "Pathway"
pathway_res$logP <- -log10(pathway_res$Raw.p)

# 绘制气泡图
p1 <- ggplot(pathway_res, aes(x = Impact, y = logP)) +
  geom_point(aes(size = Impact, fill = logP),
             color = "black", shape = 21, stroke = 0.7) +  # 黑色边框，shape支持fill 
  scale_fill_gradientn(
    colours = c("yellow", "orange", "red"),
    limits = c(0, 2),  # 控制颜色渐变范围
    oob = scales::squish) +
  scale_size(range = c(4, 8)) +
  labs(
    x = "Pathway Impact",
    y = "-log10(P)",
    color = "-log10(P)",
    size = "Impact",
    title = "Metabolic Pathway Analysis"
  ) +
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        legend.background = element_rect(fill = "white", color = NA),  # 背景优化
        legend.position = "none")
print(p1)

# 标注显著通路
library(ggrepel)

# 设置显著通路的标签（如前5个）
top_pathways <- head(pathway_res, 7)

p2 <- ggplot(pathway_res, aes(x = Impact, y = logP)) +
  geom_point(aes(size = Impact, fill = logP),
             color = "black", shape = 21, stroke = 0.7) +  # 黑色边框，shape支持fill 
  geom_text_repel(data = top_pathways, aes(label = Pathway), 
                  box.padding = 0.4,           # 标签与点的间距
                  force = 0.7,                 # 避让算法的力度（值越大，标签越分散）
                  max.overlaps = 20,           # 允许的最大重叠次数
                  max.time = 1,                # 计算避让的最大时间（秒）
                  max.iter = 1e4,              # 迭代次数上限
                  size = 4.5) +
  scale_fill_gradientn(
    colours = c("yellow", "orange", "red"),
    limits = c(0, 2),  # 控制颜色渐变范围
    oob = scales::squish) +
  scale_size(range = c(3, 12)) +
  labs(
    x = "Pathway Impact",
    y = "-log10(P)",
    color = "-log10(P)",
    size = "Impact",
    title = NULL  # "Metabolic Pathway Analysis"
  ) +
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        plot.margin = margin(t = 2, b = 2, r = 5, l = 2, unit = "mm"),  # 上下左右边距
        legend.background = element_rect(fill = "white", color = NA),  # 背景优化
        legend.position = "none",
        panel.grid.minor = element_blank(),     # 去除次级网格线
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),  # 主网格线仅保留在主刻度
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
        )
print(p2)

ggsave(
  filename = paste0("03_Result/4.Pathway analysis/OCI_AML2/High_vs_Con/Pathway_analysis.pdf"),    
  plot = p2,            # 要保存的图形对象
  device = cairo_pdf,
  dpi = 300,
  scale = 1,           # 缩放比例（相对于默认尺寸）
  width = 6,          # 图像宽度（单位：英寸或厘米，取决于 units）
  height = 6,         # 图像高度
  units = "in")          # 尺寸单位c("in", "cm", "mm", "px")

# DE ----
# 注意，data和data_anno的行名应一致
# 根据分组选择要进行差异分析的组别
library(openxlsx)
library(readr)
## 1.1 Group input ----
# 导入分组信息
data_group <- read.xlsx("./01_Data/01.MetQuant/sam_infor_combined.xlsx")
data_group <- as.data.frame(data_group)
# 去除
data_group <- data_group[-grep(""),]


## 1.2 meta matrix input ----
data_input <- read.csv("./01_Data/01.MetQuant/meta_intensity_combined.csv",row.names = 1)
data_input <- as.data.frame(data_input)
colnames(data_input) <- gsub("neg_","", colnames(data_input))   # 去除样本id多余信息
colnames(data_input) <- gsub("cas_","", colnames(data_input))
data_anno <- read.xlsx("./01_Data/01.MetQuant/meta_anno_combined.xlsx",rowNames = TRUE)

source("./02_Code/run_DE.R")
table(data_group$group)
targeted_group <- data_group[grep("6W|WT",data_group$id),]
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
                       dir = "03_result/2.DE/combined/All/") # 每次需要更改
# 统计上下调Meta个数
table(result_merge$change)
# 导出差异代谢物列表
DE_Metabolite <- read.csv('./03_Result/2.DE/combined/All/High_vs_Con/DE_results.csv')
DE_Metabolite <- DE_Metabolite[,c("Name","logFC","pvalue","qvalue","change")]
write.xlsx(DE_Metabolite, file = "./03_Result/2.DE/combined/All/High_vs_Con/DE_Metabolite_Names.xlsx")

# 注释火山图 ----
library(openxlsx)
library(ggplot2)
library(ggrepel)  # 避免标签重叠

DE_res <- read.xlsx("./03_Result/2.DE/combined/All/High_vs_Con/DE_Metabolite_Names.xlsx")
DE_res$Sig <- factor(DE_res$Sig, levels = c("up", "stable"))
max_abs_logfc <- max(abs(DE_res$logFC), na.rm = TRUE)
x_limit_right <- max_abs_logfc * 1.05
# 如果要求对称则不需要设置左界限，直接设置为 -x_limit_right 即可
x_limit_left <- min(DE_res$logFC) * 1.2

# 输出目录
output_dir <- "./03_Result/2.DE/combined/All/High_vs_Con/"

# 标记分组
group_1 <- "VR"    
group_2 <- "WT"        

# 标记阈值
logfc_threshold <- 0
pvalue_threshold <- 0.05


# 提取显著差异基因在火山图上标记
# 比如 top 10 up 
top_up <- DE_res[DE_res$change == "up", ]
# 复杂名称替换成简单名称
top_up[grep("1-Palmitoyl-Sn-Glycero-3-Phosphocholine",top_up$Name),1] <- "Lyso-PC(16:0)"
top_up[grep("5-Hydroxytryptophan",top_up$Name),1] <- "5-HTP"
not_show <- c("N~1~-(2,4-dichlorophenyl)-N~2~-(2,2-dimethoxyethyl)ethanediamide",
              "(3S,9aS)-3-benzyl-octahydro-1H-pyrido[1,2-a]pyrazin-1-one",
              "1-[5-(2-phenyleth-1-ynyl)-2-thienyl]ethan-1-one oxime",
              "(2E,4E)-N-[2-(4-hydroxyphenyl)ethyl]dodeca-2,4-dienamide",
              "Cyclohexylsulfamate",
              "Boc-beta-cyano-L-alanine",
              "LPS 19:1")
top_up <- top_up[!top_up$Name%in%not_show,]
top_up <- top_up[order(top_up$logFC,decreasing = TRUE), ][1:10, ]
# top_down <- DE_res[DE_res$Sig == "down", ][order(DE_res$P.Value), ][1:10, ]
label_data <- top_up

# 控制图例顺序
DE_res$change <- factor(DE_res$change, 
                        levels = c("up", "down", "stable"))

## plot ----
p1 <- ggplot(data = DE_res, 
             aes(x = logFC, 
                 y = -log10(pvalue))) +
  geom_point(alpha = 0.5, size = 1.5, 
             aes(color = change)) +
  ylab("-log10(P.value)")+
  scale_color_manual(
    name = "Change",
    values = c("up" = "#B30000", "stable" = "grey", "down" = "#003366"),
    labels = c("up" = paste0("Up ：", sum(DE_res$change == "up")),
               "stable" = paste0("Stable ：", sum(DE_res$change == "stable")),
               "down" = paste0("Down ：", sum(DE_res$change == "down"))))+
  #geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), lty = 4, 
             #col = "black", lwd = 0.8, alpha = 0.4) +
  geom_hline(yintercept = -log10(pvalue_threshold), lty = 4, 
             col = "black", lwd = 0.8, alpha = 0.4) +
  geom_text_repel(data = label_data,   # 标签
                  aes(label = Name), 
                  size = 3.5,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  arrow = arrow(length = unit(0.008,"npc")), # 箭头指向
                  segment.size = 0.5, min.segment.length = 0.5,            # 确保箭头显示
                  force =15, # 标签排斥力
                  max.overlaps = 20) +
  labs(title = paste0(group_2,"-",group_1)) +
  # xlim(x_limit_left, x_limit_right)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), # x轴网格线
        panel.grid.minor = element_blank(), # y轴网格线
        aspect.ratio = 1.2)
print(p1)
ggsave(filename = paste0(output_dir,"anno_volc.pdf"),
       plot = p1, device = "pdf", 
       width = 6.5, height = 5)
