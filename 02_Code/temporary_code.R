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


pathway_res <- read.csv("./03_Result/3.Enrichment/OCI_M2/High_vs_Con/pathway_results.csv")
colnames(pathway_res)[1] <- "Pathway"
pathway_res$logP <- -log10(pathway_res$Raw.p)

# 绘制气泡图
p1 <- ggplot(pathway_res, aes(x = Impact, y = logP)) +
  geom_point(aes(size = Hits, fill = logP),
             color = "black", shape = 21, stroke = 0.7) +  # 黑色边框，shape支持fill 
  scale_fill_gradientn(
    colours = c("yellow", "orange", "red"),
    limits = c(0, 2),  # 控制颜色渐变范围
    oob = scales::squish) +
  scale_size(range = c(3, 8)) +
  labs(
    x = "Pathway Impact",
    y = "-log10(P)",
    color = "-log10(P)",
    size = "Hits",
    title = "Metabolic Pathway Analysis"
  ) +
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        legend.background = element_rect(fill = "white", color = NA),  # 背景优化
        legend.position = "none")
print(p1)
ggsave(
  filename = paste0(dir_pca,"Pca_phase_Track.pdf"),    
  plot = p1,            # 要保存的图形对象
  device = cairo_pdf,
  scale = 1,           # 缩放比例（相对于默认尺寸）
  width = 7.91,          # 图像宽度（单位：英寸或厘米，取决于 units）
  height = 6.58,         # 图像高度
  units = "in")          # 尺寸单位c("in", "cm", "mm", "px")

# 标注显著通路
library(ggrepel)

# 设置显著通路的标签（如前5个）
top_pathways <- head(pathway_res, 6)

p2 <- ggplot(pathway_res, aes(x = Impact, y = logP)) +
  geom_point(aes(size = Hits, fill = logP),
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
  scale_size(range = c(3, 8)) +
  labs(
    x = "Pathway Impact",
    y = "-log10(P)",
    color = "-log10(P)",
    size = "Hits",
    title = "Metabolic Pathway Analysis"
  ) +
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        plot.margin = margin(t = 2, b = 2, r = 5, l = 2, unit = "mm"),  # 上下左右边距
        legend.background = element_rect(fill = "white", color = NA),  # 背景优化
        legend.position = "none")
print(p2)
ggsave(
  filename = paste0("03_Result/3.Enrichment/OCI_M2/High_vs_Con/Pathway_analysis.pdf"),    
  plot = p2,            # 要保存的图形对象
  device = cairo_pdf,
  scale = 1,           # 缩放比例（相对于默认尺寸）
  width = 6.14,          # 图像宽度（单位：英寸或厘米，取决于 units）
  height = 6.31,         # 图像高度
  units = "in")          # 尺寸单位c("in", "cm", "mm", "px")
