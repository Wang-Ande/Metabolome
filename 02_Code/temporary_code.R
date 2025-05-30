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
  
