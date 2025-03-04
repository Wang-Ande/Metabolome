run_DE <- function(data, data_group, data_anno=NULL, group_1, group_2, log2,
                   logfc_threshold, pvalue_threshold, qvalue_threshold = NULL,
                   test_method, paired ,dir = getwd()) {
  #  基础检查 ----
  library(fdrtool)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(mixOmics)
  if (!all(c("id", "group") %in% colnames(data_group)))
    stop("data_group must contain 'id' and 'group' columns")
  if (!all(group_1 %in% data_group$group) | !all(group_2 %in% data_group$group))
    stop("Specified groups not found in data_group")
  
  #  创建输出目录 ----
  output_dir <- file.path(dir, paste0(group_1, "_vs_", group_2))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  #  提取目标样本 ----
  samples <- data_group$id[data_group$group %in% c(group_1, group_2)]
  target_data <- data[, colnames(data) %in% samples, drop=FALSE]
  
  # log2
  if(log2){
  target_data <- log2(target_data)
  }
  
  #  分组样本获取 ----
  group_samples <- list(
    group1 = data_group$id[data_group$group == group_1],
    group2 = data_group$id[data_group$group == group_2]
  )
  
  #  检查样本存在性 ----
  group_samples <- lapply(group_samples, function(x) x[x %in% colnames(target_data)])
  if (any(lengths(group_samples) < 2)) 
    stop("At least 2 samples required in each group for t-test")
  
  #  准备结果数据框 ----
  result_df <- data.frame(
    meanA = rowMeans(target_data[, group_samples$group1], na.rm = TRUE),
    meanB = rowMeans(target_data[, group_samples$group2], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  #  computer fold change ----
  if(log2){
  result_df$logFC <- result_df$meanA - result_df$meanB
  }
  else{
    result_df$FC <- result_df$meanA / result_df$meanB
    result_df$logFC <- log2(result_df$FC)
      }
    
  # 处理特殊值
  #result_df$log2foldchange[is.na(result_df$log2foldchange)] <- 0
  #result_df$log2foldchange[is.infinite(result_df$log2foldchange)] <- 0
  
  #  选择检验方法 ----
  safe_stat_test <- function(x) {
    g1 <- x[group_samples$group1]
    g2 <- x[group_samples$group2]
    if (length(g1) < 2 | length(g2) < 2) return(NA)  # 提前检查样本数
    
    if (test_method == "t-test") {
      # t检验 （要求符合正态分布）
      tryCatch(t.test(g1, g2, paired = paired)$p.value, error = function(e) NA)
    } else if (test_method == "wilcoxon") {
      # Wilcoxon秩和检验 (小样本且没有ties时，选择exact = TRUE)
      tryCatch(wilcox.test(g1, g2, paired = paired, exact = TRUE)$p.value, error = function(e) NA)
    } else {
      stop("Invalid test method. Choose either 't-test' or 'wilcoxon'.")
    }
  }
  
  result_df$pvalue <- apply(target_data, 1, safe_stat_test)
  
  #  多重检验校正 ----
  if (requireNamespace("fdrtool", quietly = TRUE)) {
    qobj <- fdrtool::fdrtool(result_df$pvalue, statistic = "pvalue", plot = FALSE)
    result_df$qvalue <- qobj$qval
  } else {
    result_df$qvalue <- p.adjust(result_df$pvalue, method = "BH")
  }
  
  
  #  加change列 ---- 
  # 标记上下调基因，可根据需求设定阈值
  logFC = logfc_threshold
  P.Value = pvalue_threshold
  k1 <- (result_df$pvalue < P.Value) & (result_df$logFC < -logFC)
  k2 <- (result_df$pvalue < P.Value) & (result_df$logFC > logFC)
  result_df <- mutate(result_df, 
                           change = ifelse(k1, "down", 
                                           ifelse(k2, "up", "stable")))
  table(result_df$change)
  
  #  合并定量信息 ----
  # Step 1: 检查行名是否相同
  if (!all(rownames(target_data) %in% rownames(result_df))) {
    warning("Some row names in 'target_data' are not found in 'result_df'. Only matching rows will be used.")
  }
  
  # Step 2: 筛选 target_data 中存在于 result_df 的行
  target_data <- 2^target_data
  target_data_filtered <- target_data[rownames(target_data) %in% rownames(result_df), , drop = FALSE]
  
  # Step 3: 按 result_df 的行顺序合并
  result_df <- cbind(target_data_filtered[match(rownames(result_df), rownames(target_data_filtered)), ], result_df)
  
  #  合并注释信息 ----
  # Step 1: 检查行名是否相同
  if (!all(rownames(data_anno) %in% rownames(result_df))) {
    warning("Some row names in 'data_anno' are not found in 'result_df'. Only matching rows will be used.")
  }
  
  # Step 2: 筛选 data_anno 中存在于 result_df 的行
  data_anno_filtered <- data_anno[rownames(data_anno) %in% rownames(result_df), , drop = FALSE]
  
  # Step 3: 按 result_df 的行顺序合并
  result_df <- cbind(data_anno_filtered[match(rownames(result_df), rownames(data_anno_filtered)), ], result_df)
  
  #  Res Output ----
  write.csv(result_df, file = file.path(output_dir, "DE_results.csv"))
  
  
  return(result_df)
}
