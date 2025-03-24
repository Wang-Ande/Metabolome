# 加载数据
pos_data <- read.xlsx("./01_Data/01_MetQuant/meta_intensity_class_pos.xlsx")
neg_data <- read.xlsx("./01_Data/01_MetQuant/meta_intensity_class_neg.xlsx")
colnames(pos_data) <- gsub("w","W",colnames(pos_data))
colnames(neg_data) <- gsub("w","W",colnames(neg_data))
rownames(pos_data) <- pos_data$Compound_ID
rownames(neg_data) <- neg_data$Compound_ID

# 添加模式标识
pos_data$Ion_Mode <- "POS"
neg_data$Ion_Mode <- "NEG"

# view intersect
common_meta <- intersect(pos_data,neg_data)

# check colnames
colnames(pos_data) == colnames(neg_data)

# merge data
combined_data <- rbind(pos_data, neg_data)

# 去除重复代谢物
combined_data <- combined_data[!duplicated(combined_data[, c("mz", "RT")]), ]

# 导出合并后的数据
write.csv(combined_data, file = "./01_Data/01_MetQuant/meta_intensity_combined.csv")


