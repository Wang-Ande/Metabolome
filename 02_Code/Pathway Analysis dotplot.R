#' 此脚本仅用于定制 MetaboAnalyst 网站中 Pathway Analysis Module 的结果可视化
#' 为了避免错误，初学者尽量不要更改下面的变量名

# packages ----
library(ggplot2)

# pathway res input ----
pathway_res <- read.csv("your result")
colnames(pathway_res)[1] <- "Pathway"
pathway_res$logP <- -log10(pathway_res$Raw.p)

# plot ----
p1 <- ggplot(pathway_res, aes(x = Impact, y = logP)) +
  geom_point(aes(size = Hits, fill = logP),
             color = "black", shape = 21, stroke = 0.7) +  # 黑色边框，shape支持fill 
  scale_fill_gradientn(colours = c("yellow", "orange", "red"),
                       limits = c(0, 2),  # 控制颜色渐变范围
                       oob = scales::squish) +
  scale_size(range = c(3, 8)) +
  labs(x = "Pathway Impact",
       y = "-log10(P)",
       color = "-log10(P)",
       size = "Hits",
       title = "Metabolic Pathway Analysis") +
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        legend.background = element_rect(fill = "white", color = NA),  # 背景优化
        legend.position = "none")
print(p1)

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
  scale_fill_gradientn(colours = c("yellow", "orange", "red"),
                       limits = c(0, 2),       # 控制颜色渐变范围(自己调试)
                       oob = scales::squish) + # 超过渐变范围的截断为最大值
  scale_size(range = c(3, 8)) +                # 气泡的大小范围
  labs(x = "Pathway Impact", y = "-log10(P)",
       color = "-log10(P)", size = "Hits",
       title = "Metabolic Pathway Analysis") +
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5),  # 标题居中
        plot.margin = margin(t = 2, b = 2, r = 5, l = 2, unit = "mm"),  # 上下左右边距
        legend.background = element_rect(fill = "white", color = NA),  # 背景优化
        legend.position = "none") # 不显示图例，这个图不显示图例整体效果更好看
print(p2)
# output ----
ggsave(
  filename = paste0("Pathway_analysis.pdf"),    
  plot = p2,            # 要保存的图形对象
  device = cairo_pdf,
  scale = 1,             # 缩放比例（相对于默认尺寸）
  width = 6.14,          # 图像宽度（单位：英寸或厘米，取决于 units）
  height = 6.31,         # 图像高度
  units = "in")          # 尺寸单位c("in", "cm", "mm", "px")
