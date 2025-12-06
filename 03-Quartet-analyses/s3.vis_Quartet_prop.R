#!/usr/bin/env Rscript

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)


setwd("/Users/xiaojun/Nutstore Files/ILS-Sulfitobacter/intermediates/20250415/Quartet_analysis")

######### Vis for HR quartet ######
# 函数：处理2-2分布的数据
process_two_clade_2_2 <- function(data) {
  # 初始化结果数据框
  result <- data.frame(
    QUARTET_TYPE = character(),
    against_HR = numeric(),
    support_HR = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 遍历每一行数据
  for (i in 1:nrow(data)) {
    # 获取当前行的数据
    row <- data[i, ]
    clades <- strsplit(row$QUARTET_TYPE, "-")[[1]]
    
    # 计算每个topology的类型
    ab_cd_type <- ifelse(
      (clades[1] != clades[2] && clades[3] != clades[4]),
      "support",
      ifelse(
        (clades[1] == clades[2] && clades[3] == clades[4] && clades[1] != clades[3]),
        "against",
        "support"
      )
    )
    
    ac_bd_type <- ifelse(
      (clades[1] != clades[3] && clades[2] != clades[4]),
      "support",
      ifelse(
        (clades[1] == clades[3] && clades[2] == clades[4] && clades[1] != clades[2]),
        "against",
        "support"
      )
    )
    
    ad_bc_type <- ifelse(
      (clades[1] != clades[4] && clades[2] != clades[3]),
      "support",
      ifelse(
        (clades[1] == clades[4] && clades[2] == clades[3] && clades[1] != clades[2]),
        "against",
        "support"
      )
    )
    
    # 计算against和support的总数
    against_HR <- sum(
      ifelse(ab_cd_type == "against", row$COUNT.ab.cd., 0),
      ifelse(ac_bd_type == "against", row$COUNT.ac.bd., 0),
      ifelse(ad_bc_type == "against", row$COUNT.ad.bc., 0)
    )
    
    support_HR <- row$total - against_HR
    
    # 添加到结果数据框
    result <- rbind(result, data.frame(
      QUARTET_TYPE = row$QUARTET_TYPE,
      against_HR = against_HR,
      support_HR = support_HR,
      stringsAsFactors = FALSE
    ))
  }
  
  return(result)
}

# 函数：处理2-1-1分布的数据
process_three_clade <- function(data) {
  # 初始化结果数据框
  result <- data.frame(
    QUARTET_TYPE = character(),
    against_HR = numeric(),
    support_HR = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 遍历每一行数据
  for (i in 1:nrow(data)) {
    # 获取当前行的数据
    row <- data[i, ]
    clades <- strsplit(row$QUARTET_TYPE, "-")[[1]]
    
    # 计算每个topology的类型
    ab_cd_type <- ifelse(
      (clades[1] == clades[2]) || (clades[3] == clades[4]),
      "against",
      "support"
    )
    
    ac_bd_type <- ifelse(
      (clades[1] == clades[3]) || (clades[2] == clades[4]),
      "against",
      "support"
    )
    
    ad_bc_type <- ifelse(
      (clades[1] == clades[4]) || (clades[2] == clades[3]),
      "against",
      "support"
    )
    
    # 计算against和support的总数
    against_HR <- sum(
      ifelse(ab_cd_type == "against", row$COUNT.ab.cd., 0),
      ifelse(ac_bd_type == "against", row$COUNT.ac.bd., 0),
      ifelse(ad_bc_type == "against", row$COUNT.ad.bc., 0)
    )
    
    support_HR <- row$total - against_HR
    
    # 添加到结果数据框
    result <- rbind(result, data.frame(
      QUARTET_TYPE = row$QUARTET_TYPE,
      against_HR = against_HR,
      support_HR = support_HR,
      stringsAsFactors = FALSE
    ))
  }
  
  return(result)
}

# 读取并处理数据
data_2_2 <- read.table("two_clade_2_2_results.tsv", header=TRUE, sep="\t") %>%
  mutate(
    total = COUNT.ab.cd. + COUNT.ac.bd. + COUNT.ad.bc.,
    # 确保QUARTET_TYPE是字符类型
    QUARTET_TYPE = as.character(QUARTET_TYPE)
  ) %>%
  process_two_clade_2_2() %>%
  mutate(
    # 先计算总数
    total_HR = against_HR + support_HR,
    # 然后计算比例
    against_HR_prop = against_HR / total_HR,
    support_HR_prop = support_HR / total_HR
  ) %>%
  select(QUARTET_TYPE, against_HR = against_HR_prop, support_HR = support_HR_prop)

data_3 <- read.table("three_clade_results.tsv", header=TRUE, sep="\t") %>%
  mutate(
    total = COUNT.ab.cd. + COUNT.ac.bd. + COUNT.ad.bc.,
    # 确保QUARTET_TYPE是字符类型
    QUARTET_TYPE = as.character(QUARTET_TYPE)
  ) %>%
  process_three_clade() %>%
  mutate(
    # 先计算总数
    total_HR = against_HR + support_HR,
    # 然后计算比例
    against_HR_prop = against_HR / total_HR,
    support_HR_prop = support_HR / total_HR
  ) %>%
  select(QUARTET_TYPE, against_HR = against_HR_prop, support_HR = support_HR_prop)

# 打印数据以检查
print("Data 2-2:")
print(head(data_2_2))
print("Data 3:")
print(head(data_3))

# 合并数据
combined_data <- bind_rows(
  data_2_2 %>% mutate(type = "2-2"),
  data_3 %>% mutate(type = "2-1-1")
)

# 转换为长格式
data_long <- combined_data %>%
  pivot_longer(
    cols = c(against_HR, support_HR),
    names_to = "topology",
    values_to = "proportion"
  ) %>%
  mutate(
    topology = factor(
      topology,
      levels = c("against_HR", "support_HR"),
      labels = c("against HR", "support HR")
    )
  )

# 打印转换后的数据
print("Combined data:")
print(head(data_long))

# 定义颜色
colors <- c(
  "against HR" = "gray",
  "support HR" = "#fb8702"
)

# 绘制图形
p <- ggplot(data_long, aes(x = QUARTET_TYPE, y = proportion, fill = topology)) +
  geom_bar(stat = "identity", position = "stack", size = 0.2) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(labels = scales::percent) +  # 添加百分比标签
  labs(
    y = "Proportion",
    fill = "Topology"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# 保存图形
pdf(file="quartet_prop4HR.pdf", height=8, width=12)
p
dev.off()

### check quartets showing higher HR supporting
data_wide <- data_long %>%  
  pivot_wider(  
    names_from = "topology",  
    values_from = "proportion"  
  )  

colnames(data_wide)[3:4] <- c("against","support")

median(data_wide$support)

# Filter rows where any of the percentages exceed 70%  
filtered_data <- data_wide %>%  
  filter( support > 0.1 )  


#######Vis for Phylogenetic Incongruence quartet #####
# 读取数据
data <- read.table("four_clade_results.tsv", header=TRUE, sep="\t")

# 计算总数和百分比
data <- data %>%
  mutate(
    total = COUNT.ab.cd. + COUNT.ac.bd. + COUNT.ad.bc.,
    ab_cd_pct = COUNT.ab.cd. / total * 100,
    ac_bd_pct = COUNT.ac.bd. / total * 100,
    ad_bc_pct = COUNT.ad.bc. / total * 100
  )

# 转换数据为长格式
data_long <- data %>%
  select(QUARTET_TYPE, ab_cd_pct, ac_bd_pct, ad_bc_pct) %>%
  pivot_longer(
    cols = c(ab_cd_pct, ac_bd_pct, ad_bc_pct),
    names_to = "topology",
    values_to = "percentage"
  ) %>%
  mutate(
    topology = factor(
      topology,
      levels = c("ab_cd_pct", "ac_bd_pct", "ad_bc_pct"),
      labels = c("(A,B),(C,D)", "(A,C),(B,D)", "(A,D),(B,C)")
    )
  )

# 定义颜色
colors <- c(
  "(A,B),(C,D)" = "#8971d0",
  "(A,C),(B,D)" = "#ffc93c",
  "(A,D),(B,C)" = "#17b978"
)

# 绘制图形
p <- ggplot(data_long, aes(x = QUARTET_TYPE, y = percentage, fill = topology)) +
  geom_bar(stat = "identity", position = "stack", size = 0.2) +
  scale_fill_manual(values = colors) +
  labs(
    y = "Proportion (%)",
    fill = "Topology"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# 保存图形
pdf(file="quartet_prop4PI.pdf", height=8, width=12)
p
dev.off()

### check quartets showing extreme proportions
data_wide <- data_long %>%  
  pivot_wider(  
    names_from = "topology",  
    values_from = "percentage"  
  )  

colnames(data_wide)[2:4] <- c("ab_cd_pct","ac_bd_pct","ad_bc_pct")
# Filter rows where any of the percentages exceed 70%  
filtered_data <- data_wide %>%  
  filter(  
    ab_cd_pct > 70 |   
      ac_bd_pct > 70 |   
      ad_bc_pct > 70  
  )  
