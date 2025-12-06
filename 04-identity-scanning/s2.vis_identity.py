import pandas as pd  
import matplotlib.pyplot as plt  
import seaborn as sns  

# 读取 Identity 结果  
identity_file = "pop_identity_results_step.csv"  
data = pd.read_csv(identity_file)  

# 打印数据以确认格式  
print(data.head())  

# 准备数据：将组合标准化为比较标签  
data['Comparison'] = data.apply(lambda row: f"{min(row['Population A'], row['Population B'])}_vs_{max(row['Population A'], row['Population B'])}", axis=1)  

# 设置 Seaborn 样式  
sns.set(style="whitegrid", font_scale=1.2)  

# 创建 FacetGrid，以每个种群组合作为列绘制  
g = sns.FacetGrid(data, col="Comparison", col_wrap=3, height=4, sharey=False)  

# 自定义绘图函数  
def draw_plot(window_start, mean_identity, **kwargs):  
    ax = plt.gca()  # 获取当前的轴  
    sns.lineplot(x=window_start, y=mean_identity, linewidth=1.5, ax=ax, **kwargs)  # 使用当前轴  
    ax.axhline(y=99, color='red', linestyle='--', label='99% Threshold')  # 添加红线  
    ax.legend()  # 显示图例  

# 使用 map 方法绘制 Identity曲线并添加阈值线  
g.map(draw_plot, "Window Start", "Mean Identity")  

# 添加标题和标签  
#g.fig.suptitle("Identity Distribution across genome", fontsize=16)  
g.set_axis_labels("Genomic Position (Mbp)", "Mean Identity (%)")  
g.set_titles(col_template="{col_name}")  

# 调整布局  
plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为总标题留出空间  

# 保存为 PDF 文件  
plt.savefig("identity_distribution_facet_with_step.pdf", dpi=300)  

# 显示图像  
plt.show()  
