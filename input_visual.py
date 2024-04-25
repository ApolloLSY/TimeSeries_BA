import pandas as pd
import numpy as np
import pandas as pd
import numpy as np

models = {
    'age_Horvath': {'params': pd.DataFrame({'cpg': ['cpg1', 'cpg2', 'cpg3'], 'weight': [0.5, 0, -0.4]}), 'intercept': 0.2},
    'age_Hannum': {'params': pd.DataFrame({'cpg': ['cpg1', 'cpg3', 'cpg4'], 'weight': [0.3, -0.2, 0]}), 'intercept': 0.1},
    'age_Levine': {'params': pd.DataFrame({'cpg': ['cpg2', 'cpg4', 'cpg5'], 'weight': [-0.5, 0.5, 0]}), 'intercept': 0.3},
    'age_Lu': {'params': pd.DataFrame({'cpg': ['cpg1', 'cpg5', 'cpg6'], 'weight': [0, 0.4, -0.3]}), 'intercept': 0.1},
    'age_Ying_adapt': {'params': pd.DataFrame({'cpg': ['cpg2', 'cpg3', 'cpg6'], 'weight': [0.2, -0.1, 0]}), 'intercept': 0.2},
    'age_Ying_damage': {'params': pd.DataFrame({'cpg': ['cpg4', 'cpg5', 'cpg7'], 'weight': [0, -0.5, 0.4]}), 'intercept': 0.4}
}

# 假设df是你的数据集
# 为每个模型计算年龄并更新subset
subsets = {}
for model_name, model_data in models.items():
    params = model_data['params']
    intercept = model_data['intercept']
    
    # 计算年龄
    df[model_name] = predict_age(df, params, intercept).astype(int)
    
    # 过滤非零权重的cpg列
    relevant_columns = params[params['weight'] != 0]['cpg']
    subset = df[relevant_columns.tolist() + [model_name]].copy()
    subset = subset.sort_values(by=model_name)
    subset = subset.groupby(model_name).mean()
    subsets[model_name] = subset

    # 打印描述性统计
    print(f"Descriptive statistics for {model_name}:")
    print(df[model_name].describe())

# 打印每个subset的信息
for name, subset in subsets.items():




import matplotlib.pyplot as plt
import seaborn as sns

# 假设 subset_age_Horvath 是其中一个已经处理好的数据集
# 时间序列图
plt.figure(figsize=(12, 6))
for cpg in subset_age_Horvath.columns[:-1]:  # 不包括年龄列
    plt.plot(subset_age_Horvath.index, subset_age_Horvath[cpg], label=cpg)
plt.title('CpG Expression Over Age for Horvath Model')
plt.xlabel('Age')
plt.ylabel('Expression')
plt.legend()
plt.show()

# 热图
plt.figure(figsize=(10, 8))
sns.heatmap(subset_age_Horvath.T, cmap='viridis', annot=True)  # 转置以使每一行表示一个CpG位点，每一列表示一个年龄
plt.title('Heatmap of CpG Expression Over Age for Horvath Model')
plt.xlabel('Age')
plt.ylabel('CpG Sites')
plt.show()

    print(f"Updated {name} subset:")
    print(subset.head())



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# 创建3D图形
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 假设 'age', 'cpg_site', 和 'expression' 分别是年龄，CpG位点索引，和表达量数据
age = np.random.randint(20, 70, size=100)  # 示例年龄数据
cpg_site = np.random.randint(0, 100, size=100)  # 示例CpG位点数据
expression = np.random.rand(100) * 100  # 示例表达量数据

# 绘制散点图
sc = ax.scatter(age, cpg_site, expression, c=expression, cmap='viridis')

# 添加标签和标题
ax.set_xlabel('Age')
ax.set_ylabel('CpG Site Index')
ax.set_zlabel('Expression Level')
plt.colorbar(sc)  # 添加颜色条
plt.title('3D Scatter Plot of CpG Expression')

# 显示图形
plt.show()

# 如果做3D图就要用全套cpg位点且要查一下cpg位点在基因上的毗邻关系（Illumina 450K 每个CpG位点在人类基因组中的染色体和坐标） 否则位点不连续
