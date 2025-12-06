# 读取input2并创建gnm到pop的映射字典  
gnm_to_pop = {}  
with open('gnm_pop', 'r') as f:  
    for line in f:  
        gnm, pop = line.strip().split()  
        gnm_to_pop[gnm] = pop  

# 遍历input1并提取符合条件的行  
result = []  
with open('Sulfito_ANI.txt', 'r') as f:  
    for line in f:  
        gnm1, gnm2, ani, length, aln = line.strip().split('\t')  
        gnm1 = gnm1.split('.')[0]
        gnm2 = gnm2.split('.')[0]
        
        # 检查gnm1和gnm2是否分别属于pop1和pop2  
        if gnm_to_pop.get(gnm1) == 'C9' and gnm_to_pop.get(gnm2) == 'C1':  
            result.append(line)  
        # 也可以反过来，如果gnm1是pop2而gnm2是pop1  
        elif gnm_to_pop.get(gnm1) == 'C1' and gnm_to_pop.get(gnm2) == 'C9':  
            result.append(line)  

# 保存结果  
with open('output', 'w') as f:  
    f.writelines(result)  


# 读取输出文件  
ani_values = []  
with open('output', 'r') as f:  
    for line in f:  
        # 跳过空行  
        if not line.strip():  
            continue  
        columns = line.strip().split('\t')  
        # 提取ANI值（第3列）  
        if len(columns) >= 3:  
            ani = columns[2]  
            try:  
                ani_float = float(ani)  
                ani_values.append(ani_float)  
            except ValueError:  
                pass  # 跳过非数字值  

# 计算平均值  
if ani_values:  
    average_ani = sum(ani_values) / len(ani_values)  
    print(f"ANI values average: {average_ani:.4f}")  
else:  
    print("No valid ANI values found.")  
