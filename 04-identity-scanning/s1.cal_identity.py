from Bio import SeqIO  
from collections import defaultdict  
import pandas as pd  
#usage: python s1.cal_identity.py
# required input 1:core_aln_file = "C1237_concate_coreLCB.fasta"  # 核心基因组对齐文件
# required input 2:pop_file = "gnm_pop"  #基因组与clade的映射文件 
def load_population_file(pop_file):  
    """加载种群信息文件，返回{基因组: 种群}字典和种群"""  
    genome_to_pop = {}  
    populations = set()  
    with open(pop_file, 'r') as f:  
        for line in f:  
            genome, pop = line.strip().split()  
            genome_to_pop[genome] = pop  
            populations.add(pop)  
    return genome_to_pop, sorted(populations)  

def calculate_identity(sequences, window_size=1000, step_size=1000):  
    """计算每个窗口的种群间identity"""  
    identities = []  
    alignment_length = len(next(iter(sequences.values())))  

    # 滑动窗口处理  
    for start in range(0, alignment_length, step_size):  
        end = min(start + window_size, alignment_length)  

        # 提取窗口内序列  
        window_seqs = defaultdict(list)  # 使用默认字典  
        for genome, seq in sequences.items():  
            pop = genome_to_pop.get(genome)  # 使用 get 方法安全访问  
            if pop is None:  
                print(f"Warning: Genome '{genome}' not found in population file. Skipping.")  
                continue  # 如果种群未找到，则跳过  
            window_seqs[pop].append(seq[start:end])  

        # 计算不同种群组合的identity  
        for popA in window_seqs:  
            if len(window_seqs[popA]) > 1:  
                for popB in window_seqs:  
                    if popA != popB and len(window_seqs[popB]) > 1:  
                        identity_sum = 0  
                        count = 0  

                        # 逐对比对种群成员  
                        for seqA in window_seqs[popA]:  
                            for seqB in window_seqs[popB]:  
                                matching_bases = sum(1 for a, b in zip(seqA, seqB) if a == b and a != '-')  # 匹配位点数  
                                total_bases = sum(1 for a in zip(seqA, seqB) if a[0] != '-' and a[1] != '-')  # 符合率计算总碱基数  
                                if total_bases > 0:  # 排除总碱基为0的情况  
                                    identity = (matching_bases / total_bases) * 100  
                                    identity_sum += identity  
                                    count += 1  

                        # 计算平均identity  
                        if count > 0:  
                            mean_identity = identity_sum / count  
                            identities.append((popA, popB, start, end, mean_identity))  

    return identities  

if __name__ == "__main__":  
    core_aln_file = "C1237_concate_coreLCB.fasta"  # 核心基因组对齐文件
    pop_file = "gnm_pop"               # 种群信息文件  
    window_size = 1000                  # 窗口大小1 Kbp  
    step_size = 500                     # 步长500bp，实现50%重叠

    # 加载数据  
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(core_aln_file, "fasta")}  
    genome_to_pop, populations = load_population_file(pop_file)  

    # 计算identity值  
    identity_results = calculate_identity(sequences, window_size, step_size)  

    # 输出结果  
    identity_df = pd.DataFrame(identity_results, columns=["Population A", "Population B", "Window Start", "Window End", "Mean Identity"])  
    identity_df.to_csv("pop_identity_results.csv", index=False)  

    print("Identity calculation complete. Results saved to pop_identity_results.csv.")  
