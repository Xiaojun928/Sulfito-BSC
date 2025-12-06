#!/usr/bin/env python3
import sys
from collections import defaultdict

def load_clade_info(clade_file):
    """加载clade信息，返回taxon到clade的映射字典"""
    taxon_to_clade = {}
    with open(clade_file, 'r') as f:
        next(f)  # 跳过标题行
        for line in f:
            taxon, clade = line.strip().split('\t')
            taxon_to_clade[taxon] = clade
    return taxon_to_clade

def standardize_clade_combination(clades):
    """标准化clade组合，返回排序后的元组"""
    return tuple(sorted(clades))

def standardize_topology(clade1, clade2, clade3, clade4):
    """标准化topology表示，返回标准化的字符串"""
    # 将clade对排序，确保ab|cd和ba|dc被视为相同
    pair1 = tuple(sorted([clade1, clade2]))
    pair2 = tuple(sorted([clade3, clade4]))
    # 将两个对排序，确保ab|cd和cd|ab被视为相同
    return tuple(sorted([pair1, pair2]))

def process_quartets(quartet_file, taxon_to_clade):
    """处理quartet数据，返回统计结果"""
    # 使用嵌套的defaultdict来存储统计结果
    # 外层key是标准化的clade组合，内层key是标准化的topology
    stats = defaultdict(lambda: defaultdict(int))
    
    with open(quartet_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 8:  # 确保行格式正确
                continue
                
            # 提取taxon和计数
            taxons = parts[1:5]
            counts = list(map(int, parts[5:8]))
            
            # 获取对应的clade
            clades = [taxon_to_clade[t] for t in taxons]
            
            # 标准化clade组合
            std_clades = standardize_clade_combination(clades)
            
            # 处理三种可能的topology
            topologies = [
                standardize_topology(clades[0], clades[1], clades[2], clades[3]),  # ab|cd
                standardize_topology(clades[0], clades[2], clades[1], clades[3]),  # ac|bd
                standardize_topology(clades[0], clades[3], clades[1], clades[2])   # ad|bc
            ]
            
            # 累加统计
            for topology, count in zip(topologies, counts):
                stats[std_clades][topology] += count
    
    return stats

def format_topology(topology):
    """将标准化的topology转换为可读的字符串"""
    pair1, pair2 = topology
    return f"{pair1[0]}-{pair1[1]}|{pair2[0]}-{pair2[1]}"

def main():
    if len(sys.argv) != 3:
        print("Usage: python process_quartets.py <clade_file> <quartet_file>")
        sys.exit(1)
        
    clade_file = sys.argv[1]
    quartet_file = sys.argv[2]
    
    # 加载clade信息
    taxon_to_clade = load_clade_info(clade_file)
    
    # 处理quartet数据
    stats = process_quartets(quartet_file, taxon_to_clade)
    
    # 输出结果
    print("QUARTET_TYPE\tTOPOLOGY\tCOUNT")
    for clade_comb, topology_counts in sorted(stats.items()):
        clade_str = "-".join(clade_comb)
        for topology, count in sorted(topology_counts.items()):
            topo_str = format_topology(topology)
            print(f"{clade_str}\t{topo_str}\t{count}")

if __name__ == "__main__":
    main() 
