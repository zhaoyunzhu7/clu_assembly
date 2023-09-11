import networkx as nx
from Bio import SeqIO
from collections import defaultdict
import Levenshtein
from graph import compute_merged_sequence, build_graphs
import itertools


# 读取FASTA文件中的序列
def read_fasta(filename):
    sequences = []
    with open(filename, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
    return sequences


# 计算Levenshtein距离
def levenshtein_distance(seq1, seq2):
    return Levenshtein.distance(seq1, seq2)


# 将序列分组
def group_sequences(sequences, threshold):
    groups = []
    for seq in sequences:
        added = False
        for group in groups:
            for existing_seq in group:
                if levenshtein_distance(seq, existing_seq) <= threshold:
                    group.append(seq)
                    added = True
                    break
            if added:
                break
        if not added:
            groups.append([seq])
    return groups


# 构建DBG图
def build_dbg(sequences, k):
    dbg = defaultdict(int)
    edges = defaultdict(int)


    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            dbg[kmer] += 1
            if i < len(seq) - k:
                next_kmer = seq[i + 1:i + k + 1]
                edges[(kmer, next_kmer)] += 1

    # 去除频率小于阈值的节点
    filtered_dbg = {kmer: freq for kmer, freq in dbg.items() if freq >= 3}

    # 去除节点后，也要更新边的信息
    filtered_edges = {(kmer1, kmer2): freq for (kmer1, kmer2), freq in edges.items() if kmer1 in filtered_dbg and kmer2 in filtered_dbg}

    return filtered_dbg, filtered_edges




if __name__ == "__main__":
    filename = "eDNAs_all.fa"
    threshold = 5
    sequences = read_fasta(filename)
    k = 3
    # k-mer size
    #构建最小编辑距离图
    # graphs = build_graphs(sequences, threshold)
    #
    # #输出每个最小编辑距离图得到的序列
    # s = compute_merged_sequence(graphs)
    # print(s)

    #DBG划分类
    grouped_sequences = group_sequences(sequences, threshold)
    # graph_sequence =  build_graphs(grouped_sequences, threshold)
    # ss = compute_merged_sequence(graph_sequence)


    for idx, group in enumerate(grouped_sequences):
        print(f"Group {idx + 1} DBG:")
        dbg, edges = build_dbg(group, k)
        graphs_text = build_graphs(group, threshold)
        ss = compute_merged_sequence(graphs_text)
        print(f"ss:{ss}")
        # for kmer, freq in dbg.items():
        #     print(f"Node: {kmer}, Frequency: {freq}")#freq为k-mer的频率
        #
        # print("Edges:")
        # for (kmer1, kmer2), freq in edges.items():
        #     print(f"Edge: {kmer1} -> {kmer2}, Frequency: {freq}")#freq为边的频率
        # print("=" * 20)
        # Create a dictionary to store nodes with matching prefixes and their frequencies

        matching_prefixes = defaultdict(list)

        for (kmer1, kmer2), freq in edges.items():
            # 计算kmer1的最后(k-1)个碱基与kmer2的第一个(k-1)个碱基匹配的次数
            matching_count = sum(1 for i in range(1, k) if kmer1[-i:] == kmer2[:i])

            # 将kmer2追加到与匹配计数相关联的列表中
            matching_prefixes[matching_count].append(kmer2)

        # 根据出界度(有出界边的节点)组织节点
        nodes_with_outgoing_edges = set(kmer1 for (kmer1, kmer2) in edges.keys())

        for node in dbg.keys():
            if node not in nodes_with_outgoing_edges:
                # 这是一个out-degree = 0的起始节点
                print(f"Starting Node: {node}")
                continue

            # 计算节点的出度
            out_degree = sum(1 for (kmer1, kmer2) in edges.keys() if kmer1 == node)

            if out_degree == 1:
                # 如果out-degree为1，则使用匹配的前缀确定下一个节点
                matching_count = sum(1 for i in range(1, k) if node[-i:] == matching_prefixes[i][0][:i])

                # 检查是否有匹配前缀的节点
                if matching_count in matching_prefixes and len(matching_prefixes[matching_count]) > 0:
                    next_node = matching_prefixes[matching_count].pop(0)  # Remove and get the next node
                    print(f"Node: {node}, Out-Degree: {out_degree}, Next Node: {next_node}")
                else:
                    print(f"Node: {node}, Out-Degree: {out_degree}, No Next Node with Matching Prefixes")
            else:
                # 输出度> 1的节点(多选)
                print(f"Node: {node}, Out-Degree: {out_degree}, Multiple Out-Degree")

                # 确定ss中的索引以解决歧义
                current_position = 0
                current_kmer = None
                for idx, kmer in enumerate(matching_prefixes[matching_count]):
                    if ''.join(ss).startswith(kmer):
                        current_position = idx
                        current_kmer = kmer
                        break

                if current_kmer is not None:
                    next_node = current_kmer
                    print(f"Resolved Node: {next_node}")
                else:
                    print("Ambiguity in resolving the next node")

                # 将所选节点累加到路径列表中
                path = [node]
                while next_node:
                    path.append(next_node)
                    if current_position < len(matching_prefixes[matching_count]) - 1:
                        current_position += 1
                        next_node = matching_prefixes[matching_count][current_position]
                    else:
                        next_node = None
                path.append(node)
        # 将最终选择的路径打印为列表
        print("Final Path (List):", path)



