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


def calculate_outdegree(current_node, edges):
    # 初始化出度为0
    outdegree = 0

    # 遍历所有边，检查是否有以current_node为起点的边
    for edge in edges:
        if edge[0] == current_node:
            outdegree += 1

    return outdegree

def get_outgoing_neighbors(node, edges):
    # 查找通过出边连接到当前节点的所有节点
    outgoing_neighbors = [edge[1] for edge in edges if edge[0] == node]
    return outgoing_neighbors

def comparision_seq(ss_list, G):
    container = []
    container.append(ss_list[0])
    print(container)
    i = 0
    print(len(ss_list))
    while len(container) < len(ss_list):#56 56
        current_node = container[-1]
        # print(container,"container")
        # 计算当前节点的出度
        outdegree = calculate_outdegree(current_node, G.edges())
        i+=1
        # print("%%%%%%%%%%list_len",len(ss_list))
        # print(container,len(container),"第", i, "个节点", ss_list[i])
        if outdegree == 1:
            # print("11111111111111111111111111111111111111111111111111111")
            # 如果出度为1，则将当前节点的下一个节点作为container的下一个元素
            next_node = get_outgoing_neighbors(current_node, G.edges())[0]
            # 将下一个节点添加到container中，并将其设为当前节点
            container.append(next_node)
        else:
            # 如果出度不为1，则找出当前节点的出边所连节点中与ss_list[i]相同的节点
            # print("22222222222222222222222222222222222222222222222222222","current",current_node)

            next_nodes = get_outgoing_neighbors(current_node, G.edges())
            # print("不满足的next",next_nodes)
            #计算container下一个元素索引
            neighbor_matched = False
            for neighbor in next_nodes:
                # print("106行",i,neighbor,ss_list[i])
                if neighbor == ss_list[i]:
                    container.append(neighbor)
                    neighbor_matched = True
                    break
            if not neighbor_matched:
                container.append(ss_list[i])

    print(container)
    # 合并前后缀并还原为原始DNA序列
    assem_sequence = container[0]
    for segment in container[1:]:
        assem_sequence += segment[-1]
    print("组装后的序列：", assem_sequence)

    #ss_list
    ss_list_merge = ss_list[0]
    for i in ss_list[1:]:
        ss_list_merge += i[-1]
    print("聚类后的序列：",ss_list_merge)
    return container


if __name__ == "__main__":
    filename = "e.fa"
    threshold = 12
    sequences = read_fasta(filename)
    k = 3
    # k-mer sizec
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
        graphs_text = build_graphs(group, threshold)
        ss = compute_merged_sequence(graphs_text)
        # print(ss)
        ss_str = ss[0]
        # print("聚类后的序列：",ss_str)
        # for kmer, freq in dbg.items():
        #    print(f"Node: {kmer}, Frequency: {freq}")#freq为k-mer的频率

        # print("Edges:")
        # for (kmer1, kmer2), freq in edges.items():
        #     print(f"Edge: {kmer1} -> {kmer2}, Frequency: {freq}")#freq为边的频率
        # print("=" * 20)
        # 切分DNA序列并存储在列表中
        ss_list = [ss_str[i:i + k] for i in range(len(ss_str) - k + 1)]
        # print(ss_list)
        print(f"Group {idx + 1}:")
        dbg, edges = build_dbg(group, k)  # dbg是节点
        # Create a new directed graph
        G = nx.DiGraph()
        # Add nodes from dbg dictionary with their frequencies as node attributes
        for kmer, freq in dbg.items():
            G.add_node(kmer, frequency=freq)
        # Add edges from edges dictionary with their frequencies as edge attributes
        for (kmer1, kmer2), freq in edges.items():
            G.add_edge(kmer1, kmer2, frequency=freq)
        comparision_seq(ss_list, G)