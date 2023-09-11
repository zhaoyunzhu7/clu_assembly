import networkx as nx
from Bio import SeqIO
from Levenshtein import distance
from Bio import pairwise2

# 基于编辑距离（Levenshtein距离）将相似的序列聚类成图，并找到每个图中连接边数最多的节点。
# 读取eDNAs_all.fa文件并将序列存储到一个列表中
# 读取FASTA文件中的序列
def read_fasta(filename):
    sequences = []
    with open(filename, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
    return sequences

def build_graphs(sequences, threshold):
    graphs = []
    for i, seq1 in enumerate(sequences):
        graph_added = False
        for graph in graphs:
            min_dist = float('inf')
            min_dist_node = None
            for node in graph.nodes:
                dist = distance(seq1, sequences[node])
                if dist < min_dist:
                    min_dist = dist
                    min_dist_node = node
            if min_dist <= threshold:
                graph.add_node(i, sequence=seq1)
                graph.add_edge(i, min_dist_node, weight=min_dist)
                graph_added = True
                break

        if not graph_added:
            graph = nx.Graph()
            graph.add_node(i, sequence=seq1)
            graphs.append(graph)

    return graphs

# 合并每个图中的序列并输出
def compute_merged_sequence(graphs):
    merged_add = []
    for idx, graph in enumerate(graphs):#idx是graph的索引
        max_degree_node = max(graph.degree,key=lambda x: x[1])  # key=lambda x: x[1] 帮助 max() 函数在 graph.degree 字典中找到具有最大度数的节点。

        # 获取最大度数节点的序列
        max_degree_sequence = graph.nodes[max_degree_node[0]]['sequence']

        # 合并相邻节点的序列
        merged_sequence = max_degree_sequence
        for neighbor in graph.neighbors(max_degree_node[0]):
            neighbor_sequence = graph.nodes[neighbor]['sequence']
            alignments = pairwise2.align.globalxx(merged_sequence, neighbor_sequence)
            best_alignment = alignments[0]
            aligned_seq1 = best_alignment.seqA
            aligned_seq2 = best_alignment.seqB
            merged_sequence = "".join([aligned_seq1[i] if aligned_seq1[i] != "-" else aligned_seq2[i] for i in range(len(aligned_seq1))])
        # merged_add.append(merged_sequence)
        # merged_add.append(f"graph{idx+1}:{merged_sequence}")
        merged_add.append(merged_sequence)
    return merged_add


# filename = "eDNAs_all.fa"
# # 设置编辑距离阈值,<=5的记作相似
# threshold = 5
# sequences = read_fasta(filename)
# graphs = build_graphs(sequences, threshold)
# s=compute_merged_sequence(graphs)
# print(s)

