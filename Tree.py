import numpy as np
import math
from ete3 import Tree

def compare(string1, string2, no_match_c=' ', match_c='|'):
    """
    Calculate the number of differences between two sequences.
    """
    if len(string2) < len(string1):
        string1, string2 = string2, string1
    n_diff = sum(c1 != c2 for c1, c2 in zip(string1, string2))
    n_diff += len(string2) - len(string1)  # Account for length difference
    return n_diff

def load_fasta(file_path):
    """
    Load sequences from a FASTA file into a dictionary.
    """
    sequences = {}
    with open(file_path, 'r') as file:
        seqname = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                seqname = line[1:]
                sequences[seqname] = []
            else:
                sequences[seqname].append(line)
    fasta = {k: ''.join(v) for k, v in sequences.items()}
    return fasta

def global_alignment(x, y, match_score=1, mismatch_score=-1, gap=-2):
    """
    Perform global alignment using Needleman-Wunsch algorithm.
    """
    n_matrix = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            n_matrix[i, j] = match_score if x[i] == y[j] else mismatch_score

    main_matrix = np.zeros((len(x) + 1, len(y) + 1))
    for i in range(len(x) + 1):
        main_matrix[i, 0] = gap * i
    for j in range(len(y) + 1):
        main_matrix[0, j] = gap * j

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            left = main_matrix[i, j-1] + gap
            up = main_matrix[i-1, j] + gap
            diagonal = main_matrix[i-1, j-1] + n_matrix[i-1, j-1]
            main_matrix[i, j] = max(left, diagonal, up)

    align1, align2 = [], []
    xi, yj = len(x), len(y)
    while xi > 0 or yj > 0:
        if xi > 0 and yj > 0 and main_matrix[xi, yj] == main_matrix[xi-1, yj-1] + n_matrix[xi-1, yj-1]:
            align1.append(x[xi-1])
            align2.append(y[yj-1])
            xi -= 1
            yj -= 1
        elif yj > 0 and main_matrix[xi, yj] == main_matrix[xi, yj-1] + gap:
            align2.append(y[yj-1])
            align1.append('-')
            yj -= 1
        else:
            align1.append(x[xi-1])
            align2.append('-')
            xi -= 1

    return align1[::-1], align2[::-1]

def calculate_p_distance(fasta, alignment_function, comparison_function):
    """
    Calculate p-distance between sequences in the fasta dictionary.
    """
    distancedict = {}
    processed_pairs = set()  # To keep track of processed pairs
    for key1, seq1 in fasta.items():
        for key2, seq2 in fasta.items():
            if key1 != key2 and (key1, key2) not in processed_pairs and (key2, key1) not in processed_pairs:
                align1, align2 = alignment_function(seq1, seq2)
                string1, string2 = ''.join(align1), ''.join(align2)
                num_difference = comparison_function(string1, string2)
                pdistance = num_difference / len(string1)
                distancedict[f"{key1}-{key2}"] = pdistance
                processed_pairs.add((key1, key2))
    return distancedict

def juke_and_cantor_model(distancedict):
    """
    Apply Jukes-Cantor model to convert p-distance to evolutionary distance.
    """
    model_dict = {}
    for key, pvalue in distancedict.items():
        try:
            if pvalue < 0.75:  # Prevent math domain errors
                jcvalue = -3/4 * math.log(1 - 4/3 * pvalue)
            else:
                jcvalue = float('inf')  # Represents very high evolutionary distance
            model_dict[key] = jcvalue
        except ValueError:  # Handle cases where math domain error occurs
            model_dict[key] = None
    return {k: v for k, v in model_dict.items() if v is not None}

def construct_pre_newick_format(model_dict):
    """
    Construct pre-Newick format using a hierarchical clustering approach.
    """
    # Create a priority queue based on the distances
    clusters = {k: [k] for k in model_dict.keys()}
    while len(clusters) > 1:
        # Find the pair with the smallest distance
        min_pair = min(model_dict, key=model_dict.get)
        min_distance = model_dict[min_pair]

        # Merge the clusters
        cluster1, cluster2 = min_pair.split('-')
        new_cluster = f"({clusters[cluster1][0]},{clusters[cluster2][0]}):{min_distance}"
        clusters[cluster1] = [new_cluster]
        del clusters[cluster2]

        # Update the distance matrix
        model_dict = {k: v for k, v in model_dict.items() if cluster1 not in k and cluster2 not in k}

    return clusters[min(clusters.keys())][0]

def generate_newick(tree_str):
    """
    Finalize the Newick format and render the tree.
    """
    final_tree = f"({tree_str});"
    t = Tree(final_tree, format=9)
    t.render('phylogenetic_Tree.png', w=183, units="mm")
    return t

# Main execution
fasta = load_fasta('msa.fasta')
distances = calculate_p_distance(fasta, global_alignment, compare)
model = juke_and_cantor_model(distances)
pre_newick_tree = construct_pre_newick_format(model)
final_tree = generate_newick(pre_newick_tree)
print(final_tree)
