
import random

def compute_max_overlap_naive(u, v):
    """Compute the maximum overlap between the suffix of u and the prefix of v."""
    max_overlap = 0
    u_len, v_len = len(u), len(v)
    # Check all possible suffix-prefix alignments
    for i in range(1, min(u_len, v_len) + 1):
        if u[-i:] == v[:i]:
            max_overlap = i
    return max_overlap

def compute_max_overlap(u, v):
    """Compute the maximum overlap between the suffix of u and the prefix of v using KMP failure function."""
    string = v + '#' + u  # Concatenate v, a separator, and u
    lps = [0] * len(string)
    for i in range(1, len(string)):
        length = lps[i - 1]
        while length > 0 and string[i] != string[length]:
            length = lps[length - 1]
        if string[i] == string[length]:
            length += 1
        lps[i] = length
    return lps[-1]  # The last value of lps gives the maximum overlap


def build_overlap_graph(fragments, overlap_function):
    # Initialize the graph
    graph = {}
    n = len(fragments)

    # Compute overlaps for all pairs of fragments
    for i in range(n):
        for j in range(n):
            if i != j:  # No self-loop
                u, v = fragments[i], fragments[j]
                graph[(u, v)] = overlap_function(u, v)

    return graph

def generate_fragments(num_fragments, fragment_length):
    bases = ['A', 'T', 'C', 'G']
    fragments = [''.join(random.choices(bases, k=fragment_length)) for _ in range(num_fragments)]
    return fragments
