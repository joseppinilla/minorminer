""" Minorminer is not very good at post-embedding optimization because it
    searches from pre-placed neighbours. """

import minorminer
import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt

from dwave_networkx.drawing.distinguishable_colors import distinguishable_color_map as cmp

# Create problem graph
P7 = nx.Graph([(1, 2), (2, 3), (3, 4),
                (4, 5), (5, 6), (6, 7),
                (7, 8), (8, 1)])

# Use same colormap for problem graph and embedding
color = cmp(len(P7))
chain_color = [color(i-1) for i in P7]
plt.figure(0)
nx.draw(P7, node_color=chain_color, with_labels=True)

# Create target Chimera graph
C = dnx.chimera_graph(16, 16, coordinates=False)
# Miner object
mn = minorminer.miner(P7,C,verbose=4)

# Initialize with spread-out embedding
fixed_chains = {1: [0], 3: [123], 5: [2047], 7: [1927]}
bad_embedding = minorminer.find_embedding(P7, C, fixed_chains=fixed_chains)
plt.figure(1)
dnx.draw_chimera_embedding(C, bad_embedding, node_size=5)
print(mn.quality_key(bad_embedding))

# Improve using miner. Attempts at least twice (2) without improvement
embedding = mn.improve_embeddings([bad_embedding]).pop()
plt.figure(2)
dnx.draw_chimera_embedding(C, embedding, node_size=5)
print(mn.quality_key(embedding))

# Improve by initializing minorminer (more strict/attempts?)
# Attempts at least ten (10) times without improvement
embedding = minorminer.find_embedding(P7,C,initial_chains=bad_embedding,verbose=4)
plt.figure(3)
dnx.draw_chimera_embedding(C, embedding, node_size=5)
print(mn.quality_key(embedding))
