""" Minorminer is not very good at post-embedding optimization because it
    searches from pre-placed neighbours. """

import minorminer
import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt

# Create problem graph
P7 = nx.Graph([(1, 2), (2, 3), (3, 4),
                (4, 5), (5, 6), (6, 7),
                (7, 8), (8, 1)])
plt.subplot(2, 2, 1)
nx.draw(P7, with_labels=True)

# Create target Chimera graph
C = dnx.chimera_graph(16, 16, coordinates=False)

# Miner object
mn = minorminer.miner(P7,C)

# Initialize with spread-out embedding
fixed_chains = {1: [0], 3: [123], 5: [2047], 7: [1927]}
bad_embedding = minorminer.find_embedding(P7, C, fixed_chains=fixed_chains)
plt.subplot(2, 2, 2)
dnx.draw_chimera_embedding(C, bad_embedding, node_size=5)
print(mn.quality_key(bad_embedding))

# Improve using miner
embedding = mn.improve_embeddings([bad_embedding]).pop()
plt.subplot(2, 2, 3)
dnx.draw_chimera_embedding(C, embedding, node_size=5)
print(mn.quality_key(embedding))

# Improve by initializing minorminer (more strict/attempts?)
embedding = minorminer.find_embedding(P7,C,initial_chains=bad_embedding)
plt.subplot(2, 2, 4)
dnx.draw_chimera_embedding(C, embedding, node_size=5)
print(mn.quality_key(embedding))
