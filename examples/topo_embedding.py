import minorminer
import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2, subplot_kw={'aspect':'equal'})

# Trivial problem graph
plt.subplot(2, 2, 1)
G = nx.grid_2d_graph(16,16)
source_layout = {v:v for v in G}
nx.draw(G, pos=source_layout, node_size=1)

# Target graph
plt.subplot(2, 2, 2)
C = dnx.chimera_graph(16, 16, coordinates=True)
target_layout = {(i,j,u,k):(i,j) for (i,j,u,k) in C}
dnx.draw_chimera(C, node_size=1)

# For Quality Key metrics
miner = minorminer.miner(G,C)

# Minorminer
plt.subplot(2, 2, 3)
embedding = minorminer.find_embedding(G, C)
dnx.draw_chimera_embedding(C, embedding, G, node_size=1)

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
plt.title(f"Total: {total} Max: {L[0]}\n")

# Topominer
# TEMP: Prototype version of topominer maps locations directly. Sloc->Tloc.
plt.subplot(2, 2, 4)
embedding = minorminer.topo_embedding(G, C, source_layout, target_layout)
dnx.draw_chimera_embedding(C, embedding, G, node_size=1)

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
plt.title(f"Total: {total} Max: {L[0]}\n")

plt.tight_layout()
plt.show()
