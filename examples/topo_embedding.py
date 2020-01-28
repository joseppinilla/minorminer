import minorminer
import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2, subplot_kw={'aspect':'equal'})


# Trivial problem graph
plt.subplot(2, 2, 1)
G = nx.grid_2d_graph(16,16)
source_bins = {v:v for v in G}
source_layout = {k:(x/16,y/16) for k, (x,y) in source_bins.items()}
nx.draw(G, pos=source_layout, node_size=1)

# Non-trivial problem graph
# plt.subplot(2, 2, 1)
# G = nx.circular_ladder_graph(50)
# source_bins = nx.spring_layout(G)
# source_layout = {k:((x+1)/2,(y+1)/2) for k, (x,y) in source_bins.items()}
# nx.draw(G, pos=source_layout, node_size=1)

# Target graph
plt.subplot(2, 2, 2)
C = dnx.chimera_graph(16, 16, coordinates=True)
target_bins = {(i,j,u,k):(i,j) for (i,j,u,k) in C}
target_layout = {k:(x/16,y/16) for k, (x,y) in target_bins.items()}
dnx.draw_chimera(C, node_size=1)

import time
# For Quality Key metrics
miner = minorminer.miner(G,C)

# Minorminer
plt.subplot(2, 2, 3)
start = time.time()
embedding = minorminer.find_embedding(G, C)
end = time.time()
dnx.draw_chimera_embedding(C, embedding, G, node_size=1)

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
plt.title(f"Total: {total} Max: {L[0]}\n t: {end - start:.2f}s")

# Topominer
# TEMP: Prototype version of topominer maps locations directly. Sloc->Tloc.
plt.subplot(2, 2, 4)
start = time.time()
embedding = minorminer.topo_embedding(G,C,source_layout,target_layout,bins=(4,4))
end = time.time()
dnx.draw_chimera_embedding(C, embedding, G, node_size=1)

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
plt.title(f"Total: {total} Max: {L[0]}\n t: {end - start:.2f}s")

plt.tight_layout()
plt.show()
