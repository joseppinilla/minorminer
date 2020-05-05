import minorminer
import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2, subplot_kw={'aspect':'equal'})

# Trivial problem graph
G = nx.grid_2d_graph(16,16)
source_bins = {v:v for v in G}
nx.draw(G, pos=source_bins, node_size=10, ax=axs.flat[0])

# Target graph
C = dnx.chimera_graph(16, 16, coordinates=True)
target_bins = {(i,j,u,k):(i,j) for (i,j,u,k) in C}
dnx.draw_chimera(C,node_size=10,ax=axs.flat[1])

import time
# For Quality Key metrics
miner = minorminer.miner(G,C)

# Minorminer
start = time.time()
embedding = minorminer.find_embedding(G, C)
end = time.time()
dnx.draw_chimera_embedding(C, embedding, G, node_size=1, ax=axs.flat[2])

state, O, L = miner.quality_key(embedding,embedded=True)

total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
print(f"Total: {total} Max: {L[0]} t: {end - start:.2f}s")

# Topominer
# TEMP: Prototype version of topominer maps locations directly. Sloc->Tloc.
start = time.time()
embedding = minorminer.topo_embedding(G,C,source_bins,target_bins)
end = time.time()
dnx.draw_chimera_embedding(C, embedding, G, node_size=1, ax=axs.flat[3])

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
print(f"Total: {total} Max: {L[0]} t: {end - start:.2f}s")
