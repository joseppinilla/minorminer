import time
import minorminer
import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 4, subplot_kw={'aspect':'equal'})

# Target graph
C = dnx.chimera_graph(16, 16, coordinates=True)
nx.draw(C,pos=dnx.chimera_layout(C),node_size=5,ax=axs.flat[0])

# Trivial problem graph
G = nx.grid_2d_graph(16,16)
# Non-exact layout for grid i.e. has overlap
nx.draw(G,pos=nx.spring_layout(G,seed=16),node_size=5,ax=axs.flat[1])
# Exact layout
nx.draw(G,pos={v:v for v in G},node_size=5,ax=axs.flat[2])

# For quality metrics
miner = minorminer.miner(G,C)

########################### Minorminer
start = time.time()
embedding = minorminer.find_embedding(G, C)
end = time.time()
dnx.draw_chimera_embedding(C, embedding, G, node_size=5, ax=axs.flat[4])

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
print("Layout-Agnostic: Total: %d Max: %d t: %.2fs" % (total,L[0],end - start))

########################### Minorminer.layout
start = time.time()
# Measure time of finding layout as part of algorithm
source_layout = {v:v for v in G}
# source_layout = nx.spring_layout(G,seed=16)
# NOTE: target_layout must be (int,int) for now. Using Chimera tile.
target_layout = {(i,j,u,k):(i,j) for (i,j,u,k) in C}
embedding = minorminer.layout.find_embedding(G,C,(source_layout,target_layout))
end = time.time()
dnx.draw_chimera_embedding(C,embedding,G,node_size=5,ax=axs.flat[5])

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
print("Layout-Aware (MML): Total: %d Max: %d t: %.2fs" % (total,L[0],end - start))


########################### Topominer w/ non-exact layout
start = time.time()
# Measure time of finding layout as part of algorithm
source_layout = nx.spring_layout(G,seed=16)
# NOTE: target_layout must be (int,int) for now. Using Chimera tile.
target_layout = {(i,j,u,k):(i,j) for (i,j,u,k) in C}
embedding = minorminer.topo_embedding(G,C,None,target_layout)
end = time.time()
dnx.draw_chimera_embedding(C,embedding,G,node_size=5,ax=axs.flat[6])

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
print("Layout-Aware Spring: Total: %d Max: %d t: %.2fs" % (total,L[0],end - start))

########################### Topominer w/ exact layout
start = time.time()
source_layout = {v:v for v in G}
target_layout = {(i,j,u,k):(i,j) for (i,j,u,k) in C}
embedding = minorminer.topo_embedding(G,C,source_layout,target_layout)
end = time.time()
dnx.draw_chimera_embedding(C,embedding,G,node_size=5,ax=axs.flat[7])

state, O, L = miner.quality_key(embedding,embedded=True)
total = sum([L[i]*L[i+1] for i in range(0,len(L),2)])
print("Layout-Aware Exact: Total: %d Max: %d t: %.2fs" % (total,L[0],end - start))

plt.tight_layout()
plt.show()
