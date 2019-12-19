#pragma once

#include <iostream>

#include "util.hpp"
#include "graph.hpp"

namespace find_embedding {

typedef pair<int,int> locpair;
typedef map<int,locpair> locmap;
typedef map<int,vector<int>> chainmap;
typedef map<locpair,vector<int>> binmap;

/*
This is an example using location information from the input and
target graphs to create an initial mapping, aka candidates.

TODO: The real implementation would have:
  locpair as <float,float>
  Sloc would be transformed to optimize bin assignments

*/

void findCandidates(graph::input_graph &var_g, graph::input_graph &qubit_g,
                   locmap &Sloc, locmap &Tloc,
                   chainmap &candidates) {

        int num_vars = var_g.num_nodes();
        // Fill up bins with qubits with the same location. O(N)
        binmap bins;
        locpair loc;
        for (int u = 0; u < qubit_g.num_nodes(); u++) {
            loc = Tloc[u];
            bins[loc].push_back(u);
        }

        // Candidates are qubits in bin that matches
        vector<int> qubits;
        for (int u = 0; u < var_g.num_nodes(); u++) {
            loc = Sloc[u];
            qubits = bins[loc];
            for (int q : qubits) {
              candidates[u].push_back(q);
            }
        }
}
}
