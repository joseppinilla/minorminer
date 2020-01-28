#pragma once

#include <iostream>

#include "util.hpp"
#include "graph.hpp"

namespace find_embedding {

typedef pair<int,int> intpair;
typedef pair<float,float> locpair;
typedef map<int,locpair> locmap;
typedef map<int,vector<int>> chainmap;
typedef map<locpair,vector<int>> binmap;

/*
This is an example using location information from the input and
target graphs to create an initial mapping, aka candidates.

TODO:
  Sloc should be transformed to optimize bin assignments

*/

void findCandidates(graph::input_graph &var_g, graph::input_graph &qubit_g,
                   locmap &Sloc, locmap &Tloc,
                   intpair bins,
                   chainmap &candidates) {

        int num_vars = var_g.num_nodes();
        // Fill up bins with qubits with the same location. O(N)
        binmap binning;
        locpair loc;
        int x_bins = std::get<0>(bins);
        int y_bins = std::get<1>(bins);

        for (int u = 0; u < qubit_g.num_nodes(); u++) {
            loc = Tloc[u];
            float x = std::get<0>(loc)*x_bins;
            float y = std::get<1>(loc)*y_bins;
            binning[loc].push_back(u);
        }

        // Candidates are qubits in bin that matches
        vector<int> qubits;
        for (int u = 0; u < var_g.num_nodes(); u++) {
            loc = Sloc[u];

            float x = std::get<0>(loc)*x_bins;
            float y = std::get<1>(loc)*y_bins;
            qubits = binning[loc];
            for (int q : qubits) {
              candidates[u].push_back(q);
            }
        }
}
}
