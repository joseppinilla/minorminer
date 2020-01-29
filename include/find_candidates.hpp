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
                   chainmap &candidates) {

        int num_vars = var_g.num_nodes();

        // Fill up bins with qubits with the same location. O(N)
        binmap binning;
        locpair loc;

        int columns=0, rows=0;
        for (int u = 0; u < qubit_g.num_nodes(); u++) {
            loc = Tloc[u];
            int x = std::get<0>(loc);
            int y = std::get<1>(loc);

            if (x > columns) columns = x;
            if (y > rows) rows = y;
            binning[loc].push_back(u);
        }

        // Find scale of input graph locaions
        float xmin=0.0, ymin=0.0;
        float xmax=0.0, ymax=0.0;
        for (int u = 0; u < var_g.num_nodes(); u++) {
            loc = Sloc[u];
            float x = std::get<0>(loc);
            float y = std::get<1>(loc);

            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
        }

        // Candidates are qubits in bin that matches
        vector<int> qubits;
        for (int u = 0; u < var_g.num_nodes(); u++) {
            loc = Sloc[u];
            float x = std::get<0>(loc);
            float y = std::get<1>(loc);

            int binx = (int)(x-xmin)/(xmax-xmin)*(columns);
            int biny = (int)(y-ymin)/(ymax-ymin)*(rows);
            qubits = binning[pair<int,int>(binx,biny)];
            for (int q : qubits) {
              candidates[u].push_back(q);
            }
        }
}
}
