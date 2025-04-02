#pragma once

#include "graph.h"

namespace Kernel {

struct SCC {
	std::vector<Vi> cs;
	Vi c;
	SCC(const Graph &g, bool use_pi=true);
	inline void clear() { cs.clear(); c.clear(); };
};

// Remove a vertex i and recursively apply reductions LOOP, IN0, IN1, OUT0, OUT1
void removeVertex(Graph &g, int i, bool inSol);

// Apply Core reduction
// Remove vertices whose in/out adjacency is a d-clique
void core(Graph &g);

// Apply all reduction rules
void reduce(Graph &g);

// Apply all reduction rules then split in SCC
// Return a couple (vector of SCCs, solution -- vector of vertices in the solution)
std::pair<std::vector<Graph>, Vi> reduce_and_split(Graph &g);

}