#pragma once

#include "graph.h"

namespace Kernel {

struct SCC {
	std::vector<Vi> cs;
	Vi c;
	SCC(const Graph &g, bool use_pi=true);
	inline void clear() { cs.clear(); c.clear(); };
};

// Remove a vertex i and apply reductions
void removeVertex(Graph &g, int i, bool inSol);

// Remove vertices whose in/out adjacency is a clique
bool cliqueReduction(Graph &g, int u);
void cliqueReduction(Graph &g);

// Pi-edge reduction
bool PIedge(Graph &g);
// Domination reduction
bool domination(Graph &g);

// remove vertices having a single petal
void remove_single_petal(Graph &g);

// Apply deg <= 1 reduction and clique reduction
void simplify0(Graph &g);
// Plus Pi-edge reduction
void simplify1(Graph &g);

// Apply deg <= 1 reduction and clique reduction then split in SCC or try to remove single petals
// Return a couple (vector of SCCs, solution -- vector of vertices in the solution)
std::pair<std::vector<Graph>, Vi> simplify(Graph &g);

}