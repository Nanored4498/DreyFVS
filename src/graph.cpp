#include "graph.h"
#include "common.h"

#include <cassert>
#include <sstream>

using namespace std;

/////////////////
// Graph Edition
/////////////////

void Graph::add_edge(int u, int v) {
	if(!adj_out[u].insert(v)) return;
	adj_in[v].insert(u);
	++ m;
}

void Graph::remove_edge(int u, int v) {
	PROFIL_FUNC("Graph RM edge");
	adj_out[u].erase(v);
	adj_in [v].erase(u);
	-- m;
	PROFIL_RET;
}

void Graph::remove_vertex(int u) {
	PROFIL_FUNC("Graph RM vertex");
	-- n;
	assert(adj_in[u].empty() && adj_out[u].empty());
	if(u < n) {
		adj_in [u] = move(adj_in .back());
		adj_out[u] = move(adj_out.back());
		index[u] = index.back();
		inv_index->at(index[u]) = u;
		for(int v : adj_in [u]) { assert(adj_out[v].count(n)); adj_out[v].erase(n); adj_out[v].insert(u); }
		for(int v : adj_out[u]) { assert(adj_in [v].count(n)); adj_in [v].erase(n); adj_in [v].insert(u); }
	}
	adj_in .pop_back();
	adj_out.pop_back();
	index.pop_back();
	PROFIL_RET;
}

///////////
// INPUT
///////////

Graph Graph::from_istream(istream &is) {
	int n, m, t, i = -1;
	Graph g;
	string l;
	while(getline(is, l)) {
		if(!l.empty() && l[0] == '%') continue;
		istringstream iss(l);
		if(i == -1) {
			iss >> n >> m >> t;
			assert(t == 0);
			shared_ptr<Vi> inv_index = make_shared<Vi>(n);
			iota(inv_index->begin(), inv_index->end(), 0);
			g = move(Graph(n, inv_index));
		} else while(iss >> t) g.add_edge(i, t-1);
		++i;
	}
	assert(g.m == m);
	return g;
}

void Graph::print(ostream &os) const {
	os << n << ' ' << m << " 0\n";
	for(const Si &a : adj_out) {
		for(int u : a) os << u+1 << ' ';
		os << '\n';
	}
}
