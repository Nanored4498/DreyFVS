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
	adj_in [v].insert(u);
	++ m;
	moves.emplace_back(Move::Op::ADD_EDGE, u, v);
}

void Graph::remove_edge(int u, int v) {
	PROFIL_FUNC("Graph RM edge");
	adj_out[u].erase(v);
	adj_in [v].erase(u);
	-- m;
	moves.emplace_back(Move::Op::REMOVE_EDGE, u, v);
	PROFIL_RET;
}

void Graph::remove_vertex(int u) {
	PROFIL_FUNC("Graph RM vertex");
	-- n;
	assert(adj_in[u].empty() && adj_out[u].empty());
	moves.emplace_back(Move::Op::REMOVE_VERTEX, u, index[u]);
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

void Graph::undo(const std::pair<int, int> &time) {
	PROFIL_FUNC("Graph undo");
	while((int) moves.size() > time.first) {
		const Move &mv = moves.back();
		switch(mv.op) {
		case Move::Op::ADD_EDGE:
			adj_out[mv.u].erase(mv.v);
			adj_in [mv.v].erase(mv.u);
			-- m;
			break;
		case Move::Op::REMOVE_EDGE:
			adj_out[mv.u].insert(mv.v);
			adj_in [mv.v].insert(mv.u);
			++ m;
			break;
		case Move::Op::REMOVE_VERTEX:
			if(mv.u == n) {
				adj_in .emplace_back();
				adj_out.emplace_back();
				index.push_back(mv.ind);
			} else {
				adj_in .emplace_back(move(adj_in [mv.u]));
				adj_out.emplace_back(move(adj_out[mv.u]));
				index.push_back(index[mv.u]);
				inv_index->at(index[mv.u]) = n;
				assert(adj_in[mv.u].empty() && adj_out[mv.u].empty());
				index[mv.u] = mv.ind;
				for(int v : adj_in [n]) { assert(adj_out[v].count(mv.u)); adj_out[v].erase(mv.u); adj_out[v].insert(n); }
				for(int v : adj_out[n]) { assert(adj_in [v].count(mv.u)); adj_in [v].erase(mv.u); adj_in [v].insert(n); }
			}
			inv_index->at(mv.ind) = mv.u;
			++ n;
			break;
		default:
			assert(false);
		}
		moves.pop_back();
	}
	solution.resize(time.second);
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
	g.clearMoves();
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