#include "kernel.h"

#include "common.h"
#include <tuple>
#include <queue>
#include <cassert>

using namespace std;

namespace Kernel {

void removeVertex(Graph &g, int i, bool inSol) {
	static Vb isRemoved;
	static Vi stack, removed;
	if((int) isRemoved.size() < g.n) isRemoved.resize(g.n, false);
	const auto add2Stack = [&](int i, bool add) {
		isRemoved[i] = true;
		removed.push_back(i);
		stack.push_back(i);
		for(int j : g.adj_in [i]) if(j != i) g.adj_out[j].erase(i);
		for(int j : g.adj_out[i]) if(j != i) g.adj_in [j].erase(i);
		if(add) g.solution.push_back(g.index[i]);
		else for(int j : g.adj_in[i]) for(int k : g.adj_out[i]) g.add_edge(j, k);
	};
	add2Stack(i, inSol);
	do {
		dfs:
		const int i = stack.back();
		while(!g.adj_in[i].empty()) {
			const int j = *(--g.adj_in[i].end());
			g.remove_edge(j, i);
			if(!isRemoved[j]) {
				const bool self_loop = g.adj_out[j].count(j);
				if(self_loop || g.adj_in[j].size() <= 1 || g.adj_out[j].size() <= 1) {
					add2Stack(j, self_loop);
					goto dfs;
				}
			}
		}
		while(!g.adj_out[i].empty()) {
			const int j = *(--g.adj_out[i].end());
			g.remove_edge(i, j);
			if(!isRemoved[j]) {
				const bool self_loop = g.adj_out[j].count(j);
				if(self_loop || g.adj_in[j].size() <= 1 || g.adj_out[j].size() <= 1) {
					add2Stack(j, self_loop);
					goto dfs;
				}
			}
		}
		stack.pop_back();
	} while(!stack.empty());
	const int nn = g.n-removed.size();
	int r = 0;
	while(g.n > nn) {
		if(isRemoved[g.n-1]) g.remove_vertex(g.n-1);
		else {
			while(removed[r] >= nn) ++r;
			g.remove_vertex(removed[r++]);
		}
	}
	for(int i : removed) isRemoved[i] = false;
	removed.clear();
}

SCC::SCC(const Graph &g, bool use_pi): c(g.n, -1) {
	PROFIL_FUNC("SCC");
	assert(g.n);
	int next_id = 0;
	Vi id(g.n, -1), stack;
	vector<tuple<int,Si::const_iterator,int,int>> st;
	const auto add2Stack = [&](int i) {
		st.emplace_back(i, g.adj_out[i].begin(), id[i]=next_id++, stack.size());
		stack.push_back(i);
	};
	for(int i = 0; i < g.n; ++i) if(id[i] == -1) {
		add2Stack(i);
		dfs:
		auto &[i, it, min_id, start] = st.back();
		while(it != g.adj_out[i].end()) {
			const int j = *(it++);
			if(!use_pi && g.adj_in[i].count(j)) continue;
			if(id[j] == -2) continue;
			if(id[j] == -1) { add2Stack(j); goto dfs; }
			min_id = min(min_id, id[j]);
		}
		if(id[i] == min_id) {
			cs.emplace_back(stack.begin()+start, stack.end());
			stack.resize(start);
			for(int j : cs.back()) {
				c[j] = cs.size()-1;
				id[j] = -2;
			}
		}
		const int mid = min_id;
		st.pop_back();
		if(!st.empty()) {
			get<2>(st.back()) = min(get<2>(st.back()), mid);
			goto dfs;
		}
	}
	PROFIL_RET;
}

/*
 * Tests whether all cycles containing v also intersect at another vertex.
 * If yes, we can bypass v.
 */
bool has_two_petals(int v, Graph &g, Vi &parent, Vi &parent2, Vb &in_path, Vi& seen, int &count) {
	PROFIL_FUNC("Has petal");
	// First, we try to find a simple path an out_neighbor of v to v, using a BFS.
	queue<int> q;
	q.push(v); ++count;
	// Invariant: all vertices in Q have parent != -1 (except neighbors of v at init)
	while (!q.empty())
	{
		const int u = q.front(); q.pop();
		for (const int w: g.out_neighbors(u))
			if (seen[w] != count)
			{
				parent[w] = u;
				seen[w] = count;
				if(w == v) goto end_bfs;
				q.push(w);
			}
	}
	end_bfs:

	if (seen[v] != count) // There is no path from v to itself
		PROFIL_RET1(false);

	for (int u = parent[v]; u != v; u = parent[u])
		in_path[u] = true;

	// Try to find another x -> v path in the residual: 
	// we can take edges of the first path in the reverse direction only 
	q = queue<int>();
	q.push(v); ++ count;
	while (!q.empty())
	{
		const int u = q.front(); q.pop();
		for(int w: g.out_neighbors(u)) {
			if(in_path[w]) {
				w = parent[w];
				if(w == v) continue;
			}
			if(seen[w] == count) continue;
			parent2[w] = u;
			seen[w] = count;
			if(w == v) goto end_bfs2;
			q.push(w);
		}
		if(in_path[u]) {
			const int w = parent[u];
			if(w == v || seen[w] == count) continue;
			parent2[w] = u;
			seen[w] = count;
			q.push(w);
		}
	}
	end_bfs2:
	for (int u = parent[v]; u != v; u = parent[u])
		in_path[u] = false;

	PROFIL_RET1(seen[v] == count);
}

void remove_single_petal(Graph &g)
{
	PROFIL_FUNC("RM petal");
	int count = 0;
	Vi pred(g.n), pred2(g.n), seen(g.n, 0);
	Vb in_path(g.n, 0);
	bool does_something = true;
	do
	{
		does_something = false;
		for (int i = 0; i < g.n; ++i)
			if (!has_two_petals(i, g, pred, pred2, in_path, seen, count))
			{
				removeVertex(g, i, false);
				does_something = true;
			}
	}
	while (does_something);
	PROFIL_RET;
}

bool cliqueReduction(Graph &g, int u) {
	for(const Si *e : {&g.adj_in[u], &g.adj_out[u]}) {
		bool bad = false;
		for(auto it = e->begin(); it != e->end(); ++it)
			for(auto it2 = it; ++it2 != e->end();)
				if(!g.adj_in[*it].count(*it2) || !g.adj_out[*it].count(*it2))
					{ bad = true; goto endCheck; }
		endCheck:
		if(bad) continue;
		Kernel::removeVertex(g, u, false);
		return true;
	}
	return false;
}

void cliqueReduction(Graph &g) {
	PROFIL_FUNC("Clique Reduction");
	bool improve;
	do {
		improve = false;
		for(int u = g.n-1; u >= 0; --u) if(cliqueReduction(g, u)) {
			u = min(u, g.n);
			improve = true;
		}
	} while(improve && g.n);
	PROFIL_RET;
}

bool PIedge(Graph &g) {
	if(!g.n) return false;
	PROFIL_FUNC("PI edge reduc");
	SCC scc(g, false);
	if(scc.cs.size() == 1) PROFIL_RET1(false);
	Vi to_rm;
	bool simp = false;
	for(int u = 0; u < g.n; ++u) {
		for(int v : g.adj_out[u]) if(scc.c[u] != scc.c[v] && !g.adj_in[u].count(v)) to_rm.push_back(v);
		if(to_rm.empty()) continue;
		simp = true;
		for(int v : to_rm) g.remove_edge(u, v);
		to_rm.clear();
	}
	PROFIL_RET1(simp);
}

bool domination(Graph &g) {
	if(!g.n) return false;
	PROFIL_FUNC("Domination reduc");
	Vi to_rm;
	bool simp = false;
	for(int u = 0; u < g.n; ++u) {
		for(int v : g.adj_out[u]) if(!g.adj_in[u].count(v)) {
			for(int w : g.adj_in[u]) if(!g.adj_out[u].count(w) && !g.adj_in[v].count(w)) goto second;
			to_rm.push_back(v);
			continue;
			second:
			for(int w : g.adj_out[v]) if(!g.adj_in[v].count(w) && !g.adj_out[u].count(w)) goto end;
			to_rm.push_back(v);
			end:
			continue;
		}
		if(to_rm.empty()) continue;
		simp = true;
		for(int v : to_rm) g.remove_edge(u, v);
		to_rm.clear();
	}
	PROFIL_RET1(simp);
}

void simplify0(Graph &g) {
	for(int i = g.n-1; i >= 0; --i) if(g.adj_in[i].size() <= 1 || g.adj_out[i].size() <= 1) {
		removeVertex(g, i, false);
		i = min(i, g.n);
	}
	cliqueReduction(g);
}

void simplify1(Graph &g) {
	bool modif;
	do {
		simplify0(g);
		modif = false;
		modif |= PIedge(g);
		modif |= domination(g);
	} while(modif);
	if(g.n < 80) Kernel::remove_single_petal(g);
}

std::pair<std::vector<Graph>, Vi> simplify(Graph &g) {
	simplify1(g);
	if(!g.n) return {{}, g.solution};
	g.clearMoves();
	SCC scc(g);
	if(scc.cs.size() == 1) {
		if(g.n < 500) remove_single_petal(g);
		if(!g.n) return {{}, g.solution};
		Vi sol = move(g.solution);
		return {{move(g)}, move(sol)};
	}
	Vi newInd(g.n);
	std::vector<Graph> ans;
	Vi sol = move(g.solution);
	for(auto &c : scc.cs) {
		const int C = scc.c[c[0]];
		Graph g2(c.size(), g.inv_index);
		for(int i = 0; i < (int) c.size(); ++i) {
			newInd[c[i]] = i;
			g2.index[i] = g.index[c[i]];
			g.inv_index->at(g2.index[i]) = i;
		}
		for(int i = 0; i < (int) c.size(); ++i) {
			g2.adj_in [i] = move(g.adj_in [c[i]]);
			g2.adj_out[i] = move(g.adj_out[c[i]]);
			for(auto e : {&g2.adj_in [i], &g2.adj_out[i]}) {
				Vi tmp(e->begin(), e->end());
				int j = 0;
				while(j < (int) tmp.size())
					if(scc.c[tmp[j]] == C) {
						tmp[j] = newInd[tmp[j]];
						++ j;
					} else {
						tmp[j] = tmp.back();
						tmp.pop_back();
					}
				*e = Si(tmp.begin(), tmp.end());
			}
			g2.m += g2.adj_in[i].size();
		}
		auto [gs, s] = simplify(g2);
		ans.reserve(ans.size() + gs.size());
		move(gs.begin(), gs.end(), back_inserter(ans));
		sol.reserve(sol.size() + s.size());
		sol.insert(sol.end(), s.begin(), s.end());
	}
	return {move(ans), move(sol)};
}

}