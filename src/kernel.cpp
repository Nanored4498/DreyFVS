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
	const auto condAdd2Stack = [&](const int j)->void {
		if(isRemoved[j]) return;
		const bool self_loop = g.adj_out[j].count(j);
		if(self_loop || g.adj_in[j].size() <= 1 || g.adj_out[j].size() <= 1)
			add2Stack(j, self_loop);
	};
	add2Stack(i, inSol);
	do {
		const int i = stack.back();
		stack.pop_back();
		while(!g.adj_in[i].empty()) {
			const int j = *(--g.adj_in[i].end());
			g.remove_edge(j, i);
			condAdd2Stack(j);
		}
		while(!g.adj_out[i].empty()) {
			const int j = *(--g.adj_out[i].end());
			g.remove_edge(i, j);
			condAdd2Stack(j);
		}
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
static bool has_two_petals(int v, Graph &g, Vi &parent, Vi &parent2, Vb &in_path, Vi& seen, int &count) {
	PROFIL_FUNC("Has petal");
	// First, we try to find a simple path an out_neighbor of v to v, using a BFS.
	queue<int> q;
	q.push(v); ++count;
	// Invariant: all vertices in Q have parent != -1 (except neighbors of v at init)
	while (!q.empty())
	{
		const int u = q.front(); q.pop();
		for (const int w: g.adj_out[u])
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
		for(int w: g.adj_out[u]) {
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

// remove vertices having a single petal
static void remove_single_petal(Graph &g) {
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

// Apply CORE on vertex u
static bool core0(Graph &g, int u) {
	for(const Si *e : {&g.adj_in[u], &g.adj_out[u]}) {
		for(auto it = e->begin(); it != e->end(); ++it)
			for(auto it2 = it; ++it2 != e->end();)
				if(!g.adj_in[*it].count(*it2) || !g.adj_out[*it].count(*it2))
					goto notCore;
		Kernel::removeVertex(g, u, false);
		return true;
		notCore:
		continue;
	}
	return false;
}

void core(Graph &g) {
	PROFIL_FUNC("Clique Reduction");
	bool improve;
	do {
		improve = false;
		for(int u = g.n-1; u >= 0; --u) if(core0(g, u)) {
			u = min(u, g.n);
			improve = true;
		}
	} while(improve && g.n);
	PROFIL_RET;
}

// Pi-edge reduction
// Remove the acyclic edges of G-PIE
// Where PIE is the graph containing the bi-direction edges
static bool PIedge(Graph &g) {
	if(!g.n) return false;
	PROFIL_FUNC("PI edge reduc");
	SCC scc(g, false);
	if(scc.cs.size() == 1) PROFIL_RET1(false);
	Vi to_rm;
	bool simp = false;
	for(int u = 0; u < g.n; ++u) {
		for(int v : g.adj_out[u]) {
			if(scc.c[u] == scc.c[v]) continue;
			if(g.adj_in[u].count(v)) continue; // PI edge
			to_rm.push_back(v);
		}
		if(to_rm.empty()) continue;
		simp = true;
		for(int v : to_rm) g.remove_edge(u, v);
		to_rm.clear();
	}
	PROFIL_RET1(simp);
}

// Domination reduction
// Remove dominated edges
static bool domination(Graph &g) {
	if(!g.n) return false;
	PROFIL_FUNC("Domination reduc");
	Vi to_rm;
	bool simp = false;
	for(int u = 0; u < g.n; ++u) {
		for(int v : g.adj_out[u]) if(!g.adj_in[u].count(v)) { // Not Pi-edge (u, v)
			// If non Pi predecessors of u (w) are predecessors of v then (u, v) can be removed
			// w -> u -> v ==> w -> v
			for(int w : g.adj_in[u]) if(!g.adj_out[u].count(w) && !g.adj_in[v].count(w)) goto second;
			to_rm.push_back(v);
			continue;
			second:
			// If non Pi successors of v (w) are successors of u then (u, v) can be removed
			// u -> v -> w ==> u -> w
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

// Apply deg <= 1 reductions (IN0, OUT0, IN1, OUT1) and CORE' reduction
static void reduce0(Graph &g) {
	for(int i = g.n-1; i >= 0; --i) if(g.adj_in[i].size() <= 1 || g.adj_out[i].size() <= 1) {
		removeVertex(g, i, false);
		i = min(i, g.n);
	}
	core(g);
}

void reduce(Graph &g) {
	bool modif;
	do {
		reduce0(g);
		modif = false;
		modif |= PIedge(g);
		modif |= domination(g);
	} while(modif);
	if(g.n < 80) Kernel::remove_single_petal(g);
}

pair<vector<Graph>, Vi> reduce_and_split(Graph &g) {
	reduce(g);
	if(!g.n) return {{}, g.solution};
	SCC scc(g);
	if(scc.cs.size() == 1) {
		if(g.n < 500) remove_single_petal(g);
		if(!g.n) return {{}, g.solution};
		pair<vector<Graph>, Vi> ans;
		ans.second = std::move(g.solution);
		ans.first.emplace_back(std::move(g));
		return ans;
	}
	Vi newInd(g.n);
	auto ans = make_pair<vector<Graph>, Vi>(vector<Graph>{}, std::move(g.solution));
	auto &[gs, sol] = ans;
	for(auto &c : scc.cs) {
		const int C = scc.c[c[0]];
		Graph g2(c.size(), g.inv_index);
		for(int i = 0; i < (int) c.size(); ++i) {
			newInd[c[i]] = i;
			g2.index[i] = g.index[c[i]];
			g.inv_index->at(g2.index[i]) = i;
		}
		for(int i = 0; i < (int) c.size(); ++i) {
			g2.adj_in [i] = std::move(g.adj_in [c[i]]);
			g2.adj_out[i] = std::move(g.adj_out[c[i]]);
			for(Si* e : {&g2.adj_in[i], &g2.adj_out[i]}) {
				Vi tmp(e->begin(), e->end());
				for(int j = 0; j < (int) tmp.size();)
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
		auto [gs2, sol2] = reduce_and_split(g2);
		move(gs2.begin(), gs2.end(), back_inserter(gs));
		sol.insert(sol.end(), sol2.begin(), sol2.end());
	}
	return ans;
}

}