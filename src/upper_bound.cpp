#include "upper_bound.h"

#include "common.h"
#include "kernel.h"
#include <algorithm>
#include <cmath>
#include <cassert>

using namespace std;

void getUpperBound0(Graph &g, int guess) {
	constexpr int64_t MAX_OP = 5e9;

	// Greedily take the best vertex
	const double mean_val = sqrt(double(g.n)/double(g.n+g.m));
	vector<double> a(g.n, mean_val), b(g.n, mean_val);
	vector<pair<double, int>> diag(g.n);
	while(g.n) {
		assert(g.n > 1);
		int ln = ceil(log2(g.n));
		const int K = clamp((int64_t(2*g.m+g.n) * (g.n*ln)) / MAX_OP, (int64_t) 1, (int64_t) g.n-1);
		while(ln--) {
			for(int i = 0; i < g.n; ++i) {
				double s = b[i];
				for(int j : g.adj_out[i]) s += b[j];
				a[i] = 1/s;
			}
			for(int i = 0; i < g.n; ++i) {
				double s = a[i];
				for(int j : g.adj_in[i]) s += a[j];
				b[i] = 1/s;
			}
		}
		vector<pair<int, bool>> mvs(K);
		diag.resize(g.n);
		for(int i = 0; i < g.n; ++i) diag[i] = {a[i]*b[i], i};
		sort(diag.begin(), diag.end());
		const double med = diag[clamp(guess-(int)g.solution.size(), 1, g.n-2)].first;
		auto it0 = diag.begin(), it1 = --diag.end();
		for(int i = 0; i < K; ++i) 
			if(med/it0->first < it1->first/med) { mvs[i] = { g.index[it1->second], false }; --it1; }
			else { mvs[i] = { g.index[it0->second], true }; ++it0; }

		for(auto [ind, inSol] : mvs) {
			int i = g.inv_index->at(ind);
			if(i >= g.n || g.index[i] != ind) continue;
			// Remove the selected vertex i
			Kernel::removeVertex(g, i, inSol);
		}

		// Try more clique reduction and try splitting in different scc
		Kernel::cliqueReduction(g);
		if(!g.n) return;
		Kernel::SCC scc(g);
		if(scc.cs.size() < 2) continue;
		Vi newInd(g.n);
		int n2 = g.n;
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
			Kernel::simplify1(g2);
			int gu2 = max(1, guess - (int)g.solution.size());
			n2 -= g2.n;
			getUpperBound0(g2, ((int64_t)gu2 * g2.n)/(n2+g2.n));
			g.solution.insert(g.solution.end(), g2.solution.begin(), g2.solution.end());
		}
		return;
	}
}

Vi getUpperBound(const Graph &g, int guess) {
	PROFIL_FUNC("Upper bound");
	if(guess == -1) {
		Vi order(g.n);
		iota(order.begin(), order.end(), 0);
		sort(order.begin(), order.end(), [&](int u, int v){ return g.deg(u) < g.deg(v); });
		Graph g2(g, false);
		for(const int u : order) {
			const int i = g.index[u];
			const int u2 = g.inv_index->at(i);
			if(u2 >= g2.n || g2.index[u2] != i) continue;
			Kernel::removeVertex(g2, u2, true);
		}
		for(int u = 0; u < g.n; ++u) g.inv_index->at(g.index[u]) = u;
		guess = .9*(g.solution.size() + g2.solution.size());
	}
	Graph g2(g, false); getUpperBound0(g2, guess - g.solution.size());
	for(int u = 0; u < g.n; ++u) g.inv_index->at(g.index[u]) = u;
	vector<bool> del(g.n, false);
	Vi seen(g.n, 0), st; int S = 0;
	for(int i : g2.solution) del[g.inv_index->at(i)] = true;
	Vi sol;
	while(!g2.solution.empty()) {
		const int u = g.inv_index->at(g2.solution.back());
		g2.solution.pop_back();
		++S;
		st.push_back(u);
		bool found = false;
		while(!st.empty()) {
			int v = st.back(); st.pop_back();
			for(int w : g.adj_out[v])
				if(w == u) {
					found = true;
					goto end;
				} else if(!del[w] && seen[w] != S) {
					seen[w] = S;
					st.push_back(w);
				}
		}
		end:
		if(!found) del[u] = false;
		else {
			sol.push_back(g.index[u]);
			st.clear();
		}
	}
	PROFIL_RET1(sol);
}