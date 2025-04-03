#include "upper_bound.h"

#include "common.h"
#include "kernel.h"
#include <algorithm>
#include <cmath>

using namespace std;

// Construct a first solution using Sinkhorn–Knopp algorithm
static void getUpperBound0(Graph &g, int guess) {
	constexpr int64_t MAX_OP = 5e9;

	const double mean_val = sqrt(double(g.n)/double(g.n+g.m));
	vector<double> a(g.n), b(g.n, mean_val);
	vector<pair<double, int>> diag(g.n);
	while(g.n) {
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
			const int i = g.inv_index->at(ind);
			if(i >= g.n || g.index[i] != ind) continue; // vertex already removed
			Kernel::removeVertex(g, i, inSol);
		}

		// Apply CORE reduction and split in SCCs
		Kernel::core(g);
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
			Kernel::reduce(g2);
			const int64_t guess2 = max(1, guess - (int)g.solution.size());
			getUpperBound0(g2, (guess2 * g2.n)/n2);
			n2 -= g2.n;
			g.solution.insert(g.solution.end(), g2.solution.begin(), g2.solution.end());
		}
		return;
	}
}

Vi getUpperBound(const Graph &g, int guess) {
	PROFIL_FUNC("Upper bound");

	// If no guess is given, set guess to 0.9 * [solution obtained a greedy algorithm over degree]
	if(guess == -1) {
		Vi order(g.n);
		iota(order.begin(), order.end(), 0);
		sort(order.begin(), order.end(), [&](int u, int v){ return g.deg(u) < g.deg(v); });
		Graph g2(g, Graph::NoCopySol{});
		for(const int u : order) {
			const int i = g.index[u];
			const int u2 = g.inv_index->at(i);
			if(u2 >= g2.n || g2.index[u2] != i) continue; // vertex has already been removed
			Kernel::removeVertex(g2, u2, true);
		}
		for(int u = 0; u < g.n; ++u) g.inv_index->at(g.index[u]) = u;
		// TODO: g.solution.size() + 0.9 * g2.solution.size() would be better
		// but constant 0.9 may change
		guess = .9*(g.solution.size() + g2.solution.size());
	}

	// Use Sinkhorn–Knopp algorithm
	Graph g2(g, Graph::NoCopySol{});
	getUpperBound0(g2, guess - g.solution.size());
	for(int u = 0; u < g.n; ++u) g.inv_index->at(g.index[u]) = u;

	// Keep only a subset of the g2 solution
	vector<bool> del(g.n, false);
	Vi seen(g.n, 0), st;
	int S = 0;
	for(int i : g2.solution) del[g.inv_index->at(i)] = true;
	Vi sol;
	while(!g2.solution.empty()) {
		const int u = g.inv_index->at(g2.solution.back());
		g2.solution.pop_back();
		++S;
		st.push_back(u);
		while(!st.empty()) {
			int v = st.back(); st.pop_back();
			for(const int w : g.adj_out[v]) {
				if(w == u) goto cycle_found;
				if(del[w] || seen[w] == S) continue;
				seen[w] = S;
				st.push_back(w);
			}
		}
		del[u] = false;
		continue;
		cycle_found:
		sol.push_back(g.index[u]);
		st.clear();
	}
	PROFIL_RET1(sol);
}
