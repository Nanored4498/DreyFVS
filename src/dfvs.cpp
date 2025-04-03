#include "dfvs.h"

#include "upper_bound.h"
#include "common.h"
#include "kernel.h"

#include <algorithm>
#include <queue>
#include <chrono>
#include <csignal>

using namespace std;

using VVi = vector<Vi>;

static const int64_t SEED = 12012340032;
static auto rengine = std::default_random_engine(SEED);

static int SOL = 0;
static VVi best_sols;

static void showSol(int) {
	cerr << "SCORE: " << SOL << '\n';
	Vi sol;
	for(Vi &s : best_sols) move(s.begin(), s.end(), back_inserter(sol));
	printSolution(sol);
	PROFIL_SHOW();
	exit(0);
}

static int delta_swap(int u, int v, const VVi& vvn, const VVi& inv_vvn, Vi &pos, Vi &nb_back) {
	static Vi add;
	static Vb inQ;
	if(vvn.size() > add.size()) {
		add.resize(vvn.size(), 0);
		inQ.resize(vvn.size(), false);
	}
	// Swap pos of u and v (now pos[u] < pos[v])
	if(pos[u] < pos[v]) swap(u, v);
	swap(pos[u], pos[v]);

	const auto cmp = [&](int a, int b) { return pos[a] > pos[b]; };
	priority_queue<int, Vi, decltype(cmp)> q(cmp);
	const auto pq = [&](int i) { if(!inQ[i]) { q.push(i); inQ[i] = true; } };

	if(nb_back[u]) {
		for(const int w : vvn[u])
			if(!nb_back[w] && pos[u] < pos[w] && pos[w] <= pos[v])
				-- add[u];
		if(add[u] != 0) pq(u);
	} else {
		for(const int w : inv_vvn[u])
			if(pos[u] < pos[w] && pos[w] <= pos[v]) {
				++ add[w];
				pq(w);
			}
	}
	if(nb_back[v]) {
		for(const int w : vvn[v])
			if(!nb_back[w] && pos[u] < pos[w] && pos[w] < pos[v])
				++ add[v];
	} else {
		for(const int w : vvn[v])
			if(!nb_back[w] && pos[w] < pos[v] && w != u)
				++ add[v];
		for(const int w : inv_vvn[v])
			if(pos[u] < pos[w] && pos[w] < pos[v]) {
				-- add[w];
				pq(w);
			}
	}
	if(add[v] != 0) pq(v);

	int r = 0;
	while(!q.empty()) {
		const int x = q.top(); q.pop();
		const int last = nb_back[x];
		inQ[x] = false;
		nb_back[x] += add[x];
		add[x] = 0;
		if(last) {
			if(nb_back[x]) continue;
			--r;
			for(const int y : inv_vvn[x])
				if(pos[x] < pos[y]) {
					++ add[y];
					pq(y);
				}
		} else {
			if(!nb_back[x]) continue;
			++r;
			for(const int y : inv_vvn[x])
				if(pos[x] < pos[y]) {
					-- add[y];
					pq(y);
				}
		}
	}
	return r;
}

static Vi loop;
static bool find_cycle_bw_fw(const VVi& vin, const VVi& vout, const Vb& rmd, int start, bool genLoop=false) {
	static Vb reached;
	static Vi stack, prev;
	if(reached.size() < vin.size()) {
		reached.resize(vin.size(), 0);
		prev.resize(2*vin.size());
	}
	// TODO: merge q and stack
	queue<int> q;
	q.emplace(start<<1);
	q.emplace((start<<1)|1);
	reached[start] = 0b11;
	stack.push_back(start);
	while(!q.empty()) {
		const int cur = q.front();
		const int dir = cur&1;
		const Vi& neiv = dir == 0 ? vin[cur>>1] : vout[cur>>1];
		q.pop();
		for(const int n : neiv) if(!rmd[n]) {
			if(reached[n] & (1 << (dir ^ 1))) {
				for(int i : stack) reached[i] = 0;
				stack.clear();
				if(genLoop) {
					loop.clear();
					for(const int i0 : {cur,(n<<1)|(dir^1)})
						for(int i = i0; (i>>1)!=start; i = prev[i])
							loop.push_back(i>>1);
				}
				return true;
			}
			if(reached[n]) continue;
			reached[n] |= (1 << dir);
			stack.push_back(n);
			prev[(n<<1)|dir] = cur;
			q.push((n<<1)|dir);
		}
	}
	for(int i : stack) reached[i] = 0;
	stack.clear();
	return false;
}

std::vector<int> computeDFVS(Graph& g) {
	const int niters = 7000;
	const int niters2 = 50000;
	const double maxTime = 430.e9;


	// Init solution
	const chrono::high_resolution_clock::time_point st_time = chrono::high_resolution_clock::now();
	auto [gs, sol0] = Kernel::reduce_and_split(g);
	cerr << "gs size : " << gs.size() << '\n';
	int totM = 0;
	for(const Graph &g : gs) totM += g.m;
	best_sols.resize(gs.size()+1);
	best_sols.back() = std::move(sol0);
	for(int i = 0; i < (int)gs.size(); ++i) best_sols[i] = getUpperBound(gs[i], -1);
	for(const Vi &s : best_sols) SOL += s.size();
	const double ub_time = (chrono::high_resolution_clock::now() - st_time).count();

	// Setup signal handler
	struct sigaction sigact;
	sigact.sa_handler = showSol;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigaction(SIGINT, &sigact, NULL);
	sigaction(SIGTERM, &sigact, NULL);
	Vi tmp;

	for(int gi = 0; gi < (int)gs.size(); ++gi) {
		const Graph &g = gs[gi];
		Vi &sol = best_sols[gi];
		chrono::high_resolution_clock::time_point st_time = chrono::high_resolution_clock::now();
		VVi vin(g.n), vout(g.n);
		for(int i = 0; i < g.n; ++i) {
			vout[i].assign(g.adj_in[i].begin(), g.adj_in[i].end());
			vin[i].assign(g.adj_out[i].begin(), g.adj_out[i].end());
		}
		int restart = 0;
		Vi deg_v(g.n);
		Vi ids(g.n);
		iota(ids.begin(), ids.end(), 0);
		int nit = min(niters, g.m);
		int lastit;

		start:
		// starting point
		const Vi ub = restart ? getUpperBound(g, g.n > 2000 ? max(1, (int)sol.size() + int(rengine()%restart) - restart/2) : (int)sol.size()) : sol;
		// if better update solution
		if(ub.size() < sol.size()) {
			SOL -= sol.size() - ub.size();
			sol = ub;
		}
		cerr << "Ameliorations from " << SOL + ub.size() - sol.size() << endl;
		int count = ub.size();
		Vb rmd(g.n, false);
		for(const int i : ub) rmd[g.inv_index->at(i)] = true;
		const double maxT = gs.size() == 1 && g.n < 1000 ? 550.e9 : (maxTime-ub_time) * double(g.m) / totM;
		lastit = 0;
		for(int iter = 0; iter < nit; ++iter) {
			if((chrono::high_resolution_clock::now() - st_time).count() > maxT) break;

			for(int i = 0; i < g.n; ++i) {
				int a=0, b=0;
				for(int j : vin[i]) if(!rmd[j]) ++a;
				for(int j : vout[i]) if(!rmd[j]) ++b;
				if(a > b) swap(a, b);
				deg_v[i] = 3*a+2*b;
			}
			if(iter&1) shuffle(ids.begin(), ids.end(), rengine);
			else sort(ids.begin(), ids.end(), [&](int i, int j) { return deg_v[i] < deg_v[j]; });
			for(int i : ids) if(rmd[i]) {
				rmd[i] = false;
				rmd[i] = find_cycle_bw_fw(vin, vout, rmd, i, true);
				if(rmd[i]) {
					if(g.n < 15000 && (rengine()&1)) continue;
					if(rengine()&1) sort(loop.begin(), loop.end(), [&](int i, int j) { return deg_v[i] > deg_v[j]; });
					else shuffle(loop.begin(), loop.end(), rengine);
					for(const int j : loop) {
						swap(rmd[i], rmd[j]);
						if(find_cycle_bw_fw(vin, vout, rmd, i)) swap(rmd[i], rmd[j]);
						else break;
					}
				} else { // improve solution
					lastit = iter;
					if(--count < (int) sol.size()) {
						SOL -= sol.size() - count;
						cerr << "Iter " << iter << ": " << SOL << '\n';
						// TODO: Why using tmp???
						tmp.clear();
						for(int i = 0; i < g.n; ++i) if(rmd[i]) tmp.push_back(g.index[i]);
						swap(tmp, sol);
					}
				}
			}	
		}
		const auto ti = chrono::high_resolution_clock::now();
		const double elapsed = (ti - st_time).count();
		if(++restart < g.n && (maxT - elapsed) > .85 * elapsed / restart) {
			if(lastit < .6*nit) nit *= .95;
			else if(lastit > .9*nit) nit *= 1.05;
			goto start;
		}
	}

	double maxTime2 = 600.e9 - (chrono::high_resolution_clock::now() - st_time).count();
	for(int gi = 0; gi < (int)gs.size(); ++gi) {
		const Graph &g = gs[gi];
		Vi &sol = best_sols[gi];
		chrono::high_resolution_clock::time_point st_time = chrono::high_resolution_clock::now(), end_time;
		VVi vvn(g.n), inv_vvn(g.n);
		for(int i = 0; i < g.n; ++i) {
			inv_vvn[i].assign(g.adj_in[i].begin(), g.adj_in[i].end());
			vvn[i].assign(g.adj_out[i].begin(), g.adj_out[i].end());
		}
		Vi order, pos(g.n, 0), nb_backs(g.n, 0), res;
	
		order.reserve(g.n);
		Vb rmd(g.n, false);
		for(int i : sol) rmd[g.inv_index->at(i)] = true;
		for(int u = 0; u < g.n; ++u) if(!rmd[u]) {
			for(const int v : vvn[u]) if(!rmd[v]) ++ pos[u];
			if(!pos[u]) order.push_back(u);
		}
		for(int i = 0; i < (int) order.size(); ++i)
			for(const int v : inv_vvn[order[i]])
				if(!rmd[v] && (--pos[v])==0) order.push_back(v);
		reverse(order.begin(), order.end());
		const size_t s0 = order.size();
		for(int u = 0; u < g.n; ++u) if(rmd[u]) order.push_back(u);
		shuffle(order.begin()+s0, order.end(), rengine);
		for(int u = 0; u < g.n; ++u) pos[order[u]] = u;
		for(int u : order) for(int v : vvn[u])
			if(pos[v] < pos[u] && nb_backs[v] == 0)
				++nb_backs[u];
		uniform_int_distribution<int> all_vert_dist(0, g.n-1);
		cerr << "Ameliorations 2 from " << SOL << '\n';
		const double maxT = maxTime2 * double(g.m) / totM;
		for(int iter = 0; iter < niters2; ++iter) {
			if((chrono::high_resolution_clock::now() - st_time).count() > maxT) break;
			for(int u = 0; u < g.n; ++u) {
				const int v = all_vert_dist(rengine);
				if(u == v) continue;
				const int delta = delta_swap(u, v, vvn, inv_vvn, pos, nb_backs);
				if(delta < 0) {
					SOL += delta;
					cerr << "It " << iter << ": " << SOL << '\n';
					tmp.clear();
					for(int u = 0; u < g.n; ++u) if(nb_backs[u]) tmp.push_back(g.index[u]);
					swap(tmp, sol);
				} else if(delta > 0) delta_swap(u, v, vvn, inv_vvn, pos, nb_backs);
			}
		}
	}

	cerr << "SCORE: " << SOL << '\n';
	Vi sol; for(Vi &s : best_sols) move(s.begin(), s.end(), back_inserter(sol));
	return sol;
}

