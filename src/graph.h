#pragma once

#include "set.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <numeric>
#include <vector>

typedef Set<int> Si;
typedef std::vector<int> Vi;
typedef std::vector<uint8_t> Vb;

struct Graph {
	int n, m;
	std::vector<Si> adj_in, adj_out;
	Vi index;
	Vi solution;
	std::shared_ptr<Vi> inv_index;

	Graph() = default;
	Graph(int n, std::shared_ptr<Vi> &inv): n(n), m(0), adj_in(n), adj_out(n), index(n), inv_index(inv) {
		iota(index.begin(), index.end(), 0);
	}
	struct NoCopySol {};
	Graph(const Graph &g, NoCopySol): n(g.n), m(g.m), adj_in(g.adj_in), adj_out(g.adj_out),
		index(g.index), inv_index(g.inv_index) {}

    void add_edge(int u, int v);
    void remove_edge(int u, int v);
	void remove_vertex(int u);

	inline int deg(int u) const { return adj_in[u].size()+adj_out[u].size(); }

	static Graph from_istream(std::istream &is);
	inline static Graph from_cin() { return from_istream(std::cin); }
	inline static Graph from_file(const std::string &fname) { std::ifstream ifs(fname); return from_istream(ifs); }
	void print(std::ostream &os) const;
};
