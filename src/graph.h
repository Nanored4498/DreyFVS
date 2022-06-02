#pragma once

#include "set.h"

#include <iostream>
#include <fstream>
#include <numeric>
#include <memory>

typedef Set<int> Si;
typedef std::vector<int> Vi;
typedef std::vector<uint8_t> Vb;

struct Graph
{
	int n, m;
	std::vector<Si> adj_in, adj_out;
	Vi index;
	Vi solution;
	std::shared_ptr<Vi> inv_index;

	Graph() = default;
	Graph(int n, std::shared_ptr<Vi> inv): n(n), m(0), adj_in(n), adj_out(n), index(n), inv_index(inv) {
		iota(index.begin(), index.end(), 0);
	}
	Graph(const Graph &g, bool copy_sol=true): n(g.n), m(g.m), adj_in(g.adj_in), adj_out(g.adj_out),
							index(g.index), solution(copy_sol?g.solution:Vi()), inv_index(g.inv_index) {}
	Graph(Graph &&) = default;

	inline Graph& operator=(const Graph &g) {
		n = g.n; m = g.m;
		adj_in = g.adj_in; adj_out = g.adj_out;
		index = g.index;
		solution = g.solution;
		inv_index = g.inv_index;
		return *this;
	}

	inline friend void swap(Graph &a, Graph &b) {
		std::swap(a.n, b.n); std::swap(a.m, b.m);
		std::swap(a.adj_in, b.adj_in), std::swap(a.adj_out, b.adj_out);
		std::swap(a.index, b.index);
		std::swap(a.solution, b.solution);
		std::swap(a.inv_index, b.inv_index);
		std::swap(a.moves, b.moves);
	}

	void clear() {
		n = m = 0;
		adj_in.clear(); adj_in.shrink_to_fit();
		adj_out.clear(); adj_out.shrink_to_fit();
		index.clear(); index.shrink_to_fit();
		solution.clear(); solution.shrink_to_fit();
		inv_index = nullptr;
		clearMoves();
	}

    void add_edge(int u, int v);
    void remove_edge(int u, int v);
	void remove_vertex(int u);

	inline int deg(int u) const { return adj_in[u].size()+adj_out[u].size(); }
	inline const Si &in_neighbors (int u) const { return adj_in[u]; }
	inline const Si &out_neighbors(int u) const { return adj_out[u]; }

	std::pair<int, int> getTime() const { return {n_moves(), solution.size()}; };
	void undo(const std::pair<int, int> &time);
	inline int n_moves() const { return moves.size(); }
	inline void clearMoves() { moves.clear(); moves.shrink_to_fit(); }

	static Graph from_istream(std::istream &is);
	inline static Graph from_cin() { return from_istream(std::cin); }
	inline static Graph from_file(const std::string &fname) { std::ifstream ifs(fname); return from_istream(ifs); }
	void print(std::ostream &os) const;

private:
	struct Move {
		int u;
		union { int v, ind; };
		enum Op { ADD_EDGE, REMOVE_EDGE, REMOVE_VERTEX };
		Op op;
		Move(Op op, int a, int b): u(a), v(b), op(op) {}
	};
	std::vector<Move> moves;
};