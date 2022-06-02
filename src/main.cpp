#include <iostream>

#include "common.h"
#include "saumz1V9.h"

using namespace std;

int main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	Graph g = Graph::from_cin();
	PROFIL_INIT();
	auto r = pleaseSaumzDo1V9(g);
	PROFIL_SHOW();
	printSolution(r);
	return 0;
}
