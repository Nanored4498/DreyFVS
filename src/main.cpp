#include <iostream>

#include "common.h"
#include "dfvs.h"

using namespace std;

int main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	Graph g = Graph::from_cin();
	PROFIL_INIT();
	auto r = computeDFVS(g);
	PROFIL_SHOW();
	printSolution(r);
	return 0;
}
