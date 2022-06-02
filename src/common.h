#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <unordered_set>

using RNG = std::default_random_engine;

/** Helper functions to debug vectors / pairs **/
template<class U, class V>
std::ostream& operator<<(std::ostream &os, const std::pair<U, V> &p)
{
	os << "(" << p.first << ", " << p.second << ")";
	return os;
}

template<class T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v)
{
	os << "[";
	for (const T &x: v)
		os << x << ", "; 

	os << "]";
	return os;
}


template<class T>
std::ostream& operator<<(std::ostream &os, const std::unordered_set<T> &v)
{
	os << "{";
	for (const T &x: v)
		os << x << ", "; 

	os << "}";
	return os;
}

inline void printSolution(const std::vector<int> &sol) {
	for(int i : sol) std::cout << i+1 << '\n';
}

#ifdef PROFILING
#include <chrono>
#include <tuple>

namespace Profiling {
inline std::chrono::high_resolution_clock::time_point t0;
inline std::vector<std::tuple<std::string, int, std::chrono::high_resolution_clock::duration>> stats;
inline void init() {
	t0 = std::chrono::high_resolution_clock::now();
}
inline void printStat() {
	const std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	const int64_t dt = (t1 - t0).count();
	std::cerr << "=== PROFILING ===" << std::endl;
	for(const auto &[name, call, time] : stats)
		std::cerr << name << ":\r\t\t\t" << call << "\r\t\t\t\t\t" << double(100*time.count())/dt << "%" << std::endl;
	std::cerr << "=================" << std::endl;
}
}

#define PROFIL_FUNC(name)\
static int __prof_index__ = -1;\
if(__prof_index__ == -1) {\
	__prof_index__ = Profiling::stats.size();\
	Profiling::stats.emplace_back((name), 0, 0);\
}\
++std::get<1>(Profiling::stats[__prof_index__]);\
std::chrono::high_resolution_clock::time_point __time0__ = std::chrono::high_resolution_clock::now()

#define PROFIL_RET std::get<2>(Profiling::stats[__prof_index__]) += std::chrono::high_resolution_clock::now() - __time0__
#define PROFIL_RET1(ret){\
auto x = (ret);\
PROFIL_RET;\
return x;\
}

#define PROFIL_INIT() Profiling::init()
#define PROFIL_SHOW() Profiling::printStat()

#else
#define PROFIL_FUNC(name)
#define PROFIL_RET
#define PROFIL_RET1(ret) return (ret)
#define PROFIL_INIT()
#define PROFIL_SHOW()
#endif
