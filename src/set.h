#pragma once

#include <vector>
#include <cassert>

template<typename T>
struct Set {
protected:
	typedef std::pair<T, int> Node;
	typedef typename std::vector<Node>::const_iterator vec_iterator;

	struct Iterator : vec_iterator {
		typedef const T& reference;
		typedef const T* pointer;

		constexpr Iterator(const vec_iterator &it): vec_iterator(it) {}

		constexpr reference operator*() const noexcept { return vec_iterator::operator*().first; }
		constexpr pointer operator->() const noexcept { return &operator*(); }

		constexpr Iterator operator++(int) noexcept { return vec_iterator::operator++(0); }
		constexpr Iterator& operator++() noexcept { vec_iterator::operator++(); return *this; }
		constexpr Iterator operator--(int) noexcept { return vec_iterator::operator--(0); }
		constexpr Iterator& operator--() noexcept { vec_iterator::operator--(); return *this; }
	};

	std::vector<Node> v;
	std::vector<int> buckets;

public:
	typedef Iterator iterator;
	typedef Iterator const_iterator;

	Set() = default;

	template<typename InIterator>
	Set(InIterator first, InIterator last): v(), buckets(nextSize(std::distance(first, last)), -1) {
		while(first != last) insert(*(first++));
	}

	std::size_t size() const noexcept { return v.size(); }
	bool empty() const noexcept { return v.empty(); }

	iterator begin() const noexcept { return v.begin(); }
	iterator end()   const noexcept { return v.end(); }

	inline int key(const T &x) const {
		return std::hash<T>()(x) % buckets.size();
	}

	bool insert(const T &x) {
		if(buckets.empty()) buckets.assign(primes[0], -1);
		int b = key(x);
		for(int i = buckets[b]; i != -1; i = v[i].second)
			if(v[i].first == x) return false;
		if(v.size()+1 > buckets.size()) { //rehash
			rehash(nextSize(v.size()+1));
			b = key(x);
		}
		v.emplace_back(x, buckets[b]);
		buckets[b] = v.size()-1;
		return true;
	}

	void erase(const T &x) {
		if(buckets.empty()) return;
		const int b = key(x);
		if((v.size()<<2) <= buckets.size()) {
			for(int i = buckets[b]; i != -1; i = v[i].second)
				if(v[i].first == x) {
					if(i+1 != (int)v.size()) v[i] = move(v.back());
					v.pop_back();
					return rehash(nextSize(v.size()));
				}
			return;
		}
		int i = buckets[b];
		if(i == -1) return;
		if(v[i].first == x) {
			buckets[b] = v[i].second;
			return pop(i);
		}
		while(v[i].second != -1 && v[v[i].second].first != x) i = v[i].second;
		const int j = v[i].second;
		if(j == -1) return;
		v[i].second = v[j].second;
		pop(j);
	}
	inline void erase(iterator it) { erase(*it); }

	bool count(const T &x) const {
		if(buckets.empty()) return false;
		const int b = key(x);
		for(int i = buckets[b]; i != -1; i = v[i].second)
			if(v[i].first == x) return true;
		return false;
	}

protected:
	static constexpr size_t primes[] = {7, 17, 37, 79, 163, 331, 673, 1361, 2729, 5471, 10949, 21911,
										43853, 87719, 175447, 350899, 701819, 1403641, 2807303, 5614657,
										11229331, 22458671, 44917381, 89834777, 179669557, 359339171};
	
	static size_t nextSize(size_t n) {
		auto it = std::lower_bound(primes, primes+std::size(primes), n);
		assert(it != primes+std::size(primes));
		return *it;
	}

	void pop(int i) {
		if(i+1 != (int)v.size()) {
			const int b = key(v.back().first);
			int j = buckets[b];
			if(j+1 == (int)v.size()) {
				buckets[b] = i;
			} else {
				while(v[j].second+1 != (int)v.size()) j = v[j].second;
				v[j].second = i;
			}
			v[i] = move(v.back());
		}
		v.pop_back();
	}

	void rehash(size_t s) {
		buckets.assign(s, -1);
		for(int i = 0; i < (int) v.size(); ++i) {
			const int b = key(v[i].first);
			v[i].second = buckets[b];
			buckets[b] = i;
		}
	}
};
