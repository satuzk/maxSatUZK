
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "../include/Formula.hpp"
#include "../include/Sorter.hpp"
#include "../inline/Formula.hpp"

namespace maxsatuzk {

class Base {
public:
	int length() const {
		return p_vector.size();
	}
	
	void extend(int p) {
		p_vector.push_back(p);
	}

	int product() const {
		return std::accumulate(p_vector.begin(), p_vector.end(),
				1, std::multiplies<int>());
	}

	int operator[] (int index) const {
		return p_vector[index];
	}

	NumberSeq convert(int num) const {
		NumberSeq val(p_vector.size() + 1);
		val[0] = 1;
		for(int i = 0; i < p_vector.size(); i++) {
			val[i + 1] = val[i] * p_vector[i];
		}

		NumberSeq seq(p_vector.size() + 1);
		for(int i = p_vector.size(); i >= 0; i--) {
			int k = num / val[i];
			seq[i] = k;
			num -= k * val[i];
		}
		return seq;
	}
	
private:
	std::vector<int> p_vector;
};

int cost(const Base &base, const std::vector<int> &s) {
	int sum = 0;
	for(int i = 0; i < s.size(); i++) {
		NumberSeq seq = base.convert(s[i]);
		for(int j = 0; j < base.length() + 1; j++)
			sum += seq[j];
	}
	return sum;
}

// for each base b' extending a base b we have cost(b') >= partial(b)
// i.e. if the partial cost of a base exceeds the cost of the best known
// base we can prun it from the search tree
int partial(const Base &base, const std::vector<int> &s) {
	int sum = 0;
	for(int i = 0; i < s.size(); i++) {
		NumberSeq seq = base.convert(s[i]);
		for(int j = 0; j < base.length(); j++)
			sum += seq[j];
	}
	return sum;
}

// for each base b' extending a base b we have
// cost(b') >= partial(b') + heuristic(b') >= partial(b) + heuristic(b)
// i.e. if the partial cost + heuristic of a base exceed the cost of the
// best known base we can prun it from the search tree
int heuristic(const Base &base, const std::vector<int> &s) {
	int count = 0;
	for(int i = 0; i < s.size(); i++) {
		if(s[i] > base.product())
			count++;
	}
	return count;
}

void dfsBase(const Base &current, const std::vector<int> &s, int max, Base &best) {
	if(partial(current, s) + heuristic(current, s) > cost(best, s))
		return;
	
	if(cost(current, s) < cost(best, s))
		best = current;

	for(int p = 2; true; p++) {
		if(current.product() * p > max)
			break;

		Base next = current;
		next.extend(p);
		dfsBase(next, s, max, best);
	}
}

NumberSeq optimalBase(const std::vector<int> &s) {
	Base best;
	if(s.size() > 0) {
		int max = *std::max_element(s.begin(), s.end());
		while(!(best.product() * 2 > max))
			best.extend(2);
		
		dfsBase(Base(), s, max, best);
	}
	
	NumberSeq seq(best.length() + 1);
	seq[0] = 1;
	std::cout << "c Base:";
	for(int i = 0; i < best.length(); i++) {
		std::cout << " " << best[i];
		seq[i + 1] = best[i];
	}
	std::cout << std::endl;
	return seq;	
}

} // namespace maxsatuzk

