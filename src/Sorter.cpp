
#include <time.h>

#include <cstdint>
#include <cassert>
#include <memory>
#include <vector>
#include <deque>
#include <unordered_map>
#include <iostream>
#include <thread>
#include <algorithm>
#include <stdexcept>

#include "../include/Formula.hpp"
#include "../include/SolverUzk.hpp"
#include "../include/Sorter.hpp"
#include "../inline/Formula.hpp"

namespace maxsatuzk {

std::vector<int> convertBase(int64_t num, std::vector<int> &base) {
	std::vector<int> val(base.size());
	val[0] = 1;
	for(int i = 1; i < base.size(); i++) {
		val[i] = val[i - 1] * base[i];
	}

	std::vector<int> ns(base.size());
	for(int i = base.size() - 1; i >= 0; i--) {
		int64_t k = num / val[i];
		assert(k < std::numeric_limits<int>::max());
		ns[i] = k;
		num -= k * val[i];
	}
	return ns;
}

Weight countSoftWeight(InClauseSpace &f) {
	Weight count = 0;
	for(auto i = f.refsBegin(); i != f.refsEnd(); ++i)
		if((*i).getWeight() != kHardWeight)
			count += (*i).getWeight();
	return count;
}

} // namespace maxsatuzk

