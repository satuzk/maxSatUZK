
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

NumberSeq convertBase(int64_t num, NumberSeq &base) {
	NumberSeq val(base.length());
	val[0] = 1;
	for(int i = 1; i < base.length(); i++) {
		val[i] = val[i - 1] * base[i];
	}

	NumberSeq ns(base.length());
	for(int i = base.length() - 1; i >= 0; i--) {
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

