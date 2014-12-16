
#include <vector>

namespace maxsatuzk {

class NumberSeq {
public:
	NumberSeq(int length) {
		p_seq.resize(length, 0);
	}

	int &operator[] (int index) {
		return p_seq[index];
	}

	int length() {
		return p_seq.size();
	}
private:
	std::vector<int> p_seq;
};

Weight countSoftWeight(InClauseSpace &f);

NumberSeq convertBase(int64_t num, NumberSeq &base);

template<typename Solver>
void solve(InClauseSpace &in, NumberSeq &base, Weight initial_lb, Weight initial_ub);

} // namespace maxsatuzk

