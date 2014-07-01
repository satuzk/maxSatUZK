
#include <vector>

namespace maxsatuzk {

class Sorter {
public:
	Sorter(int size) {
		p_inputs.resize(size, Literal::null());
		p_outputs.resize(size, Literal::null());
	}

	Literal &input(int index) {
		return p_inputs[index];
	}
	Literal &output(int index) {
		return p_outputs[index];
	}

	int size() {
		return p_inputs.size();
	}

	std::vector<Literal> p_inputs;
	std::vector<Literal> p_outputs;
};

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

Weight countSoftWeight(ClauseSpace &f);

NumberSeq convertBase(int64_t num, NumberSeq &base);

void solve(ClauseSpace &in, NumberSeq &base, Weight initial_lb, Weight initial_ub);

} // namespace maxsatuzk

