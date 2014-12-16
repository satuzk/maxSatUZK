
#include <vector>

namespace maxsatuzk {

Weight countSoftWeight(InClauseSpace &f);

std::vector<int> convertBase(int64_t num, std::vector<int> &base);

template<typename Solver>
void solve(InClauseSpace &in, std::vector<int> &base, Weight initial_lb, Weight initial_ub);

} // namespace maxsatuzk

