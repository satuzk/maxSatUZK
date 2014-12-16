
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

#include <encodeuzk/encode.hpp>
#include <encodeuzk/static.hpp>
#include <encodeuzk/encode.inline.hpp>
#include <encodeuzk/basic.inline.hpp>
#include <encodeuzk/sorting.inline.hpp>
#include <encodeuzk/static.inline.hpp>

namespace maxsatuzk {

static const bool debugSorters = false;
static const bool debugRhs = false;

template<typename Literal>
using SorterLits = std::vector<Literal>;

template<typename Literal>
using SorterNetwork = std::vector<SorterLits<Literal>>;

template<typename VarAllocator, typename ClauseEmitter>
SorterNetwork<typename ClauseEmitter::Literal>
computeSorterNetwork(VarAllocator &allocator, ClauseEmitter &emitter,
		const std::vector<typename ClauseEmitter::Literal> &lits,
		const std::vector<Weight> &weights, const std::vector<int> &base) {
	// we need a literal that is always zero to simplify sorting
	typename ClauseEmitter::Literal null_lit = allocator.allocate().oneLiteral();
	encodeuzk::emit(emitter, { null_lit.inverse() });
	
	SorterNetwork<typename ClauseEmitter::Literal> sorters;

	for(int k = 0; k < base.size(); k++) {
		std::vector<typename ClauseEmitter::Literal> ins;

		// add carry bits from previous sorter as input
		if(k > 0) {
			for(int j = base[k] - 1; j < sorters.back().size(); j += base[k]) {
				ins.push_back(sorters.back()[j]);
				if(debugSorters)
					std::cout << "carry from " << (k - 1) << " to " << k << std::endl;
			}
		}

		for(int i = 0; i < lits.size(); i++) {
			std::vector<int> weight = convertBase(weights[i], base);
			for(int j = 0; j < weight[k]; j++) {
				ins.push_back(lits[i]);
				if(debugSorters)
					std::cout << "s[" << k << "] weight" << std::endl;
			}
		}

		std::vector<typename ClauseEmitter::Literal> outs
				= encodeuzk::computePwSort(allocator, emitter, ins, null_lit);
		std::cout << "c Size of sorter " << k << ": " << outs.size() << std::endl;
		sorters.push_back(outs);
	}

	return sorters;
}

// produces a constraint that is true iff sorter >= target
template<typename VarAllocator, typename ClauseEmitter>
typename ClauseEmitter::Literal computeSorterGe(VarAllocator &allocator, ClauseEmitter &emitter,
		const SorterLits<typename ClauseEmitter::Literal> &sorter, int target) {
	if(target == 0) {
		typename ClauseEmitter::Variable r = allocator.allocate();

		// trivial case 1: every number is >= 0
		encodeuzk::emit(emitter, { r.oneLiteral() });

		return r.oneLiteral();
	}else if(sorter.size() < target) {
		typename ClauseEmitter::Variable r = allocator.allocate();

		// trivial case 2: the sorter is not big enough to reach the limit
		encodeuzk::emit(emitter, { r.zeroLiteral() });

		return r.oneLiteral();
	}

	return sorter[target - 1];
}

// produces a constraint that is true iff sorter % divisor >= target
template<typename VarAllocator, typename ClauseEmitter>
typename ClauseEmitter::Literal computeSorterRemainderGe(VarAllocator &allocator, ClauseEmitter &emitter,
		const SorterLits<typename ClauseEmitter::Literal> &sorter, int divisor, int target) {
	if(target == 0) {
		typename ClauseEmitter::Variable r = allocator.allocate();
	
		// trivial case 1: every number is >= 0
		encodeuzk::emit(emitter, { r.oneLiteral() });

		return r.oneLiteral();
	}else if(sorter.size() < target) {
		typename ClauseEmitter::Variable r = allocator.allocate();
	
		// trivial case 2: the sorter is not big enough to reach the limit
		encodeuzk::emit(emitter, { r.zeroLiteral() });
		
		return r.oneLiteral();
	}else if(divisor <= target) {
		typename ClauseEmitter::Variable r = allocator.allocate();
	
		// trivial case 3: the modulus is not big enough to reach the limit
		encodeuzk::emit(emitter, { r.zeroLiteral() });
		
		return r.oneLiteral();
	}
	
	std::vector<typename ClauseEmitter::Literal> disjunction;
	for(int k = 0; k < sorter.size(); k += divisor) {
		if(k + target - 1 >= sorter.size())
			break;

		if(k + divisor - 1 < sorter.size()) {
			disjunction.push_back(encodeuzk::computeAnd(allocator, emitter,
					sorter[k + target - 1], sorter[k + divisor - 1].inverse()));
		}else{
			disjunction.push_back(sorter[k + target - 1]);
		}
	}
	
	return encodeuzk::computeOrN(allocator, emitter,
			disjunction.begin(), disjunction.end());
}

// generates the constraint (sorters >= rhs)
template<typename VarAllocator, typename ClauseEmitter>
typename ClauseEmitter::Literal computeSorterNetworkGe(VarAllocator &allocator, ClauseEmitter &emitter,
		const SorterNetwork<typename ClauseEmitter::Literal> &sorters,
		const std::vector<int> &base, const std::vector<int> &rhs, int i) {

	if(i == 0) {
		typename ClauseEmitter::Variable r = allocator.allocate();

		// trivial case: every number is >= 0
		encodeuzk::emit(emitter, { r.oneLiteral() });
		if(debugRhs)
			std::cout << r.oneLiteral().toNumber() << " <-> true" << std::endl;

		return r.oneLiteral();
	}

	i--;
	
	typename ClauseEmitter::Literal  p
			= computeSorterNetworkGe(allocator, emitter, sorters, base, rhs, i);
	
	typename ClauseEmitter::Literal gt;
	typename ClauseEmitter::Literal ge;
	if(i == base.size() - 1) {
		gt = computeSorterGe(allocator, emitter, sorters[i], rhs[i] + 1);
		ge = computeSorterGe(allocator, emitter, sorters[i], rhs[i]);
	}else{
		gt = computeSorterRemainderGe(allocator, emitter, sorters[i], base[i + 1], rhs[i] + 1);
		ge = computeSorterRemainderGe(allocator, emitter, sorters[i], base[i + 1], rhs[i]);
	}
	
	return encodeuzk::computeOr(allocator, emitter, gt,
			encodeuzk::computeAnd(allocator, emitter, ge, p));
}

template<typename VarAllocator, typename ClauseEmitter>
typename ClauseEmitter::Literal computeSorterNetworkGe(VarAllocator &allocator, ClauseEmitter &emitter,
		const SorterNetwork<typename ClauseEmitter::Literal> &sorters,
		const std::vector<int> &base, const std::vector<int> &rhs) {
	return computeSorterNetworkGe(allocator, emitter, sorters, base, rhs,
			sorters.size());
}

template<typename VarAllocator, typename ClauseEmitter>
void maxsatHard(VarAllocator &allocator, ClauseEmitter &emitter,
		InClauseSpace &in,
		std::unordered_map<InVariable, typename ClauseEmitter::Variable, VarHashFunc> &variable_map) {
	// copy all variables from the original formula
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		InClauseRef in_clause = *i;
		for(int j = 0; j < in_clause.length(); j++) {
			InVariable in_variable = in_clause.getLiteral(j).var();
			if(variable_map.find(in_variable) == variable_map.end())
				variable_map[in_variable] = allocator.allocate();
		}
	}
	
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		InClauseRef in_clause = *i;
		if(in_clause.getWeight() != kHardWeight)
			continue;

		// copy the literals from the original clause
		std::vector<typename ClauseEmitter::Literal> out_clause;
		for(int j = 0; j < in_clause.length(); j++) {
			InLiteral in_literal = in_clause.getLiteral(j);
			typename ClauseEmitter::Variable out_variable = variable_map[in_literal.var()];
			out_clause.push_back(in_literal.sign() > 0 ? out_variable.oneLiteral() : out_variable.zeroLiteral());
		}
		emitter.emit(out_clause.begin(), out_clause.end());
	}
}

template<typename VarAllocator, typename ClauseEmitter>
void maxsatLhs(VarAllocator &allocator, ClauseEmitter &emitter,
		InClauseSpace &in, std::vector<int> &base,
		SorterNetwork<typename ClauseEmitter::Literal> &sorters,
		std::vector<typename ClauseEmitter::Literal> &rel_lits,
		std::vector<Weight> &rel_weights,
		std::unordered_map<InVariable, typename ClauseEmitter::Variable, VarHashFunc> &variable_map) {
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		InClauseRef in_clause = *i;
		if(in_clause.getWeight() == kHardWeight)
			continue;

		// copy the literals from the original clause
		std::vector<typename ClauseEmitter::Literal> out_clause;
		for(int j = 0; j < in_clause.length(); j++) {
			InLiteral in_literal = in_clause.getLiteral(j);
			typename ClauseEmitter::Variable out_variable = variable_map[in_literal.var()];
			out_clause.push_back(in_literal.sign() > 0 ? out_variable.oneLiteral() : out_variable.zeroLiteral());
		}

		// add a relaxation literal
		typename ClauseEmitter::Literal rel_lit = allocator.allocate().oneLiteral();
		out_clause.push_back(rel_lit);
		emitter.emit(out_clause.begin(), out_clause.end());

		rel_lits.push_back(rel_lit.inverse());
		rel_weights.push_back(in_clause.getWeight());
	}

	sorters = computeSorterNetwork(allocator, emitter, rel_lits, rel_weights, base);
}

// generates a formula, that is satisfiable iff at least min soft clauses can be satisfied
template<typename VarAllocator, typename ClauseEmitter>
void maxsatRhs(VarAllocator &allocator, ClauseEmitter &emitter,
		std::vector<int> &base, int min,
		SorterNetwork<typename ClauseEmitter::Literal> &sorters) {
	std::vector<int> rhs = convertBase(min, base);
	if(debugRhs)
		for(int i = 0; i < rhs.size(); i++)
			std::cout << "rhs[" << i << "]: " << rhs[i] << std::endl;

	typename ClauseEmitter::Literal r
			= computeSorterNetworkGe(allocator, emitter, sorters, base, rhs);
	if(debugRhs)
		std::cout << "objective: " << r.toNumber() << std::endl;

	encodeuzk::forceTrue(allocator, emitter, r);;
}

template<typename Solver>
void solve(InClauseSpace &in, std::vector<int> &base, Weight initial_lb, Weight initial_ub) {
	Weight cur_lb = initial_lb; // formula is sat for v <= lb
	Weight cur_ub = initial_ub; // formula is unsat for v >= ub

	Solver solver(0);
	typename Solver::HardAllocator hard_allocator(solver);
	typename Solver::HardEmitter hard_emitter(solver);

	std::unordered_map<InVariable, typename Solver::Variable,
			VarHashFunc> variable_map;
	maxsatHard(hard_allocator, hard_emitter, in, variable_map);

	SorterNetwork<typename Solver::Literal> sorters;
	std::vector<typename Solver::Literal> rel_lits;
	std::vector<Weight> rel_weights;
	maxsatLhs(hard_allocator, hard_emitter, in, base, sorters, rel_lits, rel_weights, variable_map);

	for(auto sp = sorters.begin(); sp != sorters.end(); ++sp)
		for(int i = 0; i < sp->size(); i++)
			solver.lockVariable((*sp)[i].variable());

	std::vector<InLiteral> assignment;

	solver.solveStart();
	typename Solver::Result hard_res = Solver::Result::kBreak;
	while(hard_res == Solver::Result::kBreak)
		hard_res = solver.solveStep();

	if(hard_res == Solver::Result::kHardUnsat) {
		std::cout << "s UNSATISFIABLE" << std::endl;
		return;
	}else if(hard_res != Solver::Result::kSat) {
		std::cout << "c Error in SAT solver!" << std::endl;
		std::cout << "o UNKNOWN" << std::endl;
		return;
	}
	std::cout << "c Hard clauses are satisfiable" << std::endl;
	
	std::cout << "o " << countSoftWeight(in) - cur_lb << std::endl;

	for(auto it = variable_map.begin(); it != variable_map.end(); ++it)
		if(solver.litInModel(it->second.oneLiteral())) {
			assignment.push_back(it->first.pos());
		}else{
			assignment.push_back(it->first.neg());
	}
	solver.solveReset();
	
	while(cur_lb + 1 != cur_ub) {
		Weight target = (cur_lb + cur_ub) / 2;
		std::cout << "c Current LB: " << cur_lb << ", UB: " << cur_ub << ", target: " << target << std::endl;
		assert(target > cur_lb);
		assert(target < cur_ub);

		typename Solver::ContextEmitter rhs_emitter(solver);
		maxsatRhs(hard_allocator, rhs_emitter, base, target, sorters);

		solver.solveStart();
		typename Solver::Result res = Solver::Result::kBreak;
		while(res == Solver::Result::kBreak)
			res = solver.solveStep();
		
		if(res == Solver::Result::kSat) {
			std::cout << "c SAT solver returned satisfiable!" << std::endl;
			cur_lb = target;
	
			assignment.clear();
			for(auto it = variable_map.begin(); it != variable_map.end(); ++it)
				if(solver.litInModel(it->second.oneLiteral())) {
					assignment.push_back(it->first.pos());
				}else{
					assignment.push_back(it->first.neg());
				}
			std::cout << "o " << countSoftWeight(in) - cur_lb << std::endl;
		}else if(res == Solver::Result::kSoftUnsat) {
			std::cout << "c SAT solver returned unsatisfiable!" << std::endl;
			cur_ub = target;
		}else{
			std::cout << "c Error in SAT solver!" << std::endl;
			std::cout << "o UNKNOWN" << std::endl;
			return;
		}
		solver.solveReset();
	}
	std::cout << "s OPTIMUM FOUND" << std::endl;

	std::cout << "v";
	std::sort(assignment.begin(), assignment.end(),
		[](InLiteral a, InLiteral b) -> bool {
			return a.var().index() < b.var().index();
		});
	for(auto it = assignment.begin(); it != assignment.end(); ++it)
		std::cout << " " << it->num();
	std::cout << std::endl;
}

} // namespace maxsatuzk

