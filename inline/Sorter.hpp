
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

uint64_t statComparators = 0;

template<typename Literal>
class Sorter {
public:
	Sorter(int size) {
		p_inputs.resize(size, Literal::illegalLit());
		p_outputs.resize(size, Literal::illegalLit());
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

template<typename VarAllocator, typename ClauseEmitter>
void buildSorters(VarAllocator &allocator, ClauseEmitter &emitter,
		std::vector<typename ClauseEmitter::Literal> &lits,
		std::vector<Weight> &weights,
		NumberSeq base,
		std::vector<Sorter<typename ClauseEmitter::Literal>> &sorters) {
	// we need a literal that is always zero to simplify sorting
	typename ClauseEmitter::Literal null_lit = allocator.allocate().oneLiteral();
	encodeuzk::emit(emitter, { null_lit.inverse() });

	for(int k = 0; k < base.length(); k++) {
		std::vector<typename ClauseEmitter::Literal> ins;

		// add carry bits from previous sorter as input
		if(k > 0) {
			Sorter<typename ClauseEmitter::Literal> &p = sorters.back();
			for(int j = base[k] - 1; j < p.size(); j += base[k]) {
				ins.push_back(p.output(j));
				if(debugSorters)
					std::cout << "carry from " << (k - 1) << " to " << k << std::endl;
			}
		}

		for(int i = 0; i < lits.size(); i++) {
			NumberSeq weight = convertBase(weights[i], base);
			for(int j = 0; j < weight[k]; j++) {
				ins.push_back(lits[i]);
				if(debugSorters)
					std::cout << "s[" << k << "] weight" << std::endl;
			}
		}

		std::vector<typename ClauseEmitter::Literal> outs
				= encodeuzk::computePwSort(allocator, emitter, ins, null_lit);

//		std::cout << "Creating sorter of size " << ins.size() << std::endl;
		Sorter<typename ClauseEmitter::Literal> s(ins.size());
		for(int i = 0; i < ins.size(); ++i) {
			s.input(i) = ins[i];
			if(debugSorters)
				std::cout << "in[" << k << "][" << i << "]: " << ins[i].toNumber() << std::endl;
		}
		for(int i = 0; i < ins.size(); ++i) {
			s.output(i) = outs[i];
			if(debugSorters)
				std::cout << "out[" << k << "][" << i << "]: " << outs[i].toNumber() << std::endl;
		}
		std::cout << "c Size of sorter " << k << ": " << s.size() << std::endl;
		sorters.push_back(s);
	}
}

// produces a constraint that is true iff sorter >= target
template<typename VarAllocator, typename ClauseEmitter>
typename ClauseEmitter::Literal computeSorterGe(VarAllocator &allocator, ClauseEmitter &emitter,
		Sorter<typename ClauseEmitter::Literal> &sorter, int target) {
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

	return sorter.output(target - 1);
}

// produces a constraint that is true iff sorter % divisor >= target
template<typename VarAllocator, typename ClauseEmitter>
typename ClauseEmitter::Literal computeSorterRemainderGe(VarAllocator &allocator, ClauseEmitter &emitter,
		Sorter<typename ClauseEmitter::Literal> &sorter, int divisor, int target) {
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
					sorter.output(k + target - 1), sorter.output(k + divisor - 1).inverse()));
		}else{
			disjunction.push_back(sorter.output(k + target - 1));
		}
	}
	
	return encodeuzk::computeOrN(allocator, emitter,
			disjunction.begin(), disjunction.end());
}

// generates the constraint (sorters >= rhs)
template<typename VarAllocator, typename ClauseEmitter>
typename ClauseEmitter::Literal computeSorterNetworkGe(VarAllocator &allocator, ClauseEmitter &emitter,
		std::vector<Sorter<typename ClauseEmitter::Literal>> &sorters,
		NumberSeq base, NumberSeq rhs, int i) {

	if(i == 0) {
		typename ClauseEmitter::Variable r = allocator.allocate();

		// trivial case: every number is >= 0
		encodeuzk::emit(emitter, { r.oneLiteral() });
		if(debugRhs)
			std::cout << r.oneLiteral().toNumber() << " <-> true" << std::endl;

		return r.oneLiteral();
	}

	i--;
	Sorter<typename ClauseEmitter::Literal> &s = sorters[i];
	
	typename ClauseEmitter::Literal  p
			= computeSorterNetworkGe(allocator, emitter, sorters, base, rhs, i);
	
	typename ClauseEmitter::Literal gt;
	typename ClauseEmitter::Literal ge;
	if(i == base.length() - 1) {
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
		std::vector<Sorter<typename ClauseEmitter::Literal>> &sorters,
		NumberSeq base, NumberSeq rhs) {
	return computeSorterNetworkGe(allocator, emitter, sorters, base, rhs,
			sorters.size());
}

template<typename VarAllocator, typename ClauseEmitter>
void maxsatHard(VarAllocator &allocator, ClauseEmitter &emitter,
		InClauseSpace &in,
		std::unordered_map<InVariable, typename ClauseEmitter::Variable, VarHashFunc> &variable_map) {
	// copy all variables from the original formula
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		InClauseRef clause = *i;
		for(int j = 0; j < clause.length(); j++) {
			InVariable v = clause.getLiteral(j).var();
			if(variable_map.find(v) == variable_map.end())
				variable_map[v] = allocator.allocate();
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
		InClauseSpace &in, NumberSeq &base,
		std::vector<Sorter<typename ClauseEmitter::Literal>> &sorters,
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

	buildSorters(allocator, emitter, rel_lits, rel_weights, base, sorters);
}

// generates a formula, that is satisfiable iff at least min soft clauses can be satisfied
template<typename VarAllocator, typename ClauseEmitter>
void maxsatRhs(VarAllocator &allocator, ClauseEmitter &emitter,
		NumberSeq &base, int min,
		std::vector<Sorter<typename ClauseEmitter::Literal>> &sorters) {
	NumberSeq rhs = convertBase(min, base);
	if(debugRhs)
		for(int i = 0; i < rhs.length(); i++)
			std::cout << "rhs[" << i << "]: " << rhs[i] << std::endl;

	typename ClauseEmitter::Literal r
			= computeSorterNetworkGe(allocator, emitter, sorters, base, rhs);
	if(debugRhs)
		std::cout << "objective: " << r.toNumber() << std::endl;

	encodeuzk::forceTrue(allocator, emitter, r);;
}

template<typename Solver>
void solve(InClauseSpace &in, NumberSeq &base, Weight initial_lb, Weight initial_ub) {
	Weight cur_lb = initial_lb; // formula is sat for v <= lb
	Weight cur_ub = initial_ub; // formula is unsat for v >= ub

	Solver solver(0);
	typename Solver::HardAllocator hard_allocator(solver);
	typename Solver::HardEmitter hard_emitter(solver);

	std::unordered_map<InVariable, typename Solver::Variable,
			VarHashFunc> variable_map;
	maxsatHard(hard_allocator, hard_emitter, in, variable_map);

	std::vector<Sorter<typename Solver::Literal>> sorters;
	std::vector<typename Solver::Literal> rel_lits;
	std::vector<Weight> rel_weights;
	maxsatLhs(hard_allocator, hard_emitter, in, base, sorters, rel_lits, rel_weights, variable_map);

	for(auto sp = sorters.begin(); sp != sorters.end(); ++sp)
		for(int i = 0; i < sp->size(); i++)
			solver.lockVariable(sp->output(i).variable());

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

