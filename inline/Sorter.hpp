
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
#include <encodeuzk/mixed-radix.inline.hpp>
#include <encodeuzk/sorting.inline.hpp>
#include <encodeuzk/static.inline.hpp>

namespace maxsatuzk {

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
		encodeuzk::SorterNetwork<typename ClauseEmitter::Literal> &sorters,
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

	sorters = encodeuzk::computeSorterNetwork(allocator, emitter, rel_lits, rel_weights, base);
}

// generates a formula, that is satisfiable iff at least min soft clauses can be satisfied
template<typename VarAllocator, typename ClauseEmitter>
void maxsatRhs(VarAllocator &allocator, ClauseEmitter &emitter,
		std::vector<int> &base, int min,
		encodeuzk::SorterNetwork<typename ClauseEmitter::Literal> &sorters) {
	std::vector<int> rhs = convertBase(min, base);

	typename ClauseEmitter::Literal lit
			= encodeuzk::computeSorterNetworkGe(allocator, emitter, sorters, base, rhs);
	encodeuzk::forceTrue(allocator, emitter, lit);;
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

	encodeuzk::SorterNetwork<typename Solver::Literal> sorters;
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

