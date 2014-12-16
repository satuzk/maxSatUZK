
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

// produces a constraint that is true iff output(s) % n >= lim
template<typename VarAllocator, typename ClauseEmitter>
void makeModGe(VarAllocator &allocator, ClauseEmitter &emitter,
		std::vector<Sorter<typename ClauseEmitter::Literal>> &sorters, int i,
		int n, int lim, typename ClauseEmitter::Literal r) {
	Sorter<typename ClauseEmitter::Literal> &s = sorters[i];
	if(lim == 0) {
		// trivial case 1: every number is >= 0
		encodeuzk::emit(emitter, { r });
		if(debugRhs) {
			std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
			std::cout << r.toNumber() << " <-> true" << std::endl;
		}
		return;
	}else if(s.size() < lim) {
		// trivial case 2: the sorter is not big enough to reach the limit
		encodeuzk::emit(emitter, { r.inverse() });
		if(debugRhs) {
			std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
			std::cout << r.toNumber() << " <-> false" << std::endl;
		}
		return;
	}else if(n <= lim) {
		// trivial case 3: the modulus is not big enough to reach the limit
		encodeuzk::emit(emitter, { r.inverse() });
		if(debugRhs) {
			std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
			std::cout << r.toNumber() << " <-> false" << std::endl;
		}
		return;
	}
	
	if(debugRhs) {
		std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
		std::cout << r.toNumber() << " <-> false OR ..." << std::endl;
	}

	std::vector<typename ClauseEmitter::Literal> temp;
	for(int k = 0; k < s.size(); k += n) {
		if(k + lim - 1 >= s.size())
			break;

		typename ClauseEmitter::Literal x = allocator.allocate().oneLiteral();
		temp.push_back(x);
		
		if(k + n - 1 < s.size()) {
			// each temp literal is equivalent to a conjunction
			encodeuzk::emit(emitter, { s.output(k + lim - 1).inverse(),
					s.output(k + n - 1), x });
			encodeuzk::emit(emitter, { x.inverse(), s.output(k + lim - 1) });
			encodeuzk::emit(emitter, { x.inverse(), s.output(k + n - 1).inverse() });
			if(debugRhs) {
				std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
				std::cout << x.toNumber() << " <-> " << s.output(k + lim - 1).toNumber() << " (out) "
					<< " AND " << -(s.output(k + n - 1)).toNumber() << " (out) "  << std::endl;
			}
		}else{
			// the last temp literal is equivalent to a sorter output
			encodeuzk::emit(emitter, { s.output(k + lim - 1).inverse(), x });
			encodeuzk::emit(emitter, { x.inverse(), s.output(k + lim - 1) });
			if(debugRhs) {
				std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
				std::cout << x.toNumber() << " <-> " << s.output(k + lim - 1).toNumber() << " (out) " << std::endl;
			}
		}

		// each temp literal implies r
		encodeuzk::emit(emitter, { x.inverse(), r });
		if(debugRhs) {
			std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
			std::cout << r.toNumber() << " <-> ... OR " << x.toNumber() << std::endl;
		}
	}

	// r implies at least one temp literal
	temp.push_back(r.inverse());
	emitter.emit(temp.begin(), temp.end());
}

// produces a constraint that is true iff output(s) >= lim
template<typename VarAllocator, typename ClauseEmitter>
void makeGe(VarAllocator &allocator, ClauseEmitter &emitter,
		std::vector<Sorter<typename ClauseEmitter::Literal>> &sorters, int i,
		int lim, typename ClauseEmitter::Literal r) {
	Sorter<typename ClauseEmitter::Literal> &s = sorters[i];
	if(lim == 0) {
		/* trivial case 1: every number is >= 0 */
		encodeuzk::emit(emitter, { r });
		if(debugRhs) {
			std::cout << "s[" << i << "] >= " << lim << ": ";
			std::cout << r.toNumber() << " <-> true" << std::endl;
		}
	}else if(s.size() < lim) {
		// trivial case 2: the sorter is not big enough to reach the limit
		encodeuzk::emit(emitter, { r.inverse() });
		if(debugRhs) {
			std::cout << "s[" << i << "] >= " << lim << ": ";
			std::cout << r.toNumber() << " <-> false" << std::endl;
		}
	}else{
		encodeuzk::emit(emitter, { s.output(lim - 1).inverse(), r });
		encodeuzk::emit(emitter, { r.inverse(), s.output(lim - 1) });
		if(debugRhs) {
			std::cout << "s[" << i << "] >= " << lim << ": ";
			std::cout << r.toNumber() << " <-> " << s.output(lim - 1).inverse().toNumber() << " (out)" << std::endl;
		}
	}
}

// generates the constraint (sorters >= rhs)
template<typename VarAllocator, typename ClauseEmitter>
void buildGeRhs(VarAllocator &allocator, ClauseEmitter &emitter,
		int i, NumberSeq rhs, NumberSeq base,
		typename ClauseEmitter::Literal r, std::vector<Sorter<typename ClauseEmitter::Literal>> &sorters) {
	if(i == 0) {
		// trivial case: every number is >= 0
		encodeuzk::emit(emitter, { r });
		if(debugRhs)
			std::cout << r.toNumber() << " <-> true" << std::endl;
	}else{
		i--;
		Sorter<typename ClauseEmitter::Literal> &s = sorters[i];
		
		typename ClauseEmitter::Literal p = allocator.allocate().oneLiteral();
		buildGeRhs(allocator, emitter, i, rhs, base, p, sorters);
		
		typename ClauseEmitter::Literal gt = allocator.allocate().oneLiteral();
		typename ClauseEmitter::Literal ge = allocator.allocate().oneLiteral();
		if(i == base.length() - 1) {
			makeGe(allocator, emitter, sorters, i, rhs[i] + 1, gt);
			makeGe(allocator, emitter, sorters, i, rhs[i], ge);
		}else{
			makeModGe(allocator, emitter, sorters, i, base[i+1], rhs[i] + 1, gt);
			makeModGe(allocator, emitter, sorters, i, base[i+1], rhs[i], ge);
		}
		
		// q <-> ge & p
		typename ClauseEmitter::Literal q = allocator.allocate().oneLiteral();
		encodeuzk::emit(emitter, { q.inverse(), ge });
		encodeuzk::emit(emitter, { q.inverse(), p });
		encodeuzk::emit(emitter, { ge.inverse(), p.inverse(), q });
		if(debugRhs)
			std::cout << q.toNumber() << " <-> " << ge.toNumber() << " (ge) AND " << p.toNumber() << " (prev) " << std::endl;

		// r <-> gt | q
		encodeuzk::emit(emitter, { r.inverse(), gt, q });
		encodeuzk::emit(emitter, { gt.inverse(), r });
		encodeuzk::emit(emitter, { q.inverse(), r });
		if(debugRhs)
			std::cout << r.toNumber() << " <-> " << gt.toNumber() << " (gt) OR " << q.toNumber() << std::endl;
	}
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

	typename ClauseEmitter::Variable v = allocator.allocate();
	buildGeRhs(allocator, emitter, sorters.size(), rhs, base, v.oneLiteral(), sorters);
	if(debugRhs)
		std::cout << "objective: " << v.oneLiteral().toNumber() << std::endl;

	encodeuzk::emit(emitter, { v.oneLiteral() });
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

