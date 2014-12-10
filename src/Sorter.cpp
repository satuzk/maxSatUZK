
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
#include "../include/Solver.hpp"
#include "../include/SolverUzk.hpp"
#include "../include/Sorter.hpp"
#include "../inline/Formula.hpp"

namespace maxsatuzk {

static const bool debugSorters = false;
static const bool debugRhs = false;

uint64_t statComparators = 0;

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

void comparator(Literal x1, Literal x2, Literal y1, Literal y2,
		ClauseSpace &f) {
//	std::cout << "comparator(" << x1.num() << ", " << x2.num() << ", " << y1.num() << ", "
//		<< y2.num() << ")" << std::endl;
	// these clauses represent min(x1, x2) <= y1, max(x1, x2) <= y2
	ClauseRef c1 = f.allocate(2);
	ClauseRef c2 = f.allocate(2);
	ClauseRef c3 = f.allocate(3);
	c1.setLiteral(0, -x1);
	c1.setLiteral(1, y1);
	c2.setLiteral(0, -x2);
	c2.setLiteral(1, y1);
	c3.setLiteral(0, -x1);
	c3.setLiteral(1, -x2);
	c3.setLiteral(2, y2);

	// these clauses represent min(x1, x2) >= y1, max(x1, x2) >= y2
	ClauseRef r1 = f.allocate(2);
	ClauseRef r2 = f.allocate(2);
	ClauseRef r3 = f.allocate(3);
	r1.setLiteral(0, -y2);
	r1.setLiteral(1, x1);
	r2.setLiteral(0, -y2);
	r2.setLiteral(1, x2);
	r3.setLiteral(0, -y1);
	r3.setLiteral(1, x1);
	r3.setLiteral(2, x2);
}

template<typename Literal>
void pwSplit(VarAllocator &va, std::vector<Literal> &ins,
		std::vector<Literal> &outs_a, std::vector<Literal> &outs_b,
		Literal null_lit, ClauseSpace &f) {
	assert(ins.size() % 2 == 0);

	for(int i = 0; i < ins.size() / 2; i++) {
		Literal a = va.alloc().pos();
		Literal b = va.alloc().pos();
		outs_a.push_back(a);
		outs_b.push_back(b);
		comparator(ins[2 * i], ins[2 * i + 1], a, b, f);
	}
	statComparators += ins.size() / 2;
}

void pwMerge(VarAllocator &va, std::vector<Literal> &a, std::vector<Literal> &b,
		std::vector<Literal> &outs, Literal null_lit, ClauseSpace &f) {
	assert(a.size() > 0);
	assert(a.size() == b.size());
	if(a.size() == 1) {
		outs.push_back(a.front());
		outs.push_back(b.front());
		return;
	}else if(a.size() % 2 == 1) {
		std::vector<Literal> new_a;
		std::vector<Literal> new_b;
		for(int i = 0; i < a.size(); i++)
			new_a.push_back(a[i]);
		for(int i = 0; i < a.size(); i++)
			new_b.push_back(b[i]);
		new_a.push_back(null_lit);
		new_b.push_back(null_lit);
		pwMerge(va, new_a, new_b, outs, null_lit, f);
		outs.pop_back();
		outs.pop_back();
		return;
	}
	
	assert(a.size() % 2 == 0);
	
	std::vector<Literal> even_a;
	std::vector<Literal> even_b;
	for(unsigned int i = 0; i < a.size(); i += 2) {
		even_a.push_back(a[i]);
		even_b.push_back(b[i]);
	}

	std::vector<Literal> odd_a;
	std::vector<Literal> odd_b;
	for(unsigned int i = 1; i < a.size(); i += 2) {
		odd_a.push_back(a[i]);
		odd_b.push_back(b[i]);
	}

	std::vector<Literal> even_temps;
	std::vector<Literal> odd_temps;
	pwMerge(va, even_a, even_b, even_temps, null_lit, f);
	pwMerge(va, odd_a, odd_b, odd_temps, null_lit, f);
	assert(even_temps.size() == a.size());
	assert(odd_temps.size() == a.size());

	// number of bits that actually have to be merged
	int n_merge = a.size() - 1;

	outs.push_back(even_temps.front());
	for(int i = 0; i < n_merge; i++)
		outs.push_back(va.alloc().pos());
	for(int i = 0; i < n_merge; i++)
		outs.push_back(va.alloc().pos());
	outs.push_back(odd_temps.back());
	assert(outs.size() == 2 * a.size());

	for(unsigned int i = 0; i < n_merge; i++)
		comparator(even_temps[i + 1], odd_temps[i],
				outs[2 * i + 1], outs[2 * i + 2], f);
	statComparators += n_merge;

	// additional clauses to improve propagation
	for(int i = 0; i < outs.size() - 1; i++) {
		ClauseRef c = f.allocate(2);
		c.setLiteral(0, outs[i]);
		c.setLiteral(1, -outs[i + 1]);
	}
}

void pwSort(VarAllocator &va, std::vector<Literal> &ins,
		std::vector<Literal> &outs, Literal null_lit, ClauseSpace &f) {
	if(ins.size() == 0)
		return;
	if(ins.size() == 1) {
		outs.push_back(ins[0]);
		return;
	}
	if(ins.size() % 2 == 1) {
		std::vector<Literal> new_ins;
		for(int i = 0; i < ins.size(); i++)
			new_ins.push_back(ins[i]);
		new_ins.push_back(null_lit);
		pwSort(va, new_ins, outs, null_lit, f);
		outs.pop_back();
		return;
	}
	
	std::vector<Literal> ins_a;
	std::vector<Literal> ins_b;
	pwSplit(va, ins, ins_a, ins_b, null_lit, f);
		
	std::vector<Literal> outs_a;
	std::vector<Literal> outs_b;
	pwSort(va, ins_a, outs_a, null_lit, f);
	pwSort(va, ins_b, outs_b, null_lit, f);
	pwMerge(va, outs_a, outs_b, outs, null_lit, f);
	assert(outs.size() == ins.size());
}

void buildSorters(VarAllocator &va,
		std::vector<Literal> &lits,
		std::vector<Weight> &weights,
		NumberSeq base,
		std::vector<Sorter> &sorters,
		ClauseSpace &f) {
	// we need a literal that is always zero to simplify sorting
	Literal null_lit = va.alloc().pos();
	ClauseRef c = f.allocate(1);
	c.setLiteral(0, -null_lit);

	for(int k = 0; k < base.length(); k++) {
		std::vector<Literal> ins;

		// add carry bits from previous sorter as input
		if(k > 0) {
			Sorter &p = sorters.back();
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

		std::vector<Literal> outs;
		pwSort(va, ins, outs, null_lit, f);

//		std::cout << "Creating sorter of size " << ins.size() << std::endl;
		Sorter s(ins.size());
		for(int i = 0; i < ins.size(); ++i) {
			s.input(i) = ins[i];
			if(debugSorters)
				std::cout << "in[" << k << "][" << i << "]: " << ins[i].num() << std::endl;
		}
		for(int i = 0; i < ins.size(); ++i) {
			s.output(i) = outs[i];
			if(debugSorters)
				std::cout << "out[" << k << "][" << i << "]: " << outs[i].num() << std::endl;
		}
		std::cout << "c Size of sorter " << k << ": " << s.size() << std::endl;
		sorters.push_back(s);
	}
}

// produces a constraint that is true iff output(s) % n >= lim
void makeModGe(VarAllocator &va, std::vector<Sorter> &sorters, int i,
		int n, int lim, Literal r, ClauseSpace &f) {
	Sorter &s = sorters[i];
	if(lim == 0) {
		// trivial case 1: every number is >= 0
		ClauseRef c = f.allocate(1);
		c.setLiteral(0, r);
		if(debugRhs) {
			std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
			std::cout << r.num() << " <-> true" << std::endl;
		}
		return;
	}else if(s.size() < lim) {
		// trivial case 2: the sorter is not big enough to reach the limit
		ClauseRef c = f.allocate(1);
		c.setLiteral(0, -r);
		if(debugRhs) {
			std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
			std::cout << r.num() << " <-> false" << std::endl;
		}
		return;
	}else if(n <= lim) {
		// trivial case 3: the modulus is not big enough to reach the limit
		ClauseRef c = f.allocate(1);
		c.setLiteral(0, -r);
		if(debugRhs) {
			std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
			std::cout << r.num() << " <-> false" << std::endl;
		}
		return;
	}
	
	if(debugRhs) {
		std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
		std::cout << r.num() << " <-> false OR ..." << std::endl;
	}

	std::vector<Literal> temp;
	for(int k = 0; k < s.size(); k += n) {
		if(k + lim - 1 >= s.size())
			break;

		Literal x = va.alloc().pos();
		temp.push_back(x);
		
		if(k + n - 1 < s.size()) {
			// each temp literal is equivalent to a conjunction
			ClauseRef c1 = f.allocate(3);
			ClauseRef c2 = f.allocate(2);
			ClauseRef c3 = f.allocate(2);
			c1.setLiteral(0, -s.output(k + lim - 1));
			c1.setLiteral(1, s.output(k + n - 1));
			c1.setLiteral(2, x);
			c2.setLiteral(0, -x);
			c2.setLiteral(1, s.output(k + lim - 1));
			c3.setLiteral(0, -x);
			c3.setLiteral(1, -s.output(k + n - 1));
			if(debugRhs) {
				std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
				std::cout << x.num() << " <-> " << s.output(k + lim - 1).num() << " (out) "
					<< " AND " << -(s.output(k + n - 1)).num() << " (out) "  << std::endl;
			}
		}else{
			// the last temp literal is equivalent to a sorter output
			ClauseRef c1 = f.allocate(2);
			ClauseRef c2 = f.allocate(2);
			c1.setLiteral(0, -s.output(k + lim - 1));
			c1.setLiteral(1, x);
			c2.setLiteral(0, -x);
			c2.setLiteral(1, s.output(k + lim - 1));
			if(debugRhs) {
				std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
				std::cout << x.num() << " <-> " << s.output(k + lim - 1).num() << " (out) " << std::endl;
			}
		}

		// each temp literal implies r
		ClauseRef c = f.allocate(2);
		c.setLiteral(0, -x);
		c.setLiteral(1, r);
		if(debugRhs) {
			std::cout << "s[" << i << "] mod " << n << " >= " << lim << ": ";
			std::cout << r.num() << " <-> ... OR " << x.num() << std::endl;
		}
	}

	// r implies at least one temp literal
	ClauseRef c = f.allocate(temp.size() + 1);
	c.setLiteral(0, -r);
	for(int j = 0; j < temp.size(); j++) {
		c.setLiteral(j + 1, temp[j]);
	}
}

// produces a constraint that is true iff output(s) >= lim
void makeGe(VarAllocator &va, std::vector<Sorter> &sorters, int i,
		int lim, Literal r, ClauseSpace &f) {
	Sorter &s = sorters[i];
	if(lim == 0) {
		/* trivial case 1: every number is >= 0 */
		ClauseRef c = f.allocate(1);
		c.setLiteral(0, r);
		if(debugRhs) {
			std::cout << "s[" << i << "] >= " << lim << ": ";
			std::cout << r.num() << " <-> true" << std::endl;
		}
	}else if(s.size() < lim) {
		// trivial case 2: the sorter is not big enough to reach the limit
		ClauseRef c = f.allocate(1);
		c.setLiteral(0, -r);
		if(debugRhs) {
			std::cout << "s[" << i << "] >= " << lim << ": ";
			std::cout << r.num() << " <-> false" << std::endl;
		}
	}else{
		ClauseRef c1 = f.allocate(2);
		ClauseRef c2 = f.allocate(2);
		c1.setLiteral(0, -s.output(lim - 1));
		c1.setLiteral(1, r);
		c2.setLiteral(0, -r);
		c2.setLiteral(1, s.output(lim - 1));
		if(debugRhs) {
			std::cout << "s[" << i << "] >= " << lim << ": ";
			std::cout << r.num() << " <-> " << (-s.output(lim - 1)).num() << " (out)" << std::endl;
		}
	}
}

// generates the constraint (sorters >= rhs)
void buildGeRhs(VarAllocator &va,
		int i, NumberSeq rhs, NumberSeq base,
		Literal r, ClauseSpace &f,
		std::vector<Sorter> &sorters) {
	if(i == 0) {
		// trivial case: every number is >= 0
		ClauseRef c = f.allocate(1);
		c.setLiteral(0, r);
		if(debugRhs)
			std::cout << r.num() << " <-> true" << std::endl;
	}else{
		i--;
		Sorter &s = sorters[i];
		
		Literal p = va.alloc().pos();
		buildGeRhs(va, i, rhs, base, p, f, sorters);
		
		Literal gt = va.alloc().pos();
		Literal ge = va.alloc().pos();
		if(i == base.length() - 1) {
			makeGe(va, sorters, i, rhs[i] + 1, gt, f);
			makeGe(va, sorters, i, rhs[i], ge, f);
		}else{
			makeModGe(va, sorters, i, base[i+1], rhs[i] + 1, gt, f);
			makeModGe(va, sorters, i, base[i+1], rhs[i], ge, f);
		}
		
		// q <-> ge & p
		Literal q = va.alloc().pos();
		ClauseRef c1 = f.allocate(2);
		ClauseRef c2 = f.allocate(2);
		ClauseRef c3 = f.allocate(3);
		c1.setLiteral(0, -q);
		c1.setLiteral(1, ge);
		c2.setLiteral(0, -q);
		c2.setLiteral(1, p);
		c3.setLiteral(0, -ge);
		c3.setLiteral(1, -p);
		c3.setLiteral(2, q);
		if(debugRhs)
			std::cout << q.num() << " <-> " << ge.num() << " (ge) AND " << p.num() << " (prev) " << std::endl;

		// r <-> gt | q
		ClauseRef c4 = f.allocate(3);
		ClauseRef c5 = f.allocate(2);
		ClauseRef c6 = f.allocate(2);
		c4.setLiteral(0, -r);
		c4.setLiteral(1, gt);
		c4.setLiteral(2, q);
		c5.setLiteral(0, -gt);
		c5.setLiteral(1, r);
		c6.setLiteral(0, -q);
		c6.setLiteral(1, r);
		if(debugRhs)
			std::cout << r.num() << " <-> " << gt.num() << " (gt) OR " << q.num() << std::endl;
	}
}

Weight countSoftWeight(ClauseSpace &f) {
	Weight count = 0;
	for(auto i = f.refsBegin(); i != f.refsEnd(); ++i)
		if((*i).getWeight() != kHardWeight)
			count += (*i).getWeight();
	return count;
}

void maxsatHard(VarAllocator &va, ClauseSpace &in, ClauseSpace &res,
		std::unordered_map<Variable, Variable, VarHashFunc> &variable_map) {
	// copy all variables from the original formula
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		ClauseRef clause = *i;
		for(int j = 0; j < clause.length(); j++) {
			Variable v = clause.getLiteral(j).var();
			if(variable_map.find(v) == variable_map.end())
				variable_map[v] = va.alloc();
		}
	}
	
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		ClauseRef in_clause = *i;
		if(in_clause.getWeight() != kHardWeight)
			continue;

		int in_length = in_clause.length();
		ClauseRef tf_clause = res.allocate(in_length);
		
		// copy the literals from the original clause
		for(int j = 0; j < in_clause.length(); j++) {
			Literal in_literal = in_clause.getLiteral(j);
			Variable tf_variable = variable_map[in_literal.var()];
			tf_clause.setLiteral(j, in_literal.sign() > 0 ? tf_variable.pos() : tf_variable.neg());
		}
	}
}

void maxsatLhs(VarAllocator &va, ClauseSpace &in,
		NumberSeq &base, ClauseSpace &res,
		std::vector<Sorter> &sorters,
		std::vector<Literal> &rel_lits,
		std::vector<Weight> &rel_weights,
		std::unordered_map<Variable, Variable, VarHashFunc> &variable_map) {
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		ClauseRef in_clause = *i;
		if(in_clause.getWeight() == kHardWeight)
			continue;

		int in_length = in_clause.length();
		ClauseRef tf_clause = res.allocate(in_length + 1);
		
		// copy the literals from the original clause
		for(int j = 0; j < in_clause.length(); j++) {
			Literal in_literal = in_clause.getLiteral(j);
			Variable tf_variable = variable_map[in_literal.var()];
			tf_clause.setLiteral(j, in_literal.sign() > 0 ? tf_variable.pos() : tf_variable.neg());
		}

		// add a relaxation literal
		Literal rel_lit = va.alloc().pos();
		tf_clause.setLiteral(in_length, rel_lit);
		rel_lits.push_back(-rel_lit);
		rel_weights.push_back(in_clause.getWeight());
	}

	buildSorters(va, rel_lits, rel_weights, base, sorters, res);
}

// generates a formula, that is satisfiable iff at least min soft clauses can be satisfied
void maxsatRhs(VarAllocator &va, NumberSeq &base, int min,
		ClauseSpace &res, std::vector<Sorter> &sorters) {
	NumberSeq rhs = convertBase(min, base);
	if(debugRhs)
		for(int i = 0; i < rhs.length(); i++)
			std::cout << "rhs[" << i << "]: " << rhs[i] << std::endl;

	Variable v = va.alloc();
	buildGeRhs(va, sorters.size(), rhs, base, v.pos(), res, sorters);
	if(debugRhs)
		std::cout << "objective: " << v.pos().num() << std::endl;

	ClauseRef c = res.allocate(1);
	c.setLiteral(0, v.pos());
}

void solve(ClauseSpace &in, NumberSeq &base, Weight initial_lb, Weight initial_ub) {
	Weight cur_lb = initial_lb; // formula is sat for v <= lb
	Weight cur_ub = initial_ub; // formula is unsat for v >= ub
	
	VarAllocator va_hard;
	ClauseSpace f_lhs;
	std::unordered_map<Variable, Variable, VarHashFunc> variable_map;
	maxsatHard(va_hard, in, f_lhs, variable_map);

	VarAllocator va_lhs(va_hard.next_var(), std::numeric_limits<int>::max());
	std::vector<Sorter> sorters;
	std::vector<Literal> rel_lits;
	std::vector<Weight> rel_weights;
	maxsatLhs(va_lhs, in, base, f_lhs, sorters, rel_lits, rel_weights, variable_map);

	// estimate the number of variables required to express the rhs. this is an upper bound.
	int estimate = 1; // variable equivalent to the constraint
	estimate += 4 * sorters.size(); // four variables allocated in buildGeRhs()
	// variables allocated in makeModGe(). note that makeModGe is called two times.
	for(int i = 0; i < base.length() - 1; i++)
		estimate += 2 * (sorters[i].size() / base[i + 1] + 1);
	std::cout << "c Reserving " << estimate << " vars for the RHS" << std::endl;

	SolverUzk solver(0);
	solver.reserveVars(va_lhs.next_var() + estimate);
	solver.setupLhs(f_lhs);

	for(auto sp = sorters.begin(); sp != sorters.end(); ++sp)
		for(int i = 0; i < sp->size(); i++)
			solver.lockVariable(sp->output(i).var());
	
	// TODO: we sould not have to do this!
	for(int i = va_lhs.next_var(); i < va_lhs.next_var() + estimate; i++)
		solver.lockVariable(Variable::from_index(i));

	std::cout << "o " << countSoftWeight(in) - cur_lb << std::endl;

	std::vector<Literal> assignment;

	solver.solveStart();
	Solver::Result hard_res = Solver::Result::kBreak;
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

	for(auto it = variable_map.begin(); it != variable_map.end(); ++it)
		if(solver.litInModel(it->second.pos())) {
			assignment.push_back(it->first.pos());
		}else{
			assignment.push_back(it->first.neg());
		}
	solver.solveReset();
	
	while(cur_lb + 1 != cur_ub) {
		Weight v = (cur_lb + cur_ub) / 2;
		std::cout << "c Current LB: " << cur_lb << ", UB: " << cur_ub << ", v: " << v << std::endl;
		assert(v > cur_lb);
		assert(v < cur_ub);

		ClauseSpace f_rhs;
		VarAllocator va_rhs(va_lhs.next_var(), va_lhs.next_var() + estimate);
		maxsatRhs(va_rhs, base, v, f_rhs, sorters);

		solver.updateRhs(f_rhs);
		solver.solveStart();
		Solver::Result res = Solver::Result::kBreak;
		while(res == Solver::Result::kBreak)
			res = solver.solveStep();
		
		if(res == Solver::Result::kSat) {
			std::cout << "c SAT solver returned satisfiable!" << std::endl;
			cur_lb = v;
	
			assignment.clear();
			for(auto it = variable_map.begin(); it != variable_map.end(); ++it)
				if(solver.litInModel(it->second.pos())) {
					assignment.push_back(it->first.pos());
				}else{
					assignment.push_back(it->first.neg());
				}
			std::cout << "o " << countSoftWeight(in) - cur_lb << std::endl;
		}else if(res == Solver::Result::kSoftUnsat) {
			std::cout << "c SAT solver returned unsatisfiable!" << std::endl;
			cur_ub = v;
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
		[](Literal a, Literal b) -> bool {
			return a.var().index() < b.var().index();
		});
	for(auto it = assignment.begin(); it != assignment.end(); ++it)
		std::cout << " " << it->num();
	std::cout << std::endl;
}

} // namespace maxsatuzk

