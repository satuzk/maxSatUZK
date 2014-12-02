
#include <stdexcept>
#include "../include/Formula.hpp"
#include "../include/Solver.hpp"
#include "../include/SolverUzk.hpp"
#include "../inline/Formula.hpp"

namespace maxsatuzk {

void SolverUzk::reserveVars(int num_vars) {
	p_config.varReserve(num_vars + 1);
	for(long i = 0; i < num_vars; ++i)	
		p_config.varAlloc();

	std::cout << "c Number of generated vars: " << num_vars << std::endl;	

	// first unused variable becomes the context variable
	p_contextVar = p_config.varAlloc();
	p_config.assumptionEnable(p_contextVar.zeroLiteral());
}

void SolverUzk::setupLhs(ClauseSpace &f) {
	for(auto it = f.refsBegin(); it != f.refsEnd(); ++it) {
		std::vector<typename Config::Literal> out_clause;
		for(int i = 0; i < (*it).length(); i++) {
			Literal lit = (*it).getLiteral(i);
			int var_index = lit.var().index() - 1;
			auto solver_var = Config::Variable::fromIndex(var_index);
			if(lit.sign() > 0) {
				out_clause.push_back(solver_var.oneLiteral());
			}else{
				out_clause.push_back(solver_var.zeroLiteral());
			}
		}

		p_config.inputClause(out_clause.size(),
			out_clause.begin(), out_clause.end());
	}
}

void SolverUzk::updateRhs(ClauseSpace &f) {
	p_config.expellContaining(p_contextVar.oneLiteral());
	
	for(auto it = f.refsBegin(); it != f.refsEnd(); ++it) {
		std::vector<typename Config::Literal> out_clause;
		out_clause.push_back(p_contextVar.oneLiteral());
		for(int i = 0; i < (*it).length(); i++) {
			Literal lit = (*it).getLiteral(i);
			int var_index = lit.var().index() - 1;
			auto solver_var = Config::Variable::fromIndex(var_index);
			if(lit.sign() > 0) {
				out_clause.push_back(solver_var.oneLiteral());
			}else{
				out_clause.push_back(solver_var.zeroLiteral());
			}
		}

		p_config.inputClause(out_clause.size(),
			out_clause.begin(), out_clause.end());
	}
}

void SolverUzk::solveStart() {
	p_config.inputFinish();

	p_config.start();
}

Solver::Result SolverUzk::solveStep() {
	satuzk::TotalSearchStats stats;
	satuzk::SolveState result = satuzk::kStateUnknown;
	satuzk::search(p_config, result, stats);

	if(result == satuzk::kStateSatisfied) {
		return Solver::Result::kSat;
	}else if(result == satuzk::kStateAssumptionFail) {
		return Solver::Result::kSoftUnsat;
	}else if(result == satuzk::kStateUnsatisfiable) {
		return Solver::Result::kHardUnsat;
	}else if(result == satuzk::kStateBreak) {
		return Solver::Result::kBreak;
	}else return Solver::Result::kError;
}

void SolverUzk::solveReset() {
	p_config.reset();
}

bool SolverUzk::litInModel(Literal literal) {
	int var_index = literal.var().index() - 1;
	auto solver_var = Config::Variable::fromIndex(var_index);
	if(literal.sign() > 0) {
		//FIXME: use modelGetLiteral
		return p_config.litTrue(solver_var.oneLiteral());
	}else{
		return p_config.litTrue(solver_var.zeroLiteral());
	}
}

} // namespace maxsatuzk

