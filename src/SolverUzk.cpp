
#include <stdexcept>
#include "../include/Formula.hpp"
#include "../include/SolverUzk.hpp"
#include "../inline/Formula.hpp"

namespace maxsatuzk {

void SolverUzk::lockVariable(Variable variable) {
	p_config.lockVariable(variable);
}

void SolverUzk::solveStart() {
	p_config.inputFinish();
	
	if(p_numRuns == 0) {
		std::cout << "c [      ]  before pp:" << std::endl;
		std::cout << "c [      ]     variables: " << p_config.p_varConfig.presentCount()
				<< ", clauses: " << p_config.p_clauseConfig.numPresent()
				<< ", watch lists: " << p_config.p_varConfig.watchOverallSize() << std::endl;	
		
		p_config.occurConstruct();
		
		UnhideRunStats stat_unhide;
		bceEliminateAll(p_config);
		vecdEliminateAll(p_config);
		selfsubEliminateAll(p_config);
		for(int i = 0; i < 5; i++)
			unhideEliminateAll(p_config, true, stat_unhide);
		
		p_config.occurDestruct();
		
		std::cout << "c [      ]  after pp:" << std::endl;
		std::cout << "c [      ]     variables: " << p_config.p_varConfig.presentCount()
				<< ", clauses: " << p_config.p_clauseConfig.numPresent()
				<< ", watch lists: " << p_config.p_varConfig.watchOverallSize() << std::endl;	
	}

	p_config.start();
}

SolverUzk::Result SolverUzk::solveStep() {
	satuzk::TotalSearchStats stats;
	satuzk::SolveState result = satuzk::kStateUnknown;
	satuzk::search(p_config, result, stats);

	if(result == satuzk::kStateSatisfied) {
		return Result::kSat;
	}else if(result == satuzk::kStateAssumptionFail) {
		return Result::kSoftUnsat;
	}else if(result == satuzk::kStateUnsatisfiable) {
		return Result::kHardUnsat;
	}else if(result == satuzk::kStateBreak) {
		return Result::kBreak;
	}else return Result::kError;
}

void SolverUzk::solveReset() {
	p_config.reset();
	p_numRuns++;
}

bool SolverUzk::litInModel(Literal literal) {
		//FIXME: use modelGetLiteral
	return p_config.litTrue(literal);
}

} // namespace maxsatuzk

