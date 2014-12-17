
#include <unistd.h>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <cstring>
#include <csignal>
#include <cassert>

#include "../satuzk/inline/sys/Debug.hpp"
#include "../satuzk/inline/sys/Reporter.hpp"
#include "../satuzk/inline/sys/Performance.hpp"
#include "../satuzk/inline/util/BulkAlloc.hpp"
#include "../satuzk/inline/util/BinaryHeap.hpp"
#include "../satuzk/inline/util/Bool3.hpp"
#include "../satuzk/inline/util/SimpleVector.hpp"

#include "../satuzk/include/Antecedent.hpp"
#include "../satuzk/include/Conflict.hpp"
#include "../satuzk/include/Vars.hpp"

#include "../satuzk/inline/Antecedent.hpp"
#include "../satuzk/inline/Conflict.hpp"
#include "../satuzk/inline/Vars.hpp"
#include "../satuzk/inline/Clauses.hpp"

#include "../satuzk/inline/Propagate.hpp"
#include "../satuzk/inline/Lbd.hpp"
#include "../satuzk/inline/Learn.hpp"
#include "../satuzk/inline/Vsids.hpp"
#include "../satuzk/inline/ExtModel.hpp"
#include "../satuzk/include/Config.hpp"
#include "../satuzk/inline/Config.hpp"

#include "../satuzk/inline/simplify/BlockedClauseElim.hpp"
#include "../satuzk/inline/simplify/VarElim.hpp"
#include "../satuzk/inline/simplify/Subsumption.hpp"
#include "../satuzk/inline/simplify/Equivalent.hpp"
#include "../satuzk/inline/simplify/Unhiding.hpp"

namespace maxsatuzk {

class SolverUzk {
private:
	struct BaseDefs {
		typedef uint32_t LiteralIndex;
		typedef uint32_t ClauseIndex;
		typedef uint32_t ClauseLitIndex;
		typedef uint32_t Order;
		typedef double Activity;
		typedef uint32_t Declevel;
	};
	class Hooks {
	public:
		void onLearnedClause(satuzk::ClauseType<BaseDefs> clause) {
		}
	};

	typedef satuzk::Config<BaseDefs, Hooks> Config;
public:
	typedef Config::Literal Literal;
	typedef Config::Variable Variable;

	class HardAllocator {
	public:
		typedef Config::Literal Literal;
		typedef Config::Variable Variable;
		
		HardAllocator(SolverUzk &solver) : p_solver(solver) { }
		
		Variable allocate() {
			return p_solver.p_config.varAlloc();
		}

	private:
		SolverUzk &p_solver;
	};

	class HardEmitter {
	public:
		typedef Config::Literal Literal;
		typedef Config::Variable Variable;

		HardEmitter(SolverUzk &solver) : p_solver(solver) { }

		template<typename Iterator>
		void emit(Iterator begin, Iterator end) {
			p_solver.p_config.inputClause(end - begin, begin, end);
		}
	
	private:
		SolverUzk &p_solver;
	};
	
	class ContextEmitter {
	public:
		typedef Config::Literal Literal;
		typedef Config::Variable Variable;

		ContextEmitter(SolverUzk &solver) : p_solver(solver) {
			p_contextVar = p_solver.p_config.varAlloc();
			p_solver.p_config.lockVariable(p_contextVar);
			p_solver.p_config.assumptionEnable(p_contextVar.zeroLiteral());
		}

		~ContextEmitter() {
			p_solver.p_config.expellContaining(p_contextVar.oneLiteral());
			p_solver.p_config.expellContaining(p_contextVar.zeroLiteral());
		}

		template<typename Iterator>
		void emit(Iterator begin, Iterator end) {
			std::vector<Literal> lits(begin, end);
			lits.push_back(p_contextVar.oneLiteral());

			p_solver.p_config.inputClause(lits.size(), lits.begin(), lits.end());
		}
	
	private:
		SolverUzk &p_solver;
		Variable p_contextVar;
	};

public:
	enum class Result {
		kNone, kBreak, kError, kHardUnsat, kSoftUnsat, kSat
	};

	SolverUzk(int config_id) : p_config(Hooks(), config_id), p_numRuns(0) {
	}

	void lockVariable(Variable variable);
	void solveStart();
	Result solveStep();
	void solveReset();
	bool litInModel(Literal literal);
	
private:
	Config p_config;
	typename Config::Variable p_contextVar;
	int p_numRuns;
};

} // namespace maxsatuzk

