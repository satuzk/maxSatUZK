
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
#include "../satuzk/inline/util/BinaryStack.hpp"
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

namespace maxsatuzk {

class SolverUzk : public Solver {
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
	SolverUzk(int config_id) : p_config(Hooks(), config_id) {
	}

	void reserveVars(int n);
	virtual void setupLhs(ClauseSpace &f);
	virtual void updateRhs(ClauseSpace &f);
	virtual void solveStart();
	virtual Result solveStep();
	virtual void solveReset();
	virtual bool litInModel(Literal literal);
	
private:
	Config p_config;
	typename Config::Variable p_contextVar;
};

} // namespace maxsatuzk

