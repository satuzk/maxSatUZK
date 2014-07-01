
namespace maxsatuzk {

class Solver {
public:
	enum class Result {
		kNone, kBreak, kError, kHardUnsat, kSoftUnsat, kSat
	};

	virtual ~Solver() { };

	virtual void setupLhs(ClauseSpace &f) = 0;
	virtual void updateRhs(ClauseSpace &f) = 0;
	virtual void solveStart() = 0;
	virtual Result solveStep() = 0;
	virtual void solveReset() = 0;
	virtual bool litInModel(Literal literal) = 0;
};

} // namespace maxsatuzk

