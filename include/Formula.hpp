
#include <limits>
#include <vector>

namespace maxsatuzk {

typedef int64_t Weight;
typedef uint32_t LiteralId;

static const Weight kHardWeight = -1;

class InVariable;
class InLiteral;
class InClauseSpace;
class InClauseRef;
class InClauseLitIterator;

class InVariable {
public:
	typedef int32_t index_type;
	typedef InLiteral literal;
	
	InVariable() : p_index(0) { }
	static inline InVariable null();
	static inline InVariable from_index(int32_t index);
	inline bool operator ==(const InVariable v) const;
	inline int32_t index() const;
	inline int num() const;
	inline literal pos() const;
	inline literal neg() const;

private:
	InVariable(int32_t idx) : p_index(idx) { }
	
	int32_t p_index;
};

class InLiteral {
public:
	typedef int32_t index_type;
	typedef InVariable variable;

	InLiteral() : p_index(0) { }
	static inline InLiteral null();
	static inline InLiteral from_index(int32_t index);
	inline bool operator ==(const InLiteral l) const;
	inline int32_t index() const;
	inline int sign() const;
	inline int num() const;
	inline variable var() const;
	inline InLiteral operator -();

private:
	InLiteral(int32_t idx) : p_index(idx) { }

	int32_t p_index;
};

class InClauseSpace {
friend class InClauseRef;
friend class InClauseLitIterator;
public:
	class RefIterator : public std::iterator<std::forward_iterator_tag, InClauseRef> {
	public:
		RefIterator(InClauseSpace &space, int index)
				: p_space(space), p_index(index) {
		}
		
		InClauseRef operator* ();
		void operator++ ();
		
		bool operator== (const RefIterator &other);
		bool operator!= (const RefIterator &other);
	private:
		InClauseSpace &p_space;
		int p_index;
	};

	InClauseRef allocate(int length);
	RefIterator refsBegin();
	RefIterator refsEnd();

private:
	struct HeadStruct {
		int length;
		int pointer;
		Weight weight;
	};
	
	std::vector<HeadStruct> p_heads;
	std::vector<LiteralId> p_literals;
};

class InClauseRef {
public:
	InClauseRef(InClauseSpace &space, int index)
		: p_space(space), p_index(index) { }

	InClauseSpace &space();
	int index();
	
	int length();
	InLiteral getLiteral(int index);
	void setLiteral(int index, InLiteral literal);

	InClauseLitIterator begin();
	InClauseLitIterator end();

	Weight getWeight();
	void setWeight(Weight weight);

private:
	InClauseSpace::HeadStruct &accessHead();
	
	InClauseSpace &p_space;
	int p_index;
};

class InClauseLitIterator : public std::iterator<std::forward_iterator_tag, LiteralId> {
public:
	InClauseLitIterator(InClauseSpace &space, int clause_id, int index)
			: p_space(space), p_clauseId(clause_id), p_index(index) { }

	LiteralId &operator* ();
	
	void operator++ ();
	bool operator== (const InClauseLitIterator &other);

private:	
	InClauseSpace::HeadStruct &accessHead();
	
	InClauseSpace &p_space;
	int p_clauseId;
	int p_index;
};

struct VarAllocator {
	typedef InVariable variable;

	VarAllocator() : p_nextVar(1), p_varLimit(std::numeric_limits<int>::max()) { }
	VarAllocator(int first, int limit) : p_nextVar(first), p_varLimit(limit) { }

	uint32_t p_nextVar;
	uint32_t p_varLimit;

	int next_var() {
		return p_nextVar;
	}

	variable alloc() {
		if(p_nextVar >= p_varLimit)
			throw std::runtime_error("Reached variable limit " + std::to_string(p_varLimit));
		int index = p_nextVar++;
		return variable::from_index(index);
	}
};

class VarHashFunc {
public:
	std::size_t operator() (const InVariable v) const {
		return v.index();
	}
};

}; // namespace maxsatuzk

