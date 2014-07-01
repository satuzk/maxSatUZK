
#include <limits>
#include <vector>

namespace maxsatuzk {

typedef int64_t Weight;
typedef uint32_t LiteralId;

static const Weight kHardWeight = -1;

class Variable;
class Literal;
class ClauseSpace;
class ClauseRef;
class ClauseLitIterator;

class Variable {
public:
	typedef int32_t index_type;
	typedef Literal literal;
	
	Variable() : p_index(0) { }
	static inline Variable null();
	static inline Variable from_index(int32_t index);
	inline bool operator ==(const Variable v) const;
	inline int32_t index() const;
	inline int num() const;
	inline literal pos() const;
	inline literal neg() const;

private:
	Variable(int32_t idx) : p_index(idx) { }
	
	int32_t p_index;
};

class Literal {
public:
	typedef int32_t index_type;
	typedef Variable variable;

	Literal() : p_index(0) { }
	static inline Literal null();
	static inline Literal from_index(int32_t index);
	inline bool operator ==(const Literal l) const;
	inline int32_t index() const;
	inline int sign() const;
	inline int num() const;
	inline variable var() const;
	inline Literal operator -();

private:
	Literal(int32_t idx) : p_index(idx) { }

	int32_t p_index;
};

class ClauseSpace {
friend class ClauseRef;
friend class ClauseLitIterator;
public:
	class RefIterator : public std::iterator<std::forward_iterator_tag, ClauseRef> {
	public:
		RefIterator(ClauseSpace &space, int index)
				: p_space(space), p_index(index) {
		}
		
		ClauseRef operator* ();
		void operator++ ();
		
		bool operator== (const RefIterator &other);
		bool operator!= (const RefIterator &other);
	private:
		ClauseSpace &p_space;
		int p_index;
	};

	ClauseRef allocate(int length);
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

class ClauseRef {
public:
	ClauseRef(ClauseSpace &space, int index)
		: p_space(space), p_index(index) { }

	ClauseSpace &space();
	int index();
	
	int length();
	Literal getLiteral(int index);
	void setLiteral(int index, Literal literal);

	ClauseLitIterator begin();
	ClauseLitIterator end();

	Weight getWeight();
	void setWeight(Weight weight);

private:
	ClauseSpace::HeadStruct &accessHead();
	
	ClauseSpace &p_space;
	int p_index;
};

class ClauseLitIterator : public std::iterator<std::forward_iterator_tag, LiteralId> {
public:
	ClauseLitIterator(ClauseSpace &space, int clause_id, int index)
			: p_space(space), p_clauseId(clause_id), p_index(index) { }

	LiteralId &operator* ();
	
	void operator++ ();
	bool operator== (const ClauseLitIterator &other);

private:	
	ClauseSpace::HeadStruct &accessHead();
	
	ClauseSpace &p_space;
	int p_clauseId;
	int p_index;
};

struct VarAllocator {
	typedef Variable variable;

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
	std::size_t operator() (const Variable v) const {
		return v.index();
	}
};

}; // namespace maxsatuzk

