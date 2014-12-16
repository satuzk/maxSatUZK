
#include <limits>

namespace maxsatuzk {

//---------------------------------------------------------------------------//
// InVariable functions                                                        //
//---------------------------------------------------------------------------//

InVariable InVariable::null() {
	return InVariable(0);
}
InVariable InVariable::from_index(int32_t index) {
	return InVariable(index);
}	
bool InVariable::operator ==(const InVariable v) const {
	return p_index == v.p_index;
}
int32_t InVariable::index() const {
	return p_index;
}
int InVariable::num() const {
	return p_index;
}
InLiteral InVariable::pos() const {
	return InLiteral::from_index(p_index * 2);
}
InLiteral InVariable::neg() const {
	return InLiteral::from_index(p_index * 2 + 1);
}

//---------------------------------------------------------------------------//
// InLiteral functions                                                         //
//---------------------------------------------------------------------------//

InLiteral InLiteral::null() {
	return InLiteral(0);
}
InLiteral InLiteral::from_index(int32_t index) {
	return InLiteral(index);	
}
bool InLiteral::operator ==(const InLiteral l) const {
	return p_index == l.p_index;
}
int32_t InLiteral::index() const {
	return p_index;
}
int InLiteral::sign() const {
	if(p_index % 2 == 0)
		return 1;
	return -1;
}
int InLiteral::num() const {
	return sign() * var().num();
}
InVariable InLiteral::var() const {
	return InVariable::from_index(p_index / 2);
}
InLiteral InLiteral::operator-() {
	if(sign() > 0)
		return var().neg();
	return var().pos();
}

//---------------------------------------------------------------------------//
// definition of InClauseSpace functions                                       //
//---------------------------------------------------------------------------//

inline InClauseRef InClauseSpace::allocate(int length) {
	int index = p_heads.size();

	HeadStruct new_head;
	new_head.length = length;
	new_head.pointer = p_literals.size();

	p_heads.push_back(new_head);
	for(int i = 0; i < length; i++)
		p_literals.push_back(0);

	return InClauseRef(*this, index);
}

inline InClauseSpace::RefIterator InClauseSpace::refsBegin() {
	return RefIterator(*this, 0);
}
inline InClauseSpace::RefIterator InClauseSpace::refsEnd() {
	return RefIterator(*this, p_heads.size());
}

inline InClauseRef InClauseSpace::RefIterator::operator* () {
	return InClauseRef(p_space, p_index);
}
inline void InClauseSpace::RefIterator::operator++ () {
	p_index++;
}
inline bool InClauseSpace::RefIterator::operator== (const RefIterator &other) {
	return p_index == other.p_index;
}
inline bool InClauseSpace::RefIterator::operator!= (const RefIterator &other) {
	return !(*this == other);
}

//---------------------------------------------------------------------------//
// definition of InClauseRef functions                                         //
//---------------------------------------------------------------------------//
inline InClauseSpace &InClauseRef::space() {
	return p_space;
}
inline int InClauseRef::index() {
	return p_index;
}

inline int InClauseRef::length() {
	return accessHead().length;
}
inline InLiteral InClauseRef::getLiteral(int index) {
	return InLiteral::from_index(*std::next(begin(), index));
}
inline void InClauseRef::setLiteral(int index, InLiteral literal) {
	*std::next(begin(), index) = literal.index();
}

inline InClauseLitIterator InClauseRef::begin() {
	return InClauseLitIterator(p_space, p_index, 0);
}
inline InClauseLitIterator InClauseRef::end() {
	return InClauseLitIterator(p_space, p_index, accessHead().length);
}

inline Weight InClauseRef::getWeight() {
	return accessHead().weight;
}
inline void InClauseRef::setWeight(Weight weight) {
	accessHead().weight = weight;
}

inline InClauseSpace::HeadStruct &InClauseRef::accessHead() {
	return p_space.p_heads[p_index];
}

//---------------------------------------------------------------------------//
// definition of InClauseLitIterator functions                                 //
//---------------------------------------------------------------------------//

inline LiteralId &InClauseLitIterator::operator* () {
	return p_space.p_literals[accessHead().pointer + p_index];
}
	
inline void InClauseLitIterator::operator++ () {
	p_index++;
}
inline bool InClauseLitIterator::operator== (const InClauseLitIterator &other) {
	return p_index == other.p_index;
}

inline InClauseSpace::HeadStruct &InClauseLitIterator::accessHead() {
	return p_space.p_heads[p_clauseId];
}

}; // namespace maxsatuzk

