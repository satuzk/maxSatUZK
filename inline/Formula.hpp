
#include <limits>

namespace maxsatuzk {

//---------------------------------------------------------------------------//
// Variable functions                                                        //
//---------------------------------------------------------------------------//

Variable Variable::null() {
	return Variable(0);
}
Variable Variable::from_index(int32_t index) {
	return Variable(index);
}	
bool Variable::operator ==(const Variable v) const {
	return p_index == v.p_index;
}
int32_t Variable::index() const {
	return p_index;
}
int Variable::num() const {
	return p_index;
}
Literal Variable::pos() const {
	return Literal::from_index(p_index * 2);
}
Literal Variable::neg() const {
	return Literal::from_index(p_index * 2 + 1);
}

//---------------------------------------------------------------------------//
// Literal functions                                                         //
//---------------------------------------------------------------------------//

Literal Literal::null() {
	return Literal(0);
}
Literal Literal::from_index(int32_t index) {
	return Literal(index);	
}
bool Literal::operator ==(const Literal l) const {
	return p_index == l.p_index;
}
int32_t Literal::index() const {
	return p_index;
}
int Literal::sign() const {
	if(p_index % 2 == 0)
		return 1;
	return -1;
}
int Literal::num() const {
	return sign() * var().num();
}
Variable Literal::var() const {
	return Variable::from_index(p_index / 2);
}
Literal Literal::operator-() {
	if(sign() > 0)
		return var().neg();
	return var().pos();
}

//---------------------------------------------------------------------------//
// definition of ClauseSpace functions                                       //
//---------------------------------------------------------------------------//

inline ClauseRef ClauseSpace::allocate(int length) {
	int index = p_heads.size();

	HeadStruct new_head;
	new_head.length = length;
	new_head.pointer = p_literals.size();

	p_heads.push_back(new_head);
	for(int i = 0; i < length; i++)
		p_literals.push_back(0);

	return ClauseRef(*this, index);
}

inline ClauseSpace::RefIterator ClauseSpace::refsBegin() {
	return RefIterator(*this, 0);
}
inline ClauseSpace::RefIterator ClauseSpace::refsEnd() {
	return RefIterator(*this, p_heads.size());
}

inline ClauseRef ClauseSpace::RefIterator::operator* () {
	return ClauseRef(p_space, p_index);
}
inline void ClauseSpace::RefIterator::operator++ () {
	p_index++;
}
inline bool ClauseSpace::RefIterator::operator== (const RefIterator &other) {
	return p_index == other.p_index;
}
inline bool ClauseSpace::RefIterator::operator!= (const RefIterator &other) {
	return !(*this == other);
}

//---------------------------------------------------------------------------//
// definition of ClauseRef functions                                         //
//---------------------------------------------------------------------------//
inline ClauseSpace &ClauseRef::space() {
	return p_space;
}
inline int ClauseRef::index() {
	return p_index;
}

inline int ClauseRef::length() {
	return accessHead().length;
}
inline Literal ClauseRef::getLiteral(int index) {
	return Literal::from_index(*std::next(begin(), index));
}
inline void ClauseRef::setLiteral(int index, Literal literal) {
	*std::next(begin(), index) = literal.index();
}

inline ClauseLitIterator ClauseRef::begin() {
	return ClauseLitIterator(p_space, p_index, 0);
}
inline ClauseLitIterator ClauseRef::end() {
	return ClauseLitIterator(p_space, p_index, accessHead().length);
}

inline Weight ClauseRef::getWeight() {
	return accessHead().weight;
}
inline void ClauseRef::setWeight(Weight weight) {
	accessHead().weight = weight;
}

inline ClauseSpace::HeadStruct &ClauseRef::accessHead() {
	return p_space.p_heads[p_index];
}

//---------------------------------------------------------------------------//
// definition of ClauseLitIterator functions                                 //
//---------------------------------------------------------------------------//

inline LiteralId &ClauseLitIterator::operator* () {
	return p_space.p_literals[accessHead().pointer + p_index];
}
	
inline void ClauseLitIterator::operator++ () {
	p_index++;
}
inline bool ClauseLitIterator::operator== (const ClauseLitIterator &other) {
	return p_index == other.p_index;
}

inline ClauseSpace::HeadStruct &ClauseLitIterator::accessHead() {
	return p_space.p_heads[p_clauseId];
}

}; // namespace maxsatuzk

