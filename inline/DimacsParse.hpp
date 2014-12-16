
#include <cassert>
#include <stdexcept>
#include <limits>
#include <vector>

#include <unistd.h>

namespace maxsatuzk {

class BaseParser {
public:
	BaseParser(int fd) : p_fd(fd),
			p_readPointer(0), p_readLimit(0), p_eof(false) {
		read();
	}

	//---------------------------//
	// low level functions	     //
	//---------------------------//

	void read() {
		assert(p_readPointer == p_readLimit);
		ssize_t len = ::read(p_fd, p_readBuffer, kReadBufferSize);
		if(len == -1)
			throw std::runtime_error("Could not read input");
		if(len == 0)
			p_eof = true;
		p_readLimit = len;
		p_readPointer = 0;
	}
	void consume() {
		p_readPointer++;
		if(p_readPointer == p_readLimit)
			read();
	}
	bool atEndOfFile() {
		return p_eof;
	}
	char fetch() {
		assert(p_readPointer < p_readLimit);
		return p_readBuffer[p_readPointer];
	}
	
	//---------------------------//
	// high level functions	     //
	//---------------------------//

	bool checkSpace() {
		return fetch() == ' ' || fetch() == '\t';
	}
	bool checkBreak() {
		return fetch() == '\n' || fetch() == '\r';
	}
	
	void skipSpace() {
		while(!atEndOfFile() && checkSpace())
			consume();
	}
	void skipSpaceOrBreak() {
		while(!atEndOfFile() && (checkSpace() || checkBreak()))
			consume();
	}
	void skipUntilBreakOrEnd() {
		while(!atEndOfFile() && !checkBreak())
			consume();
		if(!atEndOfFile())
			consume();
	}
	
	void forceBreakOrEnd() {
		if(atEndOfFile())
			return;
		if(!checkBreak())
			throw std::runtime_error("Expected end-of-line");
		consume();
	}

	void readWord() {
		p_tempLength = 0;
		while(!atEndOfFile() && !checkSpace() && !checkBreak()) {
			if(p_tempLength >= kTempBufferSize)
				throw std::logic_error("Buffer limit reached");
			p_tempBuffer[p_tempLength++] = fetch();
			consume();
		}
	}	
	bool bufferMatch(const char *string) {
		for(int i = 0; string[i] != 0; i++) {
			if(i >= p_tempLength)
				return false;
			if(string[i] != p_tempBuffer[i])
				return false;
		}
		return true;
	}
	int64_t bufferGetInt() {
		int64_t num = 0;
		bool neg = false;
		int i = 0;
		
		if(p_tempLength == 0)
			throw std::runtime_error("Expected number (empty buffer)");
		if(i < p_tempLength && p_tempBuffer[i] == '-') {
			neg = true;
			i++;
		}
		while(i < p_tempLength) {
			if(!(p_tempBuffer[i] >= '0' && p_tempBuffer[i] <= '9'))
				throw std::runtime_error("Expected number (illegal char)");
			num = num * 10 + (p_tempBuffer[i] - '0');
			i++;
		}
		return neg ? -num : num;
	}
	
private:
	static const int kReadBufferSize = 16 * 1024;
	static const int kTempBufferSize = 1024;
	
	int p_fd;
	
	char p_readBuffer[kReadBufferSize];
	int p_readPointer;
	int p_readLimit;
	bool p_eof;
	
	char p_tempBuffer[kTempBufferSize];
	int p_tempLength;
};

class CnfParser {
public:
	CnfParser(int fd) : base(fd) { }

	enum CnfType {
		kNone, kStandard, kWeighted
	};
	
	void parse(InClauseSpace &f) {
		CnfType cnf = CnfType::kNone;
		int n_vars = -1, n_clauses = -1;
		int64_t top = -1;

		while(!base.atEndOfFile()) {
			base.skipSpaceOrBreak();
			if(base.atEndOfFile())
				break;
			base.readWord();
			
			if(base.bufferMatch("c")) {
				base.skipUntilBreakOrEnd();
			}else if(base.bufferMatch("p")) {
				// read the file type
				base.skipSpace();
				base.readWord();
				if(base.bufferMatch("cnf")) {
					//std::cout << "It's a cnf file" << std::endl;
					cnf = CnfType::kStandard;
				}else if(base.bufferMatch("wcnf")) {
					//std::cout << "It's a wcnf file" << std::endl;
					cnf = CnfType::kWeighted;
				}else{
					throw std::runtime_error("Illegal cnf/wcnf file");
				}
				
				// read the number of vars and clauses
				base.skipSpace();
				base.readWord();
				n_vars = base.bufferGetInt();
				base.skipSpace();
				base.readWord();
				n_clauses = base.bufferGetInt();

				base.skipSpace();
				if(!base.checkBreak()) {
					base.readWord();
					top = base.bufferGetInt();
				}
				std::cout << "c Input vars: " << n_vars << ", clauses: " << n_clauses << ", top: " << top << std::endl;
				base.skipSpace();
				base.forceBreakOrEnd();
			}else{
				// the line contains a clause
				Weight weight = 1;
				if(cnf == CnfType::kWeighted) {
					weight = base.bufferGetInt();
					assert(weight < std::numeric_limits<int>::max());
					if(weight == 0)
						throw std::runtime_error("Illegal weight value");
					base.skipSpace();
					base.readWord();
				}
				if(top != -1 && weight == top)
					weight = kHardWeight;
				
				std::vector<int> lits;
				while(true) {
					int lit = base.bufferGetInt();
					if(lit == 0) {
						base.skipSpace();
						base.forceBreakOrEnd();
						break;
					}
					lits.push_back(lit);
					base.skipSpaceOrBreak();
					base.readWord();
				}
				if(lits.size() == 0)
					throw std::runtime_error("File contains empty clause");

				InClauseRef c = f.allocate(lits.size());
				for(int i = 0; i < lits.size(); i++) {
					InVariable var = InVariable::from_index(lits[i] < 0 ? -lits[i] : lits[i]);
					c.setLiteral(i, lits[i] < 0 ? var.neg() : var.pos());
				}
				c.setWeight(weight);
			}
		}
	}
private:
	BaseParser base;
};

} // namespace maxsatuzk

