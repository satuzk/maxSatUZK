
#include <stdexcept>
#include <iostream>

#include <sys/fcntl.h>

#include "../include/Formula.hpp"
#include "../include/Sorter.hpp"
#include "../inline/Formula.hpp"
#include "../inline/DimacsParse.hpp"

namespace maxsatuzk {
	NumberSeq optimalBase(const std::vector<int> &s);
}

int main(int argc, char **argv) {
	if(argc != 2)
		throw std::runtime_error("Expected arguments: <input>");

	maxsatuzk::ClauseSpace in;
	
	int fd = open(argv[1], O_RDONLY);
	if(fd == -1)
		throw std::runtime_error("Could not open input file");
	maxsatuzk::CnfParser parser(fd);
	parser.parse(in);
	std::cout << "c Parsing finished" << std::endl;

	std::vector<int> s;
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		maxsatuzk::ClauseRef clause = *i;
		if(clause.getWeight() != maxsatuzk::kHardWeight)
			s.push_back(clause.getWeight());
	}
	maxsatuzk::NumberSeq base = maxsatuzk::optimalBase(s);
	
	maxsatuzk::solve(in, base, 0, maxsatuzk::countSoftWeight(in) + 1);
	return 0;
}

