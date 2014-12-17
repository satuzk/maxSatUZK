
#include <stdexcept>
#include <iostream>

#include <sys/fcntl.h>

#include "../include/Formula.hpp"
#include "../include/Sorter.hpp"
#include "../include/SolverUzk.hpp"
#include "../inline/Formula.hpp"
#include "../inline/Sorter.hpp"
#include "../inline/DimacsParse.hpp"

#include <encodeuzk/mixed-radix.inline.hpp>

int main(int argc, char **argv) {
	if(argc != 2)
		throw std::runtime_error("Expected arguments: <input>");

	maxsatuzk::InClauseSpace in;
	
	int fd = open(argv[1], O_RDONLY);
	if(fd == -1)
		throw std::runtime_error("Could not open input file");
	maxsatuzk::CnfParser parser(fd);
	parser.parse(in);
	std::cout << "c Parsing finished" << std::endl;

	std::vector<int> s;
	for(auto i = in.refsBegin(); i != in.refsEnd(); ++i) {
		maxsatuzk::InClauseRef clause = *i;
		if(clause.getWeight() != maxsatuzk::kHardWeight)
			s.push_back(clause.getWeight());
	}
	std::vector<int> base = encodeuzk::optimalBase(s);
	
	maxsatuzk::solve<maxsatuzk::SolverUzk>(in, base, 0, maxsatuzk::countSoftWeight(in) + 1);
	return 0;
}

