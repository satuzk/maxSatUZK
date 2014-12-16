
CPP= src/Main.cpp src/Sorter.cpp src/OptimizeBase.cpp src/SolverUzk.cpp satuzk/src/sys/Linux.cpp
LIB= -lrt

testing:
	g++ -o maxSatUZK -O3 -std=c++0x -I ../encodeUZK/include $(CPP) $(LIB)
debug:
	g++ -o maxSatUZK -g -std=c++0x -I ../encodeUZK/include $(CPP) $(LIB)

