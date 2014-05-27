CXX=c++
CXXFLAGS=-O2 -Wall -W

all: mergeperms.o perms2gz.o
	$(CXX) $(CXXFLAGS) -o perms2gz perms2gz.o $(LDFLAGS) -lboost_program_options -lgsl -lgslcblas -lz
	$(CXX) $(CXXFLAGS) -o mergeperms mergeperms.o $(LDFLAGS) -lboost_program_options -lgsl -lgslcblas -lz