CXX=c++
CXXFLAGS=-O2 -Wall -W

all: mergeperms.o perms2gz.o mergepermsh5.o
	$(CXX) $(CXXFLAGS) -o perms2gz perms2gz.o $(LDFLAGS) -lboost_program_options -lgsl -lgslcblas -lz
	$(CXX) $(CXXFLAGS) -o mergeperms mergeperms.o $(LDFLAGS) -lboost_program_options  -lz
	$(CXX) $(CXXFLAGS) -o mergepermsh5 mergepermsh5.o $(LDFLAGS) -lhdf5_cpp -lhdf5 -lboost_program_options -lz