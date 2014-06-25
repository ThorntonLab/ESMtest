CXX=c++
CXXFLAGS=-O2 -Wall -W -I.
#CXXFLAGS=-g -Wall

TARGETS=perms2h5 esmk slabunit

all: perms2h5.o esmk.o H5util.o slabunit.o
	$(CXX) $(CXXFLAGS) -o perms2h5 perms2h5.o $(LDFLAGS) -lboost_program_options -lgsl -lgslcblas -lhdf5_cpp -lhdf5
	$(CXX) $(CXXFLAGS) -o esmk esmk.o H5util.o $(LDFLAGS) -lboost_program_options -lhdf5_cpp -lhdf5
	$(CXX) $(CXXFLAGS) -o slabunit slabunit.o H5util.o $(LDFLAGS) -lboost_program_options -lhdf5_cpp -lhdf5
clean:
	rm -f $(TARGETS) *.o

bin:
	mv $(TARGETS) ../bin
