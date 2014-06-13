CXX=c++
CXXFLAGS=-O2 -Wall -W

TARGETS=perms2h5

all: perms2h5.o
	$(CXX) $(CXXFLAGS) -o perms2h5 perms2h5.o $(LDFLAGS) -lboost_program_options -lgsl -lgslcblas -lhdf5_cpp -lhdf5

clean:
	rm -f $(TARGETS) *.o

bin:
	mv $(TARGETS) ../bin
