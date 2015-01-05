CXX=c++
CXXFLAGS=-O2 -std=c++11  -L/data/apps/hdf5/1.8.11/lib -I/data/apps/hdf5/1.8.11/include -Wall -W -I.  
#CXXFLAGS=-pg -O2 -std=c++11  -L/data/apps/hdf5/1.8.11/lib -I/data/apps/hdf5/1.8.11/include -Wall -I.

TARGETS=perms2h5 esmk slabtest
all: perms2h5.o esmk.o H5util.o slabtest.o
	$(CXX) $(CXXFLAGS) -o perms2h5 perms2h5.o $(LDFLAGS) -lboost_program_options -lgsl -lgslcblas -lhdf5_cpp -lhdf5
	$(CXX) $(CXXFLAGS) -o esmk esmk.o H5util.o $(LDFLAGS) -lboost_program_options -lhdf5_cpp -lhdf5 -pthread
	$(CXX) $(CXXFLAGS) -o slabtest slabtest.o H5util.o $(LDFLAGS) -lboost_program_options -lhdf5_cpp -lhdf5 
clean:
	rm -f $(TARGETS) *.o

bin:
	mv $(TARGETS) ../bin
