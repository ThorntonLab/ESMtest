//HPC modules needed: hdf5/1.8.11 boost/1.54.0 zlib/1.2.7

/*
  libhdf5 -- the C++ interface is in this header.

  This also exposes the C interface.
*/
#include <H5Cpp.h>

//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//other boost stuff
#include <boost/bind.hpp>

//zlib is included b/c it probably makes sense to write the output as a gzip file
#include <zlib.h>

//Standard c++
#include <iostream>   //let's us print to screen
#include <string>     //strings = containers of characters
#include <fstream>    //let's us write to files
#include <vector>     //vectors = containers of other things
#include <cmath>      //The C++ version of C's math.h (puts the C functions in namespace std)
#include <algorithm>  //find, sort, etc.
#include <set>        //A set is a container, see http://www.cplusplus.com/reference/set/set/

/*
  This is a header that I wrote.

  It contains overly-simple functions
  to read in 1-dimensional data from H5 files
*/
#include <H5util.hpp>

using namespace std;
using namespace boost::program_options;
using namespace H5;

//This is a data type to hold command-line options
struct esm_options
{
  string outfile;
  int winsize,jumpsize,K;
  vector<string> infiles;
};

struct within
/*
  Function object.   
  Returns true if val is within range
  specified by (left,right)
 */
{
  typedef bool result_type;
  inline bool operator()(const int & val,
			 const int & left,
			 const int & right) const
  {
    return (val >= left && val <= right);
  }
};

/*
  get_indexes is passed a list of positions on a chromosome and
  the left and right boundaries of a window.

  The return value is the pair of indexes within pos
  whose values are >= left and <= right.

  If no position exists satisfying one of these criteria,
  the maximum value of a size_t is returned.
 */
pair<size_t,size_t> get_indexes( const vector<int> & pos,
				 const int & left,
				 const int & right );
//Parse command line options
esm_options parseargs( int argc, char ** argv );
//Ask if all the permutation files contain the same marker info
bool permfilesOK( const esm_options & O );
//Runs the esm_k test on the data


void slab_unit( const esm_options & O );

int main( int argc, char ** argv )
{
  esm_options O = parseargs(argc,argv); 
  if( !permfilesOK(O) )
    {
      cerr << "Error with permutation files\n";
      exit(10);
    }
  slab_unit(O);

  exit(0);
}

esm_options parseargs( int argc, char ** argv )
{
  esm_options rv;

  options_description desc("Calculate ESM_K p-values in sliding window");
  desc.add_options()
    ("help,h", "Produce help message")
    ("outfile,o",value<string>(&rv.outfile),"Output file name.  Format is gzipped")
    ("winsize,w",value<int>(&rv.winsize),"Window size (bp)")
    ("jumpsize,j",value<int>(&rv.jumpsize),"Window jump size (bp)")
    ("K,k",value<int>(&rv.K),"Number of markers to use for ESM_k stat in a window.  Must be > 0.")
    ;

  variables_map vm;
  store( command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm );
  notify(vm);

  if( vm.count("help") || argc == 1 )
    {
      cerr << desc << '\n';
      exit(0);
    }
  
  parsed_options parsed = 
    command_line_parser(argc, argv).options(desc).allow_unregistered().run(); 
  
  if ( argc == 1 || vm.count("help") )
    {
      cerr << desc << '\n';
      exit(0);
    }

  if (!vm.count("outfile") || !vm.count("winsize") || !vm.count("jumpsize") || !vm.count("K") )
    {
      cerr << "Too few options given.\n"
	   << desc << '\n';
      exit(10);
    }
  rv.infiles = collect_unrecognized(parsed.options, include_positional);
  if(rv.infiles.empty())
    {
      cerr << "Error: no permutation files passed to program.\n"
	   << desc << '\n';
      exit(0);
    }
  return rv;
}

//make sure each perm file contains the same markers, etc.
bool permfilesOK( const esm_options & O )
{
  if( O.infiles.empty() ) { return false; }

  vector<string> chroms_0 = read_strings(O.infiles[0].c_str(),"/Markers/chr");
  set<string> sc_0(chroms_0.begin(),chroms_0.end());
  if( sc_0.size() > 1 ) { return false; }
  vector<string> markers_0 = read_strings(O.infiles[0].c_str(),"/Markers/IDs");
  vector<int> pos_0 = read_ints(O.infiles[0].c_str(),"/Markers/pos");  

  for ( size_t i = 1 ; i < O.infiles.size() ; ++i )
    {
      vector<string> chroms_i = read_strings(O.infiles[i].c_str(),"/Markers/chr");
      set<string> sc_i(chroms_i.begin(),chroms_i.end());
      if( sc_i.size() > 1 ) { return false; }
      vector<string> markers_i = read_strings(O.infiles[i].c_str(),"/Markers/IDs");
      vector<int> pos_i = read_ints(O.infiles[i].c_str(),"/Markers/pos");
      if( markers_0 != markers_i || sc_0 != sc_i || pos_0 != pos_i )
	{
	  return false;
	}  
    }

  return true;
}

pair<size_t,size_t> get_indexes( const vector<int> & pos,
				 const int & left,
				 const int & right )
{
  vector<int>::const_iterator ci1 = find_if(pos.begin(),pos.end(),
					    boost::bind( within(),_1,left,right));
  vector<int>::const_reverse_iterator ci2 = find_if(pos.rbegin(),pos.rend(),
						    boost::bind(within(),_1,left,right) );

  if( ci1 == pos.end() && ci2 == pos.rend() )
    {
      //no snps are in window
      return make_pair( numeric_limits<size_t>::max(),
			numeric_limits<size_t>::max() );
    }
  return make_pair( ci1-pos.begin(), pos.rend()-ci2-1 );
}
				 
void slab_unit( const esm_options & O )
{
  //Step 1: read in the marker data from the first file in 0.infiles:
  //1a: the chrom labels
  vector<string> chroms_0 = read_strings(O.infiles[0].c_str(),"/Markers/chr");
  set<string> sc_0(chroms_0.begin(),chroms_0.end());
  chroms_0.clear();
  //1b: the rsID for the markers
  vector<string> markers_0 = read_strings(O.infiles[0].c_str(),"/Markers/IDs");
  //1c: the marker positions
  vector<int> pos_0 = read_ints(O.infiles[0].c_str(),"/Markers/pos");  

  //Step 2: establish the left and right boundaries of the first window
  int left = 1, right = left + O.winsize - 1;

  //Step 3: go over the current data and make sure that our helper functions are working
  

  pair<size_t,size_t> indexes = get_indexes( pos_0, left,right );
  size_t markers = indexes.second - indexes.first + 1; 

  ofstream indexfile;
  indexfile.open("indexes");
  indexfile<<"left:"<<indexes.first<<"  right:"<<indexes.second<<"  markers:"<<markers<<"\n";
  indexfile.close();
  
  
  vector<double> slab = read_doubles_slab(O.infiles[0].c_str(),"/Perms/permutations",indexes.first,markers);
   
  ofstream slabfile;
  slabfile.open(O.outfile.c_str());

  
  for ( size_t i =0; i<slab.size(); ++i)
    {
      slabfile<< slab[i] << '\n' ;
    }
  slabfile.close();

  // const int LPOS = *(pos_0.end()-1); //This is the last position in pos_0.  Equivalent to pos[pos.size()-1], but I guess  I like to complicate things.
  /*do
    {
      //3a: get the indexes in pos_0 corresponding to left- and right- most SNPs in this window
      pair<size_t,size_t> indexes = get_indexes( pos_0, left,right );
      //3b: count number of markers in this window
      size_t nmarkers = indexes.second - indexes.first + 1;
      if( indexes.first != numeric_limits<size_t>::max() ) //If there are SNPs in the window
	{
	  vector<double> chisqs;
	  
	  chisqs.resize( nmarkers * O.infiles.size())
	    
	  for( size_t i = 0 ; i < O.infiles.size() : ++i )
	    {
	      vector<double>
	      
	    }
	  /*
	    Skip windows w/no markers.

	    This is where modifications need to be made. 

	    The permutation files all need to be processed here, and then the test statistic calculated for each window.

	    For the moment, forget about LD filtering.

	    I believe that the loop will look something like this:
	    vector<double> chisqs; //store permuted values of chisq in 1 big vector.

	    for( size_t i = 0 ; i < O.infiles.size() ; ++ i )
	    {
	    //in here, open the hdf5 file, figure out how many markers are in the 
	    //interval from (indexes.first,indexes.second), and read them into chisqs.

	    //Useful things to know about std::vector<T>:
	    //Declare an empty one:
	    //vector<double> x;
	    //resize to contain 1,000 doubles
	    //x.resize( 1000 );
	    //You'll need to pass the HDF5 library a pointer to where you want to read.
	    //If that position is size_t POS, 
	    //then passing &x[POS] passes an array of doubles (double * in C) to HDF5.
	    //This array will be a sub-array of x, beginning at position POS+1.
	    }
	  
	  cerr << left << ' ' << right << ' ' << pos_0[indexes.first] << ' ' << pos_0[indexes.second] << ' ' << nmarkers << '\n';
	}
       left += O.jumpsize;
      right += O.jumpsize;
    }
  while(left <= LPOS);*/  //Possibly not a great way to terminate.  Look @ libseq for guidance
}
