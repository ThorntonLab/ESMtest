//libhdf5
#include <H5Cpp.h>

//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//other boost stuff
#include <boost/bind.hpp>

//zlib
#include <zlib.h>

//Standard c++
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>

//Headers for this project
#include <H5util.hpp>

using namespace std;
using namespace boost::program_options;
using namespace H5;

struct esm_options
{
  string outfile;
  int winsize,jumpsize,K;
  vector<string> infiles;
};

esm_options parseargs( int argc, char ** argv );
bool permfilesOK( const esm_options & O );
void run_test( const esm_options & O );

int main( int argc, char ** argv )
{
  esm_options O = parseargs(argc,argv); 
  if( !permfilesOK(O) )
    {
      cerr << "Error with permutation files\n";
      exit(10);
    }
  run_test(O);

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
      vector<string> chroms_i = read_strings(O.infiles[0].c_str(),"/Markers/chr");
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

struct within
{
  typedef bool result_type;
  inline bool operator()(const int & val,
			 const int & left,
			 const int & right) const
  {
    return (val >= left && val <= right);
  }
};
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
				 
void run_test( const esm_options & O )
{
  vector<string> chroms_0 = read_strings(O.infiles[0].c_str(),"/Markers/chr");
  set<string> sc_0(chroms_0.begin(),chroms_0.end());
  chroms_0.clear();
  vector<string> markers_0 = read_strings(O.infiles[0].c_str(),"/Markers/IDs");
  vector<int> pos_0 = read_ints(O.infiles[0].c_str(),"/Markers/pos");  

  int left = 1, right = left + O.winsize - 1;

  unsigned j=0;
  const int LPOS = *(pos_0.end()-1);
  do
    {
      pair<size_t,size_t> indexes = get_indexes( pos_0, left,right );
      size_t nmarkers = indexes.second - indexes.first + 1;
      left += O.jumpsize;
      right += O.jumpsize;
      if( indexes.first != numeric_limits<size_t>::max() )
	{
	  cerr << left << ' ' << right << ' ' << pos_0[indexes.first] << ' ' << pos_0[indexes.second] << ' ' << nmarkers << '\n';
	}
      //Iterate over perm files
    }
  while(left <= LPOS);  //not a great way to terminate.  Look @ libseq for guidance
}
