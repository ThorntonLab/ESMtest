//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//Used to check ZLIB version number during compile time
#include <boost/static_assert.hpp>

//Headers for gzip output (C language)
#include <zlib.h>

//Headers to conver chi-squared statistic into chi-squared p-value.  GNU Scientific Library (C language)
#include <gsl/gsl_cdf.h>

//standard C++ headers that we need
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

/*
  Require zlib 1.2.7 or greater.  Compilation fails if this is not true.
  The ZLIB_VERNUM is defined in zlib.h
*/
BOOST_STATIC_ASSERT(ZLIB_VERNUM >= 0x1250);

using namespace std;
using namespace boost::program_options;

struct options
/*
  This object represents the command-line options
 */
{
  string mapfile,outfile;
  vector<string> infiles;
  bool convert;
  //constructor for empty object
  options(void);
};

options::options(void) : mapfile(string()),
			 outfile(string()),
			 infiles( vector<string>() ),
			 convert(true)
{
}

//process command lines
options process_argv( int argc, char ** argv );

//Max buffer size will be 10MB uncompressed
const unsigned MBUFFER = 10*1024*1024;

int main( int argc, char ** argv )
{
  options O = process_argv( argc, argv );

 //Process the map file
  ifstream in(O.mapfile.c_str());
  if ( ! in )
    {
      cerr << "Error, " << O.mapfile
	   << " could not be read as a plain-text file for reading\n";
      exit(10);
    }
  string chrom,id;
  int pos,len;
  ostringstream mapbuffer;
  unsigned maprecords=0;
  while(! in.eof() )
    {
      in >> chrom >> id >> pos >> ws;
      ++maprecords;
      len = id.size();
      mapbuffer.write( reinterpret_cast<char*>(&len),sizeof(unsigned) );
      mapbuffer.write( id.c_str(), len*sizeof(char) );
      id = chrom.size();
      mapbuffer.write( reinterpret_cast<char*>(&len),sizeof(unsigned) );
      mapbuffer.write( chrom.c_str(), len*sizeof(char) );
      mapbuffer.write( reinterpret_cast<char*>(&pos),sizeof(int) );
    }

  gzFile gzout = gzopen( O.outfile.c_str(),"w" );
  if (gzout == NULL )
    {
      cerr << "Error, could not open " 
	   << O.outfile 
	   << " for writing\n";
      exit(10);
    }

  //write the map file info
  gzwrite( gzout, reinterpret_cast<char *>(&maprecords), sizeof(unsigned) );

  gzwrite( gzout, mapbuffer.str().c_str(), mapbuffer.str().size() );
 
  for( unsigned infile = 0 ; infile < O.infiles.size() ; ++infile )
    {
      gzFile gzin = gzopen( O.infiles[infile].c_str(),"r" );
      /*
	Increase buffer size from default of 8kb.
	This supposedly increased read speeds.
      */
      gzbuffer( gzin, 128*1024 ); 
      if (gzin == NULL )
	{
	  cerr << "Error, could not open "
	       << O.infiles[infile] 
	       << " for reading\n";
	  exit(10);
	}
      unsigned nmarkers;
      gzread( gzin, &nmarkers, sizeof(unsigned) );

      if ( nmarkers != maprecords )
	{
	  cerr << "Error, permutation file " << O.infiles[infile]
	       << " contains " << nmarkers << " markers, but "
	       << O.mapfile << " contains " << maprecords 
	       << " markers\n";
	  exit(10);
	}
      vector<double> perms(nmarkers);
      unsigned nrecs=0;
      int rv;
      do
	{
	  rv = gzread( gzin, &perms[0], nmarkers*sizeof(double) );
	  //obuff.write( reinterpret_cast< char * >(&perms[0]), nmarkers*sizeof(double) );
	  if ( rv != -1 && rv != 0)
	    {
	      gzwrite( gzout,
		       reinterpret_cast< char * >(&perms[0]), nmarkers*sizeof(double) );
	    }
	}
      while(! gzeof( gzin ) );
      gzclose(gzin);
    }
  gzclose(gzout);
}

options process_argv( int argc, char ** argv )
{
  options rv;

  options_description desc("Merge permuted statistics into table for ESM test");
  desc.add_options()
    ("help,h", "Produce help message")
    ("mapfile,m",value<string>(&rv.mapfile)->default_value(string()),"The .map file from PLINK")
    ("outfile,o",value<string>(&rv.outfile)->default_value(string()),"Output file name.  Format is binary and gzipped")
    ("convert,c",value<bool>(&rv.convert)->default_value(true),"Convert chi-squared to p-value")
    ;

  variables_map vm;
  store( command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm );
  notify(vm);

  
  parsed_options parsed = 
    command_line_parser(argc, argv).options(desc).allow_unregistered().run(); 
  
  if ( argc == 1 || vm.count("help") )
    {
      cerr << desc << '\n';
      exit(0);
    }
  rv.infiles = collect_unrecognized(parsed.options, include_positional);
  return rv;
}
