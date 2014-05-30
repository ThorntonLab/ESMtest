//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//other boost stuff
#include <boost/bind.hpp>

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
#include <algorithm>
#include <cstdlib>

//C headers
//#include <unistd.h>

using namespace std;
using namespace boost::program_options;

struct options
/*
  This object represents the command-line options
 */
{
  bool strip,convert;
  string infile,outfile;
  options(void);
};

options::options(void) : strip(true),
			 convert(true),
			 infile(string()),
			 outfile(string())
{
}

//process command lines
options process_argv( int argc, char ** argv );

int main( int argc, char ** argv )
{
  options O = process_argv( argc, argv );

  //Can we create the output file?
  gzFile ogz = gzopen(O.outfile.c_str(),"w");
  if( ogz == NULL )
    {
      cerr << "Error, could not open " << O.outfile
	   << " for writing\n";
      exit(10);
    }
  /*
    We will use C-style input, as it is faster.

    The first thing we need to figure out is how many markers were
    permute.

    The first line that comes out of PLINK 1.90a perms 
    contains the observed data
  */

  FILE * ifp = !O.infile.empty() ? fopen( O.infile.c_str(),"r" ) : stdin;

  if ( ifp == NULL )
    {
      cerr << "Error, could not open stream for reading\n";
      exit(10);
    }
  string line;
  int c;
  unsigned nspaces = 0;
  do
    {
      c = fgetc(ifp);
      if( ::isspace(c) && c != '\n') 
	{
	  ++nspaces;
	}
      line += char(c);
    }
  while( c != int('\n') ); //read in 1 char at a time until newline

  gzwrite( ogz, reinterpret_cast< char * >(&nspaces), sizeof(unsigned) );

  unsigned rep;
  double x;
  if( !O.strip ) //if we do not skip the observed data
    {
      istringstream temp(line);
      while(!temp.eof())
	{
	  temp >> x >> ws;
	  if( O.convert )
	    {
	      x = gsl_cdf_chisq_Q(x,1);
	    }
	  gzwrite( ogz, reinterpret_cast<char *>(&x),sizeof(double) );
	}
    }

  //Max buffer size will be 10MB uncompressed
  unsigned MBUFFER = 5*1024*1024/sizeof(double);
  double * buffer = new double[MBUFFER];
  size_t val = 0;
  int rv;
  while(! feof( ifp ) )
    {
      rv = fscanf(ifp,"%u",&rep);
      if(rv != -1)
	{
	  for( unsigned i = 0 ; i < nspaces ; ++i )
	    {
	      rv = fscanf(ifp,"%lf",&buffer[val]);
	      if( O.convert )
		{
		  buffer[val] = gsl_cdf_chisq_Q(buffer[val],1.);
		}
	      val++;
	      if( val >= MBUFFER )
		{
		  gzwrite( ogz,reinterpret_cast<char *>(buffer),val*sizeof(double) );
		  val=0;
		}
	    }
	}
    }
  if(val)
    {
      gzwrite( ogz,reinterpret_cast<char *>(buffer),val*sizeof(double) );
    }
  gzflush( ogz ,Z_FINISH );
  gzclose( ogz );

  exit(0);
}

options process_argv( int argc, char ** argv )
{
  options rv;

  options_description desc("Reads PLINK permuted data table from IFP.  Reformats to a binary-format data stream that is then written to a gzipped file.");
  desc.add_options()
    ("help,h", "Produce help message")
    ("nostrip","Do not strip the observed data from the file (default is to strip)")
    ("noconvert","Do not convert input into a p-value.  Default is to assume that the input is a chi^2 statistic with 1 degree of freedom")
    ("infile,i",value<string>(&rv.infile)->default_value(string()),"Input file name.  Default is to read from stdin")
    ("outfile,o",value<string>(&rv.outfile)->default_value(string()),"Output file name.  Format is binary and gzipped")
    ;

  variables_map vm;
  store( command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm );
  notify(vm);
  
  if ( argc == 1 || vm.count("help") )
    {
      cerr << desc << '\n';
      exit(0);
    }

  if( vm.count("nostrip") )
    {
      rv.strip = false;
    }

  if( vm.count("noconvert") )
    {
      rv.convert = false;
    }

  return rv;
}
