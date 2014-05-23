//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

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


using namespace std;
using namespace boost::program_options;

struct options
/*
  This object represents the command-line options
 */
{
  string mapfile,outfile;
  unsigned nperms_per_file;
  vector<string> infiles;
  bool convert;
  //constructor for empty object
  options(void);
};

options::options(void) : mapfile(string()),
			 outfile(string()),
			 nperms_per_file(0u),
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

  string line; //this is a line from the input file
  istringstream ibuffer; //We will turn each line into a buffer for reading
  ostringstream obuffer;
  double x;

  gzFile gzout = gzopen( O.outfile.c_str(),"w" );
  if (gzout == NULL )
    {
      cerr << "Error, could not open " 
	   << O.outfile 
	   << " for writing\n";
      exit(10);
    }

  for( unsigned infile = 0 ; infile < O.infiles.size() ; ++infile )
    {
      ifstream in( O.infiles[infile].c_str() );
      if( ! in )
	{
	  //Bad!
	  cerr << "Error, could not open " 
	       << O.infiles[infile]
	       << " for reading\n";
	  exit(10);
	}
      unsigned LC = 0; //line count
      while(! in.eof() )
	{
	  ++LC;
	  getline(in,line);
	  in >> ws; //ws is part of namespace std, and serves to "chomp" any extraenous whitespace

	  if( (LC==1 && infile==0) || (LC>1 && infile) )
	    /*
	      Only process the first line of a perm file (which is the observed values)
	      if we are dealing with the first file
	     */
	    {
	      istringstream ibuffer(line);
	      while(! ibuffer.eof() )
		{
		  ibuffer >> x >> ws;
		  if (O.convert) //input is chi-squared, and we will turn into p-value
		    {
		      double pv = gsl_cdf_chisq_Q(x,1.); //chi-sq w/one degree of freedom
		      obuffer.write( reinterpret_cast<char *>(&pv),sizeof(double) );
		    }
		  else //input is already a p-value, so leave it untouched
		    {
		      obuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
		    }
		  //Is buffer full, write it to file and clear it
		  if( obuffer.str().size() >= MBUFFER )
		    {
		      gzwrite( gzout, obuffer.str().c_str(), obuffer.str().size() );
		      obuffer.str( string() );//reset to an empty string
		    }
		}
	      /*
		Pro tip: you have to check for a non-empty buffer.
		This occurs when you have < MBUFFER reads near the end of the file
	      */
	      if( ! obuffer.str().empty() )
		{
		  gzwrite( gzout, obuffer.str().c_str(), obuffer.str().size() );
		  obuffer.str( string() );//reset to an empty string
		}
	    }
	}
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
    ("nperms,n",value<unsigned>(&rv.nperms_per_file)->default_value(0),"Number of permutations per input file")
    ("convert,c",value<bool>(&rv.convert)->default_value(true),"Convert chi-squared to p-value")
    ;

  parsed_options parsed = 
    command_line_parser(argc, argv).options(desc).allow_unregistered().run(); 

  if ( argc == 1 )
    {
      cerr << desc << '\n';
      exit(0);
    }
  else //check if --help/-h was input
    {
      for( unsigned i = 0 ; i < parsed.options.size() ; ++i )
	{
	  if( parsed.options[i].string_key == string("help") )
	    {
	      cerr << desc << '\n';
	      exit(0);
	    }
	}
    }
  rv.infiles = collect_unrecognized(parsed.options, include_positional);

  return rv;
}
