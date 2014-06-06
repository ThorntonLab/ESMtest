//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//Used to check ZLIB version number during compile time
#include <boost/static_assert.hpp>

//Headers for gzip output (C language)
#include <zlib.h>

//hdf5 library
#include <H5Cpp.h>

//Headers to conver chi-squared statistic into chi-squared p-value.  GNU Scientific Library (C language)
#include <gsl/gsl_cdf.h>

//standard C++ headers that we need
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>

/*
  Require zlib 1.2.5 or greater.  Compilation fails if this is not true.
  The ZLIB_VERNUM is defined in zlib.h
*/
BOOST_STATIC_ASSERT(ZLIB_VERNUM >= 0x1250);

using namespace std;
using namespace boost::program_options;
using namespace H5;

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
size_t process_mapfile( const options & O, H5File & ofile );


int main( int argc, char ** argv )
{
  options O = process_argv( argc, argv );

  //Open our hdf5 output file
  H5File ofile( O.outfile.c_str() , H5F_ACC_TRUNC );
  size_t nmarkers_mfile = process_mapfile( O, ofile );

  //Now, go through the perms
  unsigned long permno=1;     
  vector<double> perms(nmarkers_mfile);

  
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

      if ( nmarkers != nmarkers_mfile )
	{
	  cerr << "Error, permutation file " << O.infiles[infile]
	       << " contains " << nmarkers << " markers, but "
	       << O.mapfile << " contains " << nmarkers_mfile
	       << " markers\n";
	  exit(10);
	}
      unsigned nrecs=0;
      int rv;
      do
	{
	  rv = gzread( gzin, &perms[0], nmarkers*sizeof(double) );
	  if ( rv != -1 && rv != 0)
	    {
	      ostringstream gname;
	      if( nrecs == 0 )//these are the observed data
		{
		  if( infile == 0 ) //let's write them for the first file
		    {
		    }
		}
	      else
		{
		}
	      ++nrecs;
	    }
	}
      while(! gzeof( gzin ) );
      gzclose(gzin);
    }
  ofile.close();
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

size_t process_mapfile( const options & O, H5File & ofile )
/*
  Write map data into an h5 group called "Markers"
 */
{
  ifstream mapin( O.mapfile.c_str() );
  if (! mapin )
    {
      cerr << "Error, " << O.mapfile
	   << " could not be opened for reading\n";
      exit(10);
    }

  vector<string> markers,chroms;
  vector< const char * > marker_str,chrom_str;
  vector< int > vpos;
  string chrom,marker;
  int pos;
  while( !mapin.eof() )
    {
      mapin >> chrom >> marker >> pos >> ws;
      chroms.push_back( chrom );
      markers.push_back( marker );
      vpos.push_back( pos );
    }

  for( unsigned i = 0 ; i < markers.size() ; ++i )
    {
      marker_str.push_back( markers[i].c_str() );
      chrom_str.push_back( chroms[i].c_str() );
    }

  ofile.createGroup("/Markers");

  DSetCreatPropList cparms;
  hsize_t chunk_dims[1] = {markers.size()};
  hsize_t maxdims[1] = {markers.size()};

  cparms.setChunk( 1, chunk_dims );
  cparms.setDeflate( 6 ); //compression level makes a big differences in large files!  Default is 0 = uncompressed.
  
  DataSpace dataspace(1,chunk_dims, maxdims);

  H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE); 

  DataSet marker_dset = ofile.createDataSet("/Markers/IDs",
					    datatype,
					    dataspace,
					    cparms);

  marker_dset.write(marker_str.data(), datatype );

  DataSet chrom_dset = ofile.createDataSet("/Markers/chr",
					    datatype,
					    dataspace,
					    cparms);

  chrom_dset.write(chrom_str.data(), datatype );

  DataSet pos_dset = ofile.createDataSet("/Markers/pos",
					 PredType::NATIVE_INT,
					 dataspace,
					 cparms);

  pos_dset.write( vpos.data(),
		  PredType::NATIVE_INT );

  return markers.size();
}
