//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//other boost stuff
#include <boost/bind.hpp>

#include <H5Cpp.h>
//Headers for gzip output (C language)
//#include <zlib.h>

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
using namespace H5;

struct options
/*
  This object represents the command-line options
 */
{
  bool strip,convert;
  string mapfile,infile,outfile;
  options(void);
};

options::options(void) : strip(true),
			 convert(true),
			 infile(string()),
			 outfile(string())
{
}

options process_argv( int argc, char ** argv );
size_t process_mapfile( const options & O, H5File & ofile );
void process_perms( const options & O, size_t nmarkers, H5File & ofile );

int main( int argc, char ** argv )
{
  options O = process_argv( argc, argv );

  //Create output file
  H5File ofile( O.outfile.c_str() , H5F_ACC_TRUNC );

  size_t nmarkers = process_mapfile( O, ofile );
  process_perms( O, nmarkers, ofile );
  ofile.close();
  exit(0);
  //OK, write this marker info to an H5 file

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
		  cerr << buffer[val] << " -> ";
		  buffer[val] = gsl_cdf_chisq_Q(buffer[val],1.);
		  cerr << buffer[val] << '\n';
		}
	      val++;
	      if( val >= MBUFFER )
		{
		  val=0;
		}
	    }
	}
    }
  if(val)
    {
    }
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
    ("mapfile,m",value<string>(&rv.mapfile)->default_value(string()),".map file name.")
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

  bool bad_input = false;

  if (! vm.count("mapfile") )
    {
      cerr << "Error, no map file name specified\n";
      bad_input = true;
    }

  if ( bad_input )
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

void process_perms( const options & O, size_t nmarkers, H5File & ofile )
{
    FILE * ifp = !O.infile.empty() ? fopen( O.infile.c_str(),"r" ) : stdin;

    if ( ifp == NULL )
      {
	cerr << "Error, input stream could not be opened.\n";
	exit(10);
      }

    vector< double > data(nmarkers);
    int repno;
    fscanf(ifp,"%d",&repno);
    //The first line is the observed data
    for( size_t i = 0 ; i < nmarkers ; ++i )
      {
	int rv = fscanf(ifp,"%lf",&data[i]);
	if(O.convert)
	  {
	    data[i] = gsl_cdf_chisq_Q(data[i],1.);
	  }
      }

    ofile.createGroup("/Perms");

    DSetCreatPropList cparms;
    hsize_t chunk_dims[1] = {nmarkers};
    hsize_t maxdims[1] = {nmarkers};
    
    cparms.setChunk( 1, chunk_dims );
    cparms.setDeflate( 6 ); //compression level makes a big differences in large files!  Default is 0 = uncompressed.
  
    DataSpace * dataspace = new DataSpace(1,chunk_dims, maxdims);

    DataSet * d = new DataSet(ofile.createDataSet("/Perms/observed",
						  PredType::NATIVE_DOUBLE,
						  *dataspace,
						  cparms));

    d->write( data.data(), PredType::NATIVE_DOUBLE );

    delete dataspace;
    delete d;

    //ok, now we write a big matrix of the permuted values
    hsize_t chunk_dims2[2] = {1, nmarkers};
    hsize_t maxdims2[2] = {H5S_UNLIMITED,nmarkers};
    hsize_t datadims[2] = {0,nmarkers};
    hsize_t offsetdims[2] = {0,0};
    hsize_t recorddims[2] = {1,nmarkers};
    cparms.setChunk( 2, chunk_dims2 );
    cparms.setDeflate( 6 );

    dataspace = new DataSpace(2, datadims, maxdims2);
    d = new DataSet(ofile.createDataSet("/Perms/permutations",
					PredType::NATIVE_DOUBLE,
					*dataspace,
					cparms));

    DataSpace memspace(2,recorddims);
    while (fscanf(ifp,"%d",&repno) != -1 )
      {
	//Read in the data...
	for( size_t i = 0 ; i < nmarkers ; ++i )
	  {
	    int rv = fscanf(ifp,"%lf",&data[i]);
	    if( rv == 0 || rv == -1 )
	      {
		cerr << "Error, input stream ended before expected...\n";
		exit(10);
	      }
	    if(O.convert)
	      {
		data[i] = gsl_cdf_chisq_Q(data[i],1.);
	      }
	  }
	++datadims[0];
	offsetdims[0]=datadims[0]-1;
	d->extend( datadims );
	*dataspace = d->getSpace();
	dataspace->selectHyperslab(H5S_SELECT_SET, recorddims , offsetdims);
	d->write( data.data(), PredType::NATIVE_DOUBLE, memspace,*dataspace );
      }
}
