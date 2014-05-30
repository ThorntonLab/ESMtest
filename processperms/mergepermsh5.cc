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
  Require zlib 1.2.7 or greater.  Compilation fails if this is not true.
  The ZLIB_VERNUM is defined in zlib.h
*/
//BOOST_STATIC_ASSERT(ZLIB_VERNUM >= 0x1270);

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

//Max buffer size will be 10MB uncompressed
const unsigned MBUFFER = 10*1024*1024;

int main( int argc, char ** argv )
{
  options O = process_argv( argc, argv );

  //Open our hdf5 output file
  H5File *file = new H5File( O.outfile.c_str() , H5F_ACC_TRUNC );

  //Create group for marker info
  Group * group = new Group(file->createGroup("/Markers"));

  //Process the map file
  ifstream in(O.mapfile.c_str());
  if ( ! in )
    {
      cerr << "Error, " << O.mapfile
	   << " could not be read as a plain-text file for reading\n";
      exit(10);
    }
  string chrom,id;
  int pos,longest_chrom=0,longest_id=0;
  ostringstream mapbuffer;
  int maprecords=0;//this will also be the dim of our data space
  vector<unsigned> indexes,poss;
  vector<string> ids,chroms;
  while(! in.eof() )
    {
      in >> chrom >> id >> pos >> ws;
      longest_chrom= max(longest_chrom,int(chrom.length()));
      longest_id = max(longest_id,int(id.length()));
      ids.push_back( id );
      poss.push_back(pos );
      chroms.push_back(chrom);
      indexes.push_back( maprecords++ );
    }
  //H5 data spaces
  hsize_t dims[1];
  dims[0] = maprecords;

  DSetCreatPropList ds_creatplist;  // create dataset creation prop list
  ds_creatplist.setChunk( 1, dims );  // then modify it for compression
  ds_creatplist.setDeflate( 9 );

  DataSpace * dataspace = new DataSpace(1,dims);

  DataSet * dataset = new DataSet(file->createDataSet("/Markers/positions",
						      PredType::NATIVE_INT,
						      *dataspace,
						      ds_creatplist));
  dataset->write( &*poss.begin(),PredType::NATIVE_INT );

  delete dataset;

  dataset = new DataSet(file->createDataSet("/Markers/indexes",
					    PredType::NATIVE_INT,
					    *dataspace,
					    ds_creatplist));

  dataset->write( &*indexes.begin(), 	    
		  PredType::NATIVE_INT );

  delete dataset;

  const char * stringtemp[maprecords];
  for( unsigned i = 0 ; i < maprecords ; ++i )
    {
      stringtemp[i] = ids[i].c_str();
    }

  //delete dataspace;

  //dataspace = new  new DataSpace(1,dims);
  dataset = new DataSet(file->createDataSet("/Markers/ids",
					    //H5Tcopy(H5T_C_S1),
					    StrType(H5::PredType::C_S1, H5T_VARIABLE),
					    //PredType::C_S1,
					    //PredType::NATIVE_CHAR,
					    *dataspace,
					    ds_creatplist));

  cerr << stringtemp[0] << ' ' << stringtemp[1] << '\n';
  dataset->write(stringtemp,
		 StrType(H5::PredType::C_S1, H5T_VARIABLE));
		 //PredType::C_S1);//PredType::NATIVE_CHAR);
  //for( 

  delete dataset;

  dataset = new DataSet(file->createDataSet("/Markers/chromos",
					    //H5Tcopy(H5T_C_S1),
					    StrType(H5::PredType::C_S1, H5T_VARIABLE),
					    //PredType::C_S1,
					    //PredType::NATIVE_CHAR,
					    *dataspace,
					    ds_creatplist));

 for( unsigned i = 0 ; i < maprecords ; ++i )
    {
      stringtemp[i] = chroms[i].c_str();
    }

  dataset->write(stringtemp,
		 StrType(H5::PredType::C_S1, H5T_VARIABLE));

  delete group;

  //Now, go through the perms

  group = new Group(file->createGroup("/Perms"));
  unsigned long permno=1;
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
	  if ( rv != -1 && rv != 0)
	    {
	      ostringstream gname;
	      gname << "/Perms/perm" << permno++;

	      dataset = new DataSet(file->createDataSet(gname.str().c_str(),
							PredType::NATIVE_DOUBLE,
							*dataspace,
							ds_creatplist));
	      dataset->write( &*perms.begin(),	
			      PredType::NATIVE_DOUBLE );
	      delete dataset;
	    }
	}
      while(! gzeof( gzin ) );
      gzclose(gzin);
    }
  
  delete group;
  delete dataspace;

  file->close();
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
