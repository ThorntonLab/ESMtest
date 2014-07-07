//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//other boost stuff
#include <boost/bind.hpp>

#include <H5Cpp.h>


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


using namespace std;
using namespace boost::program_options;
using namespace H5;

struct options
/*
  This object represents the command-line options
 */
{
  bool strip,convert;
  string bimfile,infile,outfile;
  size_t nrecords;
  options(void);
};

options::options(void) : strip(true),
			 convert(true),
			 infile(string()),
			 outfile(string()),
			 nrecords(1)
{
}

options process_argv( int argc, char ** argv );
size_t process_bimfile( const options & O, H5File & ofile );
void process_perms( const options & O, size_t nmarkers, H5File & ofile );

int main( int argc, char ** argv )
{
  options O = process_argv( argc, argv );

  //Create output file
  H5File ofile( O.outfile.c_str() , H5F_ACC_TRUNC );
  size_t nmarkers = process_bimfile( O, ofile );
  process_perms( O, nmarkers, ofile );
  ofile.close();
  exit(0);
}

options process_argv( int argc, char ** argv )
{
  options rv;

  options_description desc("Reads PLINK permuted data table from IFP.  Reformats to a HDF5 file.");
  desc.add_options()
    ("help,h", "Produce help message")
    ("nostrip","Do not strip the observed data from the file (default is to strip)")
    ("noconvert","Do not convert input into a p-value.  Default is to assume that the input is a chi^2 statistic with 1 degree of freedom")
    ("bim,b",value<string>(&rv.bimfile)->default_value(string()),"The bim file (map file for binary PLINK data)")
    ("infile,i",value<string>(&rv.infile)->default_value(string()),"Input file name containing permutations.  Default is to read from stdin")
    ("outfile,o",value<string>(&rv.outfile)->default_value(string()),"Output file name.  Format is HDF5")
    ("nrecords,n",value<size_t>(&rv.nrecords)->default_value(1),"Number of records to buffer.")
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

  if (! vm.count("bim") )
    {
      cerr << "Error, no bim file name specified\n";
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

size_t process_bimfile( const options & O, H5File & ofile )
/*
  Write map data into an h5 group called "Markers"
 */
{
  ifstream bimin( O.bimfile.c_str() );
  if (! bimin )
    {
      cerr << "Error, " << O.bimfile
	   << " could not be opened for reading\n";
      exit(10);
    }

  vector<string> markers,chroms;
  vector< const char * > marker_str,chrom_str;
  vector< int > vpos;
  string chrom,marker,dummy,line;
  int pos;
  while( !bimin.eof() )
    {
      bimin >> chrom >> marker >> dummy >> pos >> ws;
      getline( bimin, line );
      bimin >> ws;
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

    //vector< double > data(10*nmarkers);
    // double ** data = new double *[O.nrecords];
    // for( size_t i = 0 ; i < O.nrecords ; ++i )
    //   {
    // 	data[i] = new double[nmarkers];
    //   }
    //double * data = new double(O.nrecords*nmarkers);
    vector<double> data(O.nrecords*nmarkers);
    int repno;
    fscanf(ifp,"%d",&repno);
    //The first line is the observed data
    //for( size_t i = 0 ; i < nmarkers-1 ; ++i )
    for( size_t i = 0 ; i < nmarkers ; ++i )
      {
	int rv = fscanf(ifp,"%lf",&data[i]);
	if(O.convert)
	  {
	    data[i] = (data[i]!=1.) ? -log10(gsl_cdf_chisq_Q(data[i],1.)) : 0.;
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
    hsize_t chunk_dims2[2] = {10, nmarkers};
    hsize_t maxdims2[2] = {H5S_UNLIMITED,nmarkers};
    hsize_t datadims[2] = {0,nmarkers};
    hsize_t offsetdims[2] = {0,0};
    hsize_t recorddims[2] = {O.nrecords,nmarkers};
    cparms.setChunk( 2, chunk_dims2 );
    cparms.setDeflate( 6 );

    //dataspace = new DataSpace(2, recorddims, maxdims2);
    DataSpace fspace(2,datadims,maxdims2);
    d = new DataSet(ofile.createDataSet("/Perms/permutations",
					PredType::NATIVE_DOUBLE,
					fspace,
					cparms));

    DataSpace memspace(2,recorddims);
    //ofile.close();exit(10);
    //while (fscanf(ifp,"%d",&repno) != -1 )
    while(!feof(ifp))
      {
	int rv = 1;
	size_t I = 0;
	for( size_t i = 0 ; (rv!=-1&&rv!=0) && i < O.nrecords ; ++i )
	  {
	    repno = -1;
	    rv=fscanf(ifp,"%d",&repno);
	    if( rv == 0 || rv == -1 || feof(ifp ) )
	      {
		return; //we have hit the end of the file
	      }
	    //for( size_t j = 0 ; j < nmarkers-1 ; ++j,++I )
	    for( size_t j = 0 ; j < nmarkers ; ++j,++I )
	      {
		rv = fscanf(ifp,"%lf",&data[I]);
		if(rv==0||rv==-1)
		  {
		    cerr << "Error, input stream ended before expected...\n";
		    ofile.close();
		    exit(10);
		  }
		if(O.convert)
		  {
		    data[I]=-log10(gsl_cdf_chisq_Q(data[I],1.));
		  }
	      }
	  }
	datadims[0] += O.nrecords;
	d->extend( datadims );
	DataSpace * dspace = new DataSpace( d->getSpace() );
	dspace->selectHyperslab(H5S_SELECT_SET,recorddims,offsetdims);
	//d->write(data.data(), PredType::NATIVE_DOUBLE,memspace,*dspace);
	d->write(data.data(), PredType::NATIVE_DOUBLE,memspace,*dspace);
	delete dspace;
	offsetdims[0] += O.nrecords;
      }
}
