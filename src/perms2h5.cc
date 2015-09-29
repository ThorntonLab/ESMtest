//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//other boost stuff
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

#include <H5Cpp.h>
#include <ESMH5type.hpp>

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
#include <cctype>
#include <cstring>
#include <typeinfo>

using namespace std;
using namespace boost::program_options;
using namespace H5;

struct options
/*
  This object represents the command-line options
 */
{
  bool strip,convert,verbose,compression,dbprec,nochunk;
  string bimfile,ldfile,infile,outfile;
  size_t nrecords,ccache,cmarkers;
  options(void);
};

options::options(void) : strip(true),
			 convert(true),
			 verbose(false),
			 compression(false),
			 dbprec(false),
			 nochunk(false),
			 bimfile(string()),
			 ldfile(string()),
			 infile(string()),
			 outfile(string()),
			 nrecords(1),
			 ccache(5),
			 cmarkers(50)
{
}

options process_argv( int argc, char ** argv );
size_t process_bimfile( const options & O, H5File & ofile );
void process_ldfile( const options & O, H5File & ofile );
void process_perms( const options & O, size_t nmarkers, H5File & ofile );
void firstprime ( size_t & num );
int main( int argc, char ** argv )
{
  options O = process_argv( argc, argv );
  
  //Create output file
  size_t cache_bytes = O.ccache*1024*1024; //param in mb -> b
  size_t chunk_dat = O.nrecords*O.cmarkers;//number of total entries in a chunk
  size_t chunk_bytes = chunk_dat*4;  
  if ( O.dbprec ) 
    {
      size_t chunk_bytes = chunk_dat*8;	
    }
  size_t nc_cache =  cache_bytes/chunk_bytes + 1 ; 
  if ( nc_cache == 1 )
    {
      cerr << "Raw data cache size smaller than chunk size, either increase cache, decrease chunk size or disable chunking entirely" << '\n';
      exit(0);
    }
  //find the first prime number which is bigger than 10*nc_cache
  size_t rdcc = 10*nc_cache;

  firstprime(rdcc);
  
  FileAccPropList fapl;
  fapl.setCache(23,rdcc,cache_bytes,0) ;
  //Number of elements in meta data cache(small prime)
  //Number of elements in the raw data chunk cache size (should be prime number 10-100 * NCHUNKS/cache)
  //O.ccache command arg set to 5MB as default( H5 default = 1MB)
  //Preemption policy set to no preemption
  H5File ofile( O.outfile.c_str() , H5F_ACC_TRUNC,H5P_DEFAULT,fapl );
  size_t nmarkers = process_bimfile( O, ofile );
  process_perms( O, nmarkers, ofile );
    if ( O.verbose )
    {
      cerr << "I finished processing perms" <<"\n";
    }

  process_ldfile( O, ofile );
  if ( O.verbose )
    {
      cerr << "I finished processing LD" <<"\n";
    }

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
    ("linkage,l",value<string>(&rv.ldfile)->default_value(string()),"The LD file (pairwise r^2 from PLINK)")
    ("infile,i",value<string>(&rv.infile)->default_value(string()),"Input file name containing permutations.  Default is to read from stdin")
    ("outfile,o",value<string>(&rv.outfile)->default_value(string()),"Output file name.  Format is HDF5")
    ("nrecords,n",value<size_t>(&rv.nrecords)->default_value(1),"Number of records to buffer.")
    ("ccache,a",value<size_t>(&rv.ccache)->default_value(5),"Raw data chunk cache in mega bytes(will be converted to bytes), default = 5MB")
    ("cmarkers,m",value<size_t>(&rv.cmarkers)->default_value(50),"Number of markers in a chunk")
    ("compression,c","Gzip level 6 + shuffle compression")
    ("nochunk","Chunked storage, default is true, false=contiguous")
    ("dbprec","ESM base type set to double")
    ("verbose,v","Write process info to STDERR")
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

  if( vm.count("verbose") )
    {
      rv.verbose = true;
    }
  if( vm.count("compression") )
    {
      rv.compression = true;
    }
  if (vm.count("nochunk"))
    {
      rv.nochunk = true;
    }
  if (vm.count("dbprec"))
    {
      rv.dbprec = true;
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
  if ( O.verbose )
    {
      cerr << ".bim file contains " << markers.size() << " markers\n";
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

  if (O.compression )
    {
      cparms.setDeflate( 6 ); //compression level makes a big differences in large files!  Default is 0 = uncompressed.
    }
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

    if ( O.verbose )
      {
	cerr << "Starting to process perms for "
	     << nmarkers << " markers\n";
      }
    
    vector<ESMBASE> data(O.nrecords*nmarkers);
    int repno;
    fscanf(ifp,"%d",&repno);

    char * buffer = new char[100];
    //The first line is the observed data
    for( size_t i = 0 ; i < nmarkers ; ++i )
      {
	int rv = fscanf(ifp,"%f",&data[i]);
	if(O.convert)
	  {
	    data[i] = (data[i]!=1.) ? -log10(gsl_cdf_chisq_Q(data[i],1.)) : 0.;
	  }
      }

    if ( O.verbose )
      {
	cerr << "Read in observed data.\n" << endl;
      }
    ofile.createGroup("/Perms");

    DSetCreatPropList cparms;
    hsize_t chunk_dims[1] = {nmarkers};
    hsize_t maxdims[1] = {nmarkers};
    
    cparms.setChunk( 1, chunk_dims );
    if ( O.compression)
      {
	cparms.setShuffle();//try this for better compression permformance
	cparms.setDeflate( 6 ); //compression level makes a big differences in large files!  Default is 0 = uncompressed.
      }
  
    DataSpace * dataspace = new DataSpace(1,chunk_dims, maxdims);

    DataSet * d = new DataSet(ofile.createDataSet("/Perms/observed",
						  PredType::NATIVE_FLOAT,
						  *dataspace,
						  cparms));

    d->write( data.data(), PredType::NATIVE_FLOAT );

    delete dataspace;
    delete d;

    //ok, now we write a big matrix of the permuted values
    hsize_t chunk_dims2[2] = {O.nrecords,O.cmarkers};//{10, nmarkers};
    hsize_t maxdims2[2] = {H5S_UNLIMITED,nmarkers};
    hsize_t datadims[2] = {0,nmarkers};
    hsize_t offsetdims[2] = {0,0};
    hsize_t recorddims[2] = {O.nrecords,nmarkers};

    cparms.setChunk( 2, chunk_dims2 );

    DataSpace fspace(2,datadims,maxdims2);
    d = new DataSet(ofile.createDataSet("/Perms/permutations",
					PredType::NATIVE_FLOAT,
					fspace,
					cparms));

    
    size_t RECSREAD = 0;
    while(!feof(ifp))
      {
	int rv = 1;
	size_t I = 0;
	RECSREAD = 0;
	for( size_t i = 0 ; !feof(ifp) && i < O.nrecords ; ++i,++RECSREAD )
	  {
	    if ( O.verbose )
	      {
		cerr << i << ": ";
	      }
	    repno = -1;
	    rv=fscanf(ifp,"%d",&repno);
	    if( rv == 0 || rv == -1 || feof(ifp ) )
	      {
	    	if( O.verbose ) cerr << "return 1\n";
	    	return; //we have hit the end of the file
	      }
	    
	    for( size_t j = 0 ; j < nmarkers ; ++j,++I )
	      {
		if(O.verbose) cerr << j << ',' << I << ' ';
		rv = fscanf(ifp,"%f",&data[I]);
	
		if(O.convert)
		  {
		    data[I]= (data[I] != 1.) ? -log10(gsl_cdf_chisq_Q(data[I],1.)) : 0.;
		  }
	      }
	    if(O.verbose) cerr<< endl;
	  }
	datadims[0] += RECSREAD;
	recorddims[0] = RECSREAD;
	DataSpace memspace(2,recorddims);
	d->extend( datadims );
	DataSpace * dspace = new DataSpace( d->getSpace() );
	dspace->selectHyperslab(H5S_SELECT_SET,recorddims,offsetdims);
	d->write(data.data(), PredType::NATIVE_FLOAT,memspace,*dspace);
	delete dspace;
	offsetdims[0] += O.nrecords;
      }
}

void process_ldfile( const options & O, H5File & ofile )
{
  if ( O.verbose )
    {
      cerr << "I am processing LD" <<"\n";
    }
  ifstream ldf (O.ldfile.c_str());
  if ( O.verbose )
    {
      cerr << "I opened LD file" <<"\n";
    }
  string line;
  vector<string> lineVector,snpA,snpB;
  vector < const char * > snpA_str, snpB_str;
  vector<ESMBASE> rsq;
  ESMBASE r2;
  string::size_type sz;
  if (ldf.is_open())
    {
      if ( O.verbose )
    {
      cerr << "LD file is open" <<"\n";
    }
      while(getline(ldf,line))
	{

	  boost::split(lineVector,line,boost::is_any_of(" "),boost::token_compress_on);
	  snpA.push_back(lineVector.at(3));
	  if ( O.verbose )
    {
      cerr << "SNP A = " << lineVector.at(3).c_str() << "\n";
      cerr << "snpA type is  " << typeid(lineVector.at(3).c_str()).name()<<"\n";
    }
	  snpB.push_back(lineVector.at(6));
	  r2 = atof( lineVector.at(7).c_str());
	  rsq.push_back(r2);
	
	}
        if ( O.verbose )
    {
      cerr << "I finished reading ld" <<"\n";
    }
      ldf.close();
        if ( O.verbose )
    {
      cerr << "I closed LD" <<"\n";
    }
    }
  else cerr <<"Unable to read LD file"<< '\n';
    for( unsigned i = 0 ; i < snpA.size() ; ++i )
    {
      snpA_str.push_back( snpA[i].c_str() );
      snpB_str.push_back( snpA[i].c_str() );
    }
  ofile.createGroup("/LD");
  if ( O.verbose )
    {
      cerr << "I made the group" <<"\n";
    }
  DSetCreatPropList cparms;
  hsize_t chunk_dims[1] = {snpA_str.size()};
  hsize_t maxdims[1] = {snpA_str.size()};

  cparms.setChunk( 1, chunk_dims );
  cparms.setDeflate( 6 ); //compression level makes a big differences in large files!  Default is 0 = uncompressed.
  
  DataSpace dataspace(1,chunk_dims, maxdims);

  H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE); 

  if ( O.verbose )
    {
      cerr << "I setup the space" <<"\n";
    }
  DataSet snpA_dset = ofile.createDataSet("/LD/snpA",
					    datatype,
					    dataspace,
					    cparms);
if ( O.verbose )
    {
      cerr << "snpA type is  " << typeid(snpA).name()<<"\n";
      //cerr << "datatype is " << datatype <<"\n";
    }
  snpA_dset.write(snpA_str.data(), datatype );
  if ( O.verbose )
    {
      cerr << "wrote SNP A" <<"\n";
    }
  DataSet snpB_dset = ofile.createDataSet("/LD/snpB",
					    datatype,
					    dataspace,
					    cparms);

  snpB_dset.write(snpB_str.data(), datatype );
  if ( O.verbose )
    {
      cerr << "wrote SNP B" <<"\n";
    }
  DataSet rsq_dset = ofile.createDataSet("/LD/rsq",
					 H5::PredType::NATIVE_FLOAT,
					 dataspace,
					 cparms);


  rsq_dset.write( rsq.data(),H5::PredType::NATIVE_FLOAT );
  if ( O.verbose )
    {
      cerr << "wrote Rsq" <<"\n";
    }
}

void firstprime (size_t & num)
{
  size_t count = 0;
  bool prime = false;
  while ( prime == false )
    {
      count = 0;

      for (size_t i=2;i<num;i++)
        {

          if (num%i==0)
            {
              count++;
            }
        }
      if ( count == 0 )
        {
          prime = true;
        }
      else
        {
          num++;
        }
    }
}
