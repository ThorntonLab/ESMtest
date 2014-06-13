//Command line parsing using boost (C++)
#include <boost/program_options.hpp>

//Used to check ZLIB version number during compile time
#include <boost/static_assert.hpp>

//Headers for gzip output (C language)
#include <zlib.h>

//hdf5 library
#include <H5Cpp.h>

//Headers to conver chi-squared statistic into chi-squared p-value.  GNU Scientific Library (C language)
//#include <gsl/gsl_cdf.h>

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
//BOOST_STATIC_ASSERT(ZLIB_VERNUM >= 0x1250);

//for reading in variable-length strings
static const size_t MAXSTRINGSIZE=1000;

using namespace std;
using namespace boost::program_options;
using namespace H5;

struct options
/*
  This object represents the command-line options
 */
{
  string outfile;
  vector<string> infiles;
  bool convert;
  //constructor for empty object
  options(void);
};

options::options(void) : outfile(string()),
			 infiles( vector<string>() ),
			 convert(false)
{
}

struct markerdata
{
  vector<string> ids,chroms;
  vector<int> pos;
  markerdata( const vector<string> & i,
	      const vector<string> & c,
	      const vector<int> p ) : ids(i), 
				      chroms(c),
				      pos(p)
  {
  }
  markerdata() : ids( vector<string>() ),
		 chroms( vector<string>() ),
		 pos( vector<int>() )
  {
  }
};

//process command lines
options process_argv( int argc, char ** argv );
//Reads in the marker data from a file, and writes them to output file
markerdata process_markers( const options & O, H5File ofile );
//Reads in a data set of variable-length strings, where each string has length <= MAXSTRINGSIZE
vector< string > read_strings( const char * filename, const char * dsetname );
vector<int> read_ints( const char * filename, const char * dsetname );
vector<double> read_doubles(const char * filename, const char * dsetname );
void write_strings( const std::vector<string> & data,
		    const char * dsetname,
		    H5File ofile );
void write_ints ( const vector<int> & data ,
		  const char * dsetname,
		  H5File ofile );
void write_doubles ( const vector<double> & data ,
		     const char * dsetname,
		     H5File ofile );
void make_permset(const size_t & nmarkers, H5File ofile);

void copy_perms( const char * source,
		 H5File ofile );

int main( int argc, char ** argv )
{
  options O = process_argv( argc, argv );

  //Open our hdf5 output file
  H5File ofile( O.outfile.c_str() , H5F_ACC_TRUNC );

  markerdata md = process_markers( O, ofile );
  vector<double> observed = read_doubles( O.infiles[0].c_str(),"/Perms/observed");
  ofile.createGroup("/Perms");
  write_doubles(observed,"/Perms/observed",ofile);
  make_permset( md.ids.size(), ofile);
  for( unsigned i = 0 ; i < O.infiles.size() ; ++i )
    {
      copy_perms(O.infiles[i].c_str(), ofile);
    }
  ofile.close();
}

options process_argv( int argc, char ** argv )
{
  options rv;

  options_description desc("Merge permuted statistics into table for ESM test");
  desc.add_options()
    ("help,h", "Produce help message")
    ("outfile,o",value<string>(&rv.outfile)->default_value(string()),"Output file name.  Format is binary and gzipped")
    ("convert,c",value<bool>(&rv.convert)->default_value(false),"Convert chi-squared to p-value")
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

vector< string > read_strings( const char * filename, const char * dsetname )
{
  hid_t file = H5Fopen( filename,H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t dset = H5Dopen (file, dsetname, H5P_DEFAULT);
  hid_t filetype = H5Dget_type (dset);

  hsize_t dims[1] = {MAXSTRINGSIZE};
  hid_t space = H5Dget_space (dset);
  int ndims = H5Sget_simple_extent_dims (space, dims, NULL);
  char **rdata = (char **) malloc (dims[0] * sizeof (char *));  
    
  hid_t memtype = H5Tcopy (H5T_C_S1);
  herr_t status = H5Tset_size (memtype, H5T_VARIABLE);

  status = H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);

  vector<string> rv;
  for (unsigned i=0; i<dims[0]; i++)
    {
      rv.push_back( string(rdata[i]) );
      free( rdata[i] );
    }
  free(rdata);

  H5Dclose(dset);
  H5Sclose(space);
  H5Tclose(filetype);
  H5Tclose(memtype);
  H5Fclose(file);
  return rv;
}

vector<int> read_ints( const char * filename, const char * dsetname )
{
  H5File ifile( filename, H5F_ACC_RDONLY );
  DataSet ds( ifile.openDataSet(dsetname) );
  DataSpace dsp(ds.getSpace());
  int rank_j = dsp.getSimpleExtentNdims();
  hsize_t dims_out[rank_j];
  int ndims = dsp.getSimpleExtentDims( dims_out, NULL);
  vector<int> receiver(dims_out[0]); //allocate memory to receive
  IntType intype = ds.getIntType();
  ds.read( &receiver[0], intype );
  return receiver;

  /*
  //C methods work -- for comparison
  hid_t file = H5Fopen( filename,H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t dset = H5Dopen (file, dsetname, H5P_DEFAULT);
  hid_t filetype = H5Dget_type (dset);

  hsize_t dims[1] = {MAXSTRINGSIZE};
  hid_t space = H5Dget_space (dset);
  int ndims;
  ndims = H5Sget_simple_extent_dims (space, dims, NULL);

  vector<int> rv( dims[0] );

  herr_t status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rv[0] );//srv.data());

  H5Fclose(dset);
  H5Dclose(file);
  return rv;
  */
}

vector<double> read_doubles( const char * filename, const char * dsetname )
{
  H5File ifile( filename, H5F_ACC_RDONLY );
  DataSet ds( ifile.openDataSet(dsetname) );
  DataSpace dsp(ds.getSpace());
  int rank_j = dsp.getSimpleExtentNdims();
  hsize_t dims_out[rank_j];
  int ndims = dsp.getSimpleExtentDims( dims_out, NULL);
  vector<double> receiver(dims_out[0]); //allocate memory to receive
  IntType intype = ds.getIntType();
  ds.read( &receiver[0], intype );
  return receiver;
}

void write_strings( const std::vector<string> & data,
		    const char * dsetname,
		    H5File ofile )
{
  vector< const char * > temp;
  for( unsigned i = 0; i < data.size() ; ++i )
    {
      temp.push_back( data[i].c_str() );
    }

  DSetCreatPropList cparms;
  hsize_t chunk_dims[1] = {data.size()};
  hsize_t maxdims[1] = {data.size()};
  
  cparms.setChunk( 1, chunk_dims );
  cparms.setDeflate( 6 );

  DataSpace dataspace(1,chunk_dims, maxdims);
  
  H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE); 
  
  DataSet dset = ofile.createDataSet(dsetname,
				     datatype,
				     dataspace,
				     cparms);
  
  dset.write(temp.data(), datatype );
}

void write_ints ( const vector<int> & data ,
		  const char * dsetname,
		  H5File ofile )
{
  DSetCreatPropList cparms;
  hsize_t chunk_dims[1] = {data.size()};
  hsize_t maxdims[1] = {data.size()};
  
  cparms.setChunk( 1, chunk_dims );
  cparms.setDeflate( 6 );

  DataSpace dataspace(1,chunk_dims, maxdims);
  
  DataSet dset = ofile.createDataSet(dsetname,
				     PredType::NATIVE_INT,
				     dataspace,
				     cparms);
  
  dset.write(data.data(), PredType::NATIVE_INT );
}

void write_doubles ( const vector<double> & data ,
		     const char * dsetname,
		     H5File ofile )
{
  DSetCreatPropList cparms;
  hsize_t chunk_dims[1] = {data.size()};
  hsize_t maxdims[1] = {data.size()};
  
  cparms.setChunk( 1, chunk_dims );
  cparms.setDeflate( 6 );

  DataSpace dataspace(1,chunk_dims, maxdims);
  
  DataSet dset = ofile.createDataSet(dsetname,
				     PredType::NATIVE_DOUBLE,
				     dataspace,
				     cparms);
  
  dset.write(data.data(), PredType::NATIVE_DOUBLE );
}

void make_permset(const size_t & nmarkers, H5File ofile)
{
 //Make a data set in the out file for the merged perms
  DSetCreatPropList cparms;

  hsize_t chunk_dims[2] = {10,nmarkers};
  hsize_t maxdims[2] = {H5S_UNLIMITED,nmarkers};
  
  cparms.setChunk( 2, chunk_dims );
  cparms.setDeflate( 6 );

  chunk_dims[0]=0;
  DataSpace dataspace(2,chunk_dims, maxdims);
  
  DataSet dset = ofile.createDataSet("/Perms/permutations",
				     PredType::NATIVE_DOUBLE,
				     dataspace,
				     cparms);
}

markerdata process_markers( const options & O, H5File ofile )
{
  if ( O.infiles.empty() ) return markerdata();

  //Get the marker data
  vector< string > markers = read_strings(O.infiles[0].c_str(),"/Markers/IDs");
  vector< string > chroms = read_strings(O.infiles[0].c_str(),"/Markers/chr");
  vector< int > pos = read_ints(O.infiles[0].c_str(),"/Markers/pos");

  //create group in ofile
  ofile.createGroup("/Markers");
  //Now write it to output file
  write_strings( markers, "/Markers/IDs", ofile );
  write_strings( chroms, "/Markers/chr", ofile );
  write_ints( pos, "/Markers/pos", ofile );


  return markerdata(markers,chroms,pos);
}

void copy_perms( const char * source,
		 H5File ofile )
{
  H5File ifile(source, H5F_ACC_RDONLY);

  //Set up properties of the source data set
  DataSet sdset( ifile.openDataSet("/Perms/permutations") );
  DataSpace sdspace( sdset.getSpace() );
  hsize_t source_dim[2];
  sdspace.getSimpleExtentDims( source_dim, NULL);
  hsize_t hslab_dims[2] = {1,source_dim[1]};
  vector<double> buffer( source_dim[1] );
  DataSpace memspace(2,hslab_dims),memspace2(memspace);

  //Set up dest
  DataSet ddset = ofile.openDataSet("/Perms/permutations");
  DataSpace ddspace( ddset.getSpace() );
  hsize_t dest_dim[2];
  ddspace.getSimpleExtentDims( dest_dim, NULL);

  hsize_t hslab_offset[2] = {dest_dim[0],0};
  hsize_t oset_size[2] = {dest_dim[0]+1,dest_dim[1]};
  cerr<< oset_size[0] << ' ' << oset_size[1] << '\n';
  for( hsize_t i = 0 ; i < source_dim[0] ; ++i,++hslab_offset[0],++oset_size[0] )
    {
      //read in the data from source
      sdspace.selectHyperslab(H5S_SELECT_SET, hslab_dims, hslab_offset);
      sdset.read(buffer.data(),PredType::NATIVE_DOUBLE,memspace,sdspace);

      ddset.extend(oset_size);
      ddspace = DataSpace(ddset.getSpace());
      ddspace.selectHyperslab(H5S_SELECT_SET,hslab_dims,hslab_offset);
      ddset.write( buffer.data(), PredType::NATIVE_DOUBLE, memspace2,ddspace);
    }
  cerr<< oset_size[0] << ' ' << oset_size[1] << '\n';
  ofile.close();exit(100);
}
