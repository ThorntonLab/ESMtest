#include <H5util.hpp>
#include <cstdlib>

using namespace std;
using namespace H5;

static const size_t MAXSTRINGSIZE=1000;

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

vector<double> read_doubles_slab( const char * filename, 
				  const char * dsetname,
				  const size_t & start,
				  const size_t & len)
{
  H5File ifile( filename, H5F_ACC_RDONLY );
  DataSet ds( ifile.openDataSet(dsetname) );
  DataSpace dsp(ds.getSpace());
  int rank_j = dsp.getSimpleExtentNdims();
  hsize_t dims_out[rank_j];
  int ndims = dsp.getSimpleExtentDims( dims_out, NULL);
  
  //*Define the hyperslab in the dataset; see readdata.cpp in 
  // the HDF5 group c++ API
  hsize_t offset[2];
  hsize_t count[2];
  offset[0]= 0;
  offset[1]= start;
  count[0] = dims_out[0];
  count[1]= len;
  //should select a hyperslab which ds.read can reference
  dsp.selectHyperslab(H5S_SELECT_SET,count,offset);
  
  //define memspace

  hsize_t dimsm[2];
  dimsm[0]=dims_out[0];
  dimsm[1]=len;
  DataSpace memspace(2,dimsm);

  //define hyperslab in memory...this is done so you can go from
  //high dimensions to lower ones while still in H5

  hsize_t offset_out[2]; //mem offset
  hsize_t count_out[2]; //

  offset_out[0] = 0;
  offset_out[1] = 0;
  count_out[0]  = dims_out[0];
  count_out[1]  = len;
  
  memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );
  

  vector<double> receiver(dims_out[0]*len); //allocate memory to receive
  IntType intype = ds.getIntType();
  ds.read( &receiver[0], intype, memspace, dsp);
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
