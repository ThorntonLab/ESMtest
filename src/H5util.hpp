#ifndef __H5util_HPP__
#define __H5util_HPP__

#include <H5Cpp.h>
#include <vector>
#include <string>
#include <ESMH5type.hpp>

std::vector< std::string > read_strings( const char * filename, 
					      const char * dsetname );

std::vector<int> read_ints( const char * filename, 
			    const char * dsetname );

std::vector<ESMBASE> read_doubles(const char * filename, 
				 const char * dsetname );

std::vector<ESMBASE> read_doubles_slab(const char * filename,
				      const char * dsetname,
				      const size_t & start,
				      const size_t & len);

void write_strings( const std::vector<std::string> & data,
			 const char * dsetname,
			 H5::H5File ofile );

void write_ints ( const std::vector<int> & data ,
		  const char * dsetname,
		  H5::H5File ofile );

void write_doubles ( const std::vector<double> & data ,
		     const char * dsetname,
		     H5::H5File ofile );

#endif
