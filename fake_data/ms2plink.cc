#include <Sequence/SimData.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdlib>

using namespace std;
using namespace Sequence;

int main( int argc, char ** argv )
{
  int argn = 1;
  const char * pedfile = argv[argn++];
  const char * mapfile = argv[argn++];

  SimData d;
  d.fromfile(stdin);

  ostringstream pedbuffer,mapbuffer;

  //buffer the ped file
  bool sex = 0;
  for( unsigned i = 0 ; i < d.size() ; i += 2 )
    {
      pedbuffer << 0 << ' ' //family ID
		<< (i/2+1) << ' ' //individual ID
		<< 0 << ' ' //paternal ID
		<< 0 << ' ' //maternal ID
		<< (sex+1) << ' ' //sex ID, 1 = male, 2 = female
		<< ( (i >= d.size()/2) + 1 ) << ' '; //phenotype.  1 = control, 2 = case
      sex = !sex;
      for(unsigned j=0;j<d.numsites(); ++j)
	{
	  pedbuffer << ((d[i][j] == '0') ? 'A' : 'T') << ' ';
	  pedbuffer << ((d[i+1][j] == '0') ? 'A' : 'T')<< ' ';
	}
      pedbuffer << '\n';
    }

  //make MAP file buffer
  unsigned SNP=1;
  for( SimData::const_pos_iterator p = d.pbegin() ; p !=d.pend() ; ++p )
    {
      mapbuffer << 1 << ' ' //chromo
		<< "rs" << SNP++ << ' ' //SNP label
		<< int( double(100000)**p) << '\n'; //position
    }

  ofstream pedout(pedfile),mapout(mapfile);
  pedout << pedbuffer.str();
  mapout << mapbuffer.str();

  pedout.close();
  mapout.close();
  exit(10);
}
