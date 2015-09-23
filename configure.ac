AC_PREREQ(2.63)
AC_INIT([ESM],[0.1.0],[https://github.com/ThorntonLab/ESMtest])
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX
AC_CONFIG_MACRO_DIR([m4])

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
AC_LANG(C++)
AC_C_INLINE
AC_FUNC_MALLOC
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T


AC_CHECK_HEADER(Sequence/SimData.hpp, , [AC_MSG_ERROR([Sequence/SimData.hpp not found. libsequence (www.molpopgen.org) >= 1.7.8 required.])])
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5 is requred])])
AC_CHECK_HEADER(H5Cpp.h,,[AC_MSG_ERROR([H5Cpp.h not found.  HDF5 >= 1.8.11 is requred])])
AC_CHECK_HEADER(boost/program_options.hpp,[AC_DEFINE([HAVE_BOOST_PROGRAM_OPTIONS],[1],[Is boost program options header found?])],[AC_MSG_ERROR([boost program options requested but boost/program_options.hpp not found])])
AC_CHECK_HEADER(boost/bind.hpp,[AC_DEFINE([HAVE_BOOST_BIND],[1],[Is boost bind header found?])],[AC_MSG_ERROR([boost bind requested but boost/bind.hpp not found])])
AC_CHECK_HEADER(boost/algorithm/string.hpp,[AC_DEFINE([HAVE_BOOST_ALGORITHM_STRING],[1],[Is boost algorithm string header found?])],[AC_MSG_ERROR([boost algorithm string requested but boost/algorithm/string.hpp not found])])
AC_CHECK_HEADER(boost/unordered_map.hpp,[AC_DEFINE([HAVE_BOOST_UNORDERED_MAP],[1],[Is boost unordered map header found>])],[AC_MSG_ERROR([boost map requested but boost/unordered_map.hpp not foun])])
AC_CHECK_HEADER(gsl/gsl_cdf.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_cdf.h not found.  Please install the GNU Scientific Library])])


dnl check for C run-time libraries
AC_CHECK_LIB([z],gzungetc,,[AC_MSG_ERROR([zlib run time library not found])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR([gslcblas run-time library not found])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([gsl run-time library not found])])
AC_CHECK_LIB([hdf5],[H5Dget_type],,[AC_MSG_ERROR([HDF5 run-time library not found])])

dnl check for C++ run-time libraries
AC_LANG_SAVE
	AC_LANG_CPLUSPLUS
	AC_REQUIRE_CPP
	AC_CHECK_LIB(sequence, main,,[AC_MSG_ERROR([libsequence runtime library does not seem to be here])])	    
	AC_CHECK_LIB(boost_system,main,,[AC_MSG_WARN([boost containers requested but the boost_system runtime library does not seem to be here])])
	AC_CHECK_LIB(boost_program_options,main,,[AC_MSG_WARN([boost program options requested but the boost_program_option runtime library does not seem to be here])])
	AC_CHECK_LIB([hdf5_cpp],main,,[AC_MSG_ERROR([HDF5 C++ run-time library not found])])
AC_LANG_RESTORE

AC_OUTPUT
