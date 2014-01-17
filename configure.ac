AC_INIT([context], [1.0], [hickey@mcb.mcgill.ca]) 
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_PROG_RANLIB
AM_PROG_AR
#AC_CHECK_LIB([gsl], [gsl_linalg_exponential_ss], [], [echo "Error! gsl not found.  Make sure GSL is installed"
#        exit -1
#        ])
#AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [echo "Error! gslcblas not found.  Make sure GSL is installed"
#        exit -1
#        ])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
		  Makefile
		  src/Makefile
])
AC_OUTPUT

