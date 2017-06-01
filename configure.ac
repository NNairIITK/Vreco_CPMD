#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([VRECO_CPMD], [10.3.0], [nnair@iitk.ac.in])

# All temporary building files are put in this directory.
# This must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)

# Use a config.h file to store macro definitions, required for CUDA.
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([config/m4])

# Initialize
AM_INIT_AUTOMAKE


# Checks for programs.
AC_PROG_FC([], [Fortran 90])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# Define Fortran for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

AC_FC_WRAPPERS()
AX_FC_DEFS()
AX_FC_OPENMP()


# Set the FCFLAGS and FFLAGS variable
if test -z "$FFLAGS" ; then
  AC_SUBST(FFLAGS, "-O2")
fi
if test -n "$FCFLAGS_f90" ; then
  FCFLAGS="$FCFLAGS $FCFLAGS_f90"
fi
if test -z "$FCFLAGS" ; then
  AC_SUBST(FCFLAGS, "-O2")
fi

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
