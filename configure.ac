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

echo "
VRECO_CPMD configuration summary:
=================================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 flags:           $FCFLAGS
  Fortran77 compiler:        $F77
  Fortran77 flags:           $FFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIBS_EXT_DEPS $LIBS

Build:
  Dynamic library:           $ax_build_dynamic

Options:
  With OpenMP:               $ax_fc_openmp_msg
  With dgemmsy support:      $ax_have_dgemmsy
  With libyaml:              $ax_have_yaml
   | C flags:                $LIB_YAML_CFLAGS
   | linked libraries:       $LIB_YAML_LIBS

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Static libraries:          $libdir
  Fortran modules:           $moduledir
  Documentation:             $docdir
"

if test x"$ac_use_mpi" = x"warn" -o x"$withscalapack" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
  if test x"$withscalapack" = x"warn" ; then
    echo "  Scalapack binding:         disabled by tests"
  fi

fi
