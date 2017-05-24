#   Name of the fortran compiler
#-------------PGI --------------
#FC= pgf90 -O2 
#-------------PGI OpenMP --------------
#FC= pgf90 -O2 -mp
#-------------G95----------------
#FC= g95 -O2 -funroll-loops
#-------------GFORTRAN-----------
#FC= gfortran -Wall -O2 -funroll-loops
#-------------GFORTRAN OpenMP ---
#FC= gfortran -Wall -O2 -fopenmp -Wall -funroll-loops
#-------------IFORT--------------
#FC= ifort -i-static -O2 -unroll 
#-------------IFORT OpenMP-------
FC= ifort -i-static -O2 -unroll -openmp 
#-------------IBM-RISC/PPC-------
#FC= xlf95 -O2 -qarch=auto 
#-------------IBM-RISC/PPC OpenMP-------
#FC= xlf95_r -O2 -qarch=auto -qsmp=omp 
#-------------Sun Studio Linux ---------
#FC= sunf95 -O2 
#-------------Sun Studio Linux OpenMP---------
#FC= sunf95 -xopenmp=parallel -xO3 
#################################

default: Vreco_CPMD.x mtd_restart_extract.x

mtd_restart_extract.x : mtd_restart_extract.f90
	$(FC) -o  $@  $^

Vreco_CPMD.x : Vreco_CPMD.f90
	$(FC) -o  $@  $^

#
gtar: clean
	@( d=`date +%Y%m%d` ; b=`basename $$PWD` ; cd .. ;		   \
	tar -c --exclude \*,v --exclude \*.bak --exclude \*~ --exclude CVS \
	--exclude \*.o --exclude \*.a --exclude \*.x --exclude \*.log	   \
        --exclude \*.out --exclude \*.prj --exclude \*.chk 		   \
	--exclude \*.orig --exclude \*.rej -zvvf $$b-$$d.tar.gz $$b &&	   \
        echo successfully created $$b-$$d-tar.gz ; cd $$b )

diff:
	rcsdiff -u *,v

clean :
	rm -f *.mod *.o *.L *~ Vreco_CPMD.x mtd_restart_extract.x
