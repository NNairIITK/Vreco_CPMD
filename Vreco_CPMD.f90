! written by nisanth.nair@theochem.rub.de
! please report bugs and send comments to nisanth.nair@theochem.rub.de or cpmd-list@cpmd.org
! Version 10.3 OpenMP parallelization and clean up by Axel Kohlmeyer
MODULE kinds
  IMPLICIT NONE
  
  PUBLIC
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)
  CHARACTER (LEN=4), PARAMETER :: version = '10.4'
END MODULE kinds
!
MODULE pot_data
  USE kinds

  IMPLICIT NONE
  INTEGER         :: it,nd,dyn_frq
  LOGICAL         :: shft,tube
  REAL (KIND=dp)   :: rcut

  INTERFACE get_potential
     MODULE PROCEDURE get_potential
  END INTERFACE

  INTERFACE toproject
     MODULE PROCEDURE toproject
  END INTERFACE

  PRIVATE :: get_potential_dyn, get_potential_reg

CONTAINS
!    
  SUBROUTINE get_potential(so,scoll,is_torsion,scal,ds,s,w,prnt_dyn,ndyn,tmp,v)
    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(IN)           :: scoll(nd), ndyn
    REAL (KIND=dp), INTENT(IN  )  :: so(it,nd),scal(it,nd),ds(nd),s(nd),w(it)
    REAL (KIND=dp), INTENT(INOUT) :: v,tmp(ndyn)
    LOGICAL, INTENT(IN)           :: is_torsion(nd), prnt_dyn

    IF (prnt_dyn) THEN
       CALL get_potential_dyn(so,scoll,is_torsion,scal,ds,s,w,ndyn,tmp,v)
    ELSE
       CALL get_potential_reg(so,scoll,is_torsion,scal,ds,s,w,v)
    END IF
  END SUBROUTINE get_potential

  SUBROUTINE get_potential_dyn(so,scoll,is_torsion,scal,ds,s,w,ndyn,tmp,v)
      
    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(IN)           :: scoll(nd), ndyn
    REAL (KIND=dp), INTENT(IN)    :: so(it,nd),scal(it,nd),ds(nd),s(nd),w(it)
    REAL (KIND=dp), INTENT(INOUT) :: v,tmp(ndyn)
    LOGICAL, INTENT(IN)           :: is_torsion(nd)
    
    ! local variables
    INTEGER                       :: i, ii, k, imax
    REAL (KIND=dp)                :: ef1, ef2, ef3, sp, vtmp
    REAL (KIND=dp), ALLOCATABLE   :: dif1(:), dif2(:), s1(:), s2(:)
    REAL (KIND=dp), PARAMETER     :: pi=4.0_dp*ATAN(1.0_dp)
    
    v    = 0.0_dp

!$OMP parallel private(ii,i,k,imax,dif1,s1,s2,dif2,ef1,ef2,ef3,sp,vtmp) shared(tmp,v,s,scal,ds,w)
    ALLOCATE(dif1(nd),s1(nd),s2(nd),dif2(nd))
    
    DO ii=1,ndyn  ! loop over all blocks of mtd time steps

       vtmp = 0.0_dp
       imax = ii*dyn_frq
       IF (imax >= it) imax = it-1

!$OMP do schedule(static) 
       DO i=(ii-1)*dyn_frq+1,imax
          DO k=1,nd
             s1(k) = so(i,scoll(k))
             s2(k) = so(i+1,scoll(k))
          END DO
       
          dif1(1:nd) = s(1:nd)-s1(1:nd)
          dif1(1:nd) = dif1(1:nd)/scal(i,1:nd)
       
          DO k=1,nd
            IF (is_torsion(k)) THEN
              IF(dif1(k).GT.pi) THEN
                dif1(k) = dif1(k) - 2.0D0*pi
              ELSEIF(dif1(k).LT.-pi)THEN
                dif1(k) = dif1(k) + 2.0D0*pi
              ENDIF
            ENDIF
          ENDDO


          ef1 = DOT_PRODUCT(dif1(1:nd),dif1(1:nd))
       
          IF(.not.(shft .and. (SQRT(ef1) >= rcut*ds(i)) ) ) THEN
             
             ef1 = ef1/2.0_dp/ds(i)/ds(i)
             ef1 = EXP(-ef1) ! exponential factor 1
             
             dif1(1:nd) = s2(1:nd)-s1(1:nd)
             dif2(1:nd) = s(1:nd)-s1(1:nd)
             
             dif1(1:nd) = dif1(1:nd)/scal(i,1:nd)
             dif2(1:nd) = dif2(1:nd)/scal(i,1:nd)
          
             ef2 = DOT_PRODUCT(dif1(1:nd),dif2(1:nd))
             ef2 = ef2*ef2
             ef3 = 0.0_dp
             IF (shft) ef3 = EXP(-0.5_dp*(rcut)**2.0_dp)
          
             sp = DOT_PRODUCT(dif2(1:nd),dif2(1:nd))
             sp = 2.0_dp*sp*sp
          
             ef2 = EXP(-ef2/sp) ! exponential factor 2
             IF (.not.tube) ef2 = 1.0_dp
          
             vtmp = vtmp+w(i)*(ef1-ef3)*ef2      ! sum over all time steps in this block
          END IF

       END DO
!$OMP atomic
       v = v+vtmp
!$OMP barrier
!$OMP master
       tmp(ii) = v 
!$OMP end master

    END DO
    DEALLOCATE(dif1,dif2,s1,s2)
!$OMP end parallel

  END SUBROUTINE get_potential_dyn
!
!
  SUBROUTINE get_potential_reg(so,scoll,is_torsion,scal,ds,s,w,v)
      
    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(IN)           :: scoll(nd)
    REAL (KIND=dp), INTENT(IN)    :: so(it,nd),scal(it,nd),ds(nd),s(nd),w(it)
    REAL (KIND=dp), INTENT(INOUT) :: v
    LOGICAL, INTENT(IN)           :: is_torsion(nd)
    
    ! local variables
    INTEGER                       :: i, k
    REAL (KIND=dp)                :: ef1, ef2, ef3, sp
    REAL (KIND=dp), ALLOCATABLE   :: dif1(:), dif2(:), s1(:), s2(:)
    REAL (KIND=dp), PARAMETER     :: pi=4.0_dp*ATAN(1.0_dp)
    integer:: isss
    data isss /0/
    save ::  isss
    
    v    = 0.0_dp

!$OMP parallel private(dif1,s1,s2,dif2,i,k,ef1,ef2,ef3,sp) reduction(+:v)
    ALLOCATE(dif1(nd),s1(nd),s2(nd),dif2(nd))

!$OMP do schedule(static)
    DO i=1,it-1  ! loop over all mtd time steps
       
       DO k=1,nd
          s1(k) = so(i,scoll(k))
          s2(k) = so(i+1,scoll(k))
       END DO
       
       dif1(1:nd) = s(1:nd)-s1(1:nd)
       isss=isss+1
       DO k=1,nd
         IF (is_torsion(k)) THEN
           !if(isss.lt.2) print *, " k =", k, " is dihedral", "dif1(k)=", dif1(k), pi
           IF(dif1(k).GT.pi) THEN
             dif1(k) = dif1(k) - 2.0_dp*pi
             !if(isss.lt.2) print *, " dif1 > pi ", k, "dif1=", dif1(k), " pi =", pi
           ELSEIF(dif1(k).LT.-pi)THEN
             dif1(k) = dif1(k) + 2.0_dp*pi
             !if(isss.lt.2) print *, " dif1 < -pi", k, "dif1=", dif1(k), " pi =", pi
           ENDIF
         ENDIF
       ENDDO
       dif1(1:nd) = dif1(1:nd)/scal(i,1:nd)
       
       ef1 = DOT_PRODUCT(dif1(1:nd),dif1(1:nd))
       
       IF(.not.(shft .and. (SQRT(ef1) >= rcut*ds(i)) ) ) THEN
          
          ef1 = ef1/2.0_dp/ds(i)/ds(i)
          ef1 = EXP(-ef1) ! exponential factor 1
          
          dif1(1:nd) = s2(1:nd)-s1(1:nd)
          dif2(1:nd) = s(1:nd)-s1(1:nd)
          
          dif1(1:nd) = dif1(1:nd)/scal(i,1:nd)
          dif2(1:nd) = dif2(1:nd)/scal(i,1:nd)
          
          ef2 = DOT_PRODUCT(dif1(1:nd),dif2(1:nd))
          ef2=ef2*ef2
          ef3=0.0_dp
          
          IF (shft) ef3 = EXP(-0.5_dp*(rcut)**2.0_dp)
          
          sp = DOT_PRODUCT(dif2(1:nd),dif2(1:nd))
          sp = 2.0_dp*sp*sp
          
          ef2 = EXP(-ef2/sp) ! exponential factor 2
          
          IF (.not.tube) ef2 = 1.0_dp
          
          v = v+w(i)*(ef1-ef3)*ef2      ! sum over all time
       END IF
    END DO
    DEALLOCATE(dif1,dif2,s1,s2)
!$OMP end parallel    

  END SUBROUTINE get_potential_reg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION toproject(iii,nproj,cv_proj) RESULT (doproject)
    IMPLICIT NONE
    ! arguments
    INTEGER, INTENT(IN) :: iii, nproj
    INTEGER, INTENT(IN) :: cv_proj(nproj)
    LOGICAL :: doproject
    !local variables
    INTEGER :: i
    !
    doproject=.false.
    DO i=1,nproj
       IF (cv_proj(i) == iii) doproject=.true.
    END DO
  END FUNCTION toproject
  
END MODULE pot_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Vreco

USE kinds
USE pot_data 

IMPLICIT NONE

REAL (KIND=dp) :: v, sp, la, vmax, maxw, minw, ed(5)
REAL (KIND=dp) :: cube_tmp(6),cube_rad,norm_proj,v_proj
INTEGER        :: i, j, k, l, iii, ios, il, ix, iy, iz, kk, ll
INTEGER        :: ncoll, dtype, ndyn, eunit, nproj, mtdstp, lproj, iproj, ntorsion
LOGICAL        :: prnt_dyn, walker, read_int, cube_cv, cube_only
LOGICAL        :: proj, cube_cvmdck, cube_colvar, cube_full, reduce
LOGICAL, ALLOCATABLE :: is_torsion(:)

REAL (KIND=dp), ALLOCATABLE :: tmp(:), smin(:), smax(:), dg(:),         &
     s(:), s1(:), s2(:), dif1(:), dif2(:), so(:,:), ds(:), w(:), g(:),  &
     scal(:,:), maxcol(:), mincol(:), scal0(:), cube_value(:,:,:)

INTEGER, ALLOCATABLE :: scoll(:), n(:),cv_proj(:),dims(:),cv_red(:), int_tmp(:)
INTEGER, EXTERNAL    :: OMP_GET_NUM_THREADS
INTEGER              :: nthreads

CHARACTER(LEN=80)    :: pfmt0, pfmt1, pfmt2, pfmt4

CHARACTER(LEN=8), PARAMETER :: eunit_string(3) = (/ 'a.u     ', 'kJ/mol  ',  'kcal/mol' /)
REAL (KIND=dp),   PARAMETER :: econv(3)        = (/  1.0_dp,    2625.500_dp,  627.5132632235_dp /)
REAL (KIND=dp),   PARAMETER :: zero =0.0_dp

WRITE(*,'(t2,a)')  '**************************************************************'
WRITE(*,'(t2,a)')  '**     Reconstructing Metadynamics Free Energy Surface for  **'
WRITE(*,'(t2,a)')  '**                        CPMD                              **'
WRITE(*,'(t2,a,a,a)')  '**                     Version ',version,'                         **'
WRITE(*,'(t2,a)')  '**************************************************************'
WRITE(*,'(t2,a)')  '**                                                          **'
WRITE(*,'(t2,a)')  '**     Required files :                                     **'
WRITE(*,'(t2,a)')  '**                      colvar_mtd                          **'
WRITE(*,'(t2,a)')  '**                      parvar_mtd                          **' 
WRITE(*,'(t2,a)')  '**                      enevar_mtd [optional]               **'
WRITE(*,'(t2,a)')  '**                      cvmdck_mtd [optional]               **'
WRITE(*,'(t2,a)')  '**     Output files   :                                     **'
WRITE(*,'(t2,a)')  '**                      V.final.out                         **'
WRITE(*,'(t2,a)')  '**                      V.dynamic.out [optional]            **'
WRITE(*,'(t2,a)')  '**                      walker.out    [optional]            **'
WRITE(*,'(t2,a)')  '**                      V.final.cube  [optional]            **'
WRITE(*,'(t2,a)')  '**                      V.colvar.cube [optional]            **'
WRITE(*,'(t2,a)')  '**                      V.cvmdck.cube [optional]            **'
WRITE(*,'(t2,a)')  '**                                                          **'
WRITE(*,'(t2,a)')  '**     Help           : see the README file                 **'
WRITE(*,'(t2,a)')  '**     Contact        : nisanth.nair@theochem.rub.de        **'
WRITE(*,'(t2,a,/)')'**************************************************************'

proj=.false.
reduce=.false.
nproj=0
shft=.false.
tube=.false.
walker=.false.
prnt_dyn=.false.
rcut=zero
dyn_frq=1
eunit=1
cube_cv=.false.
cube_only=.false.
dtype=1

cube_cvmdck=.false.
cube_full=.false.
cube_colvar=.false.

read_int=.false.
nd=0
ncoll=0
mtdstp=0

nd=0
READ(5,'(a)',err=333,end=333) pfmt0(1:80)
IF (INDEX(pfmt0,'NCV') /= 0) THEN
   IF (INDEX(pfmt0,'RECON') /= 0) THEN
      READ(5,*,err=333,end=333)ncoll,nd
      ALLOCATE(scoll(nd))
      DO i=1,nd
         READ(5,*,err=333,end=333)scoll(i)
      END DO
   ELSE 
      READ(5,*,err=333,end=333)ncoll
      nd=ncoll
      ALLOCATE(scoll(nd))
      DO i=1,nd
         scoll(i)=i
      END DO
   END IF
ELSE
   STOP '*** First Keyword should be NCV ***'
END IF

!allocating
ALLOCATE(smin(nd))   !Grid minimum
ALLOCATE(smax(nd))   !Grid maximum
ALLOCATE(dg(nd))     !Grid width
ALLOCATE(n(nd))      !Number of grid points in each dimensions
ALLOCATE(dif1(nd))   !Displacements of collecitve coordinates 
ALLOCATE(dif2(nd))   ! -do- 
ALLOCATE(s(nd))      !grid value for which potential is calculated over all mtd time steps
ALLOCATE(s1(nd))     !gaussian center at time t
ALLOCATE(s2(nd))     !gaussian center at time t+1
ALLOCATE(g(nd))      !grid value
ALLOCATE(maxcol(nd)) !max coll. variable
ALLOCATE(mincol(nd)) !min coll. variable
ALLOCATE(int_tmp(nd)) 
ALLOCATE(is_torsion(nd)) !min coll. variable
is_torsion(1:nd)=.false.

read_input: DO j=1,999 !input file
   IF (j == 999) stop '*** END is missing in the inputfile ***'
   !Reading grid data; Calculate total number of grid points
   READ(5,'(a)',err=333,end=333) pfmt0(1:80) 
   IF (INDEX(pfmt0,'END') /= 0) THEN
      EXIT read_input
   ELSE IF(INDEX(pfmt0,'GRIDS') /= 0) THEN
      DO i=1,nd
         READ(5,*,err=333,end=333) smin(i),smax(i),dg(i)
      END DO
   ELSE IF (INDEX(pfmt0,'HILL') /= 0) THEN
      IF (INDEX(pfmt0,'SPHER') /= 0) THEN
         tube=.false. ! use spherical Gaussians
      ELSE IF (INDEX(pfmt0,'TUBE') /= 0) THEN
         tube=.true.  ! use tubular  Gaussians 
      END IF
      IF (INDEX(pfmt0,'RCUT') /= 0) THEN
         shft=.true.
         READ(5,*,err=333,end=333) rcut
      END IF
   ELSE IF (INDEX(pfmt0,'PRINT') /= 0) THEN
      IF (INDEX(pfmt0,'DYN') /= 0) THEN
         prnt_dyn=.true. 
         IF (INDEX(pfmt0,'WAL') /= 0) walker=.true.
         READ(5,*,err=333,end=333) dyn_frq
      END IF
   ELSE IF (INDEX(pfmt0,'OUT') /= 0) THEN
      IF (INDEX(pfmt0,'NOGRID') /= 0) THEN
         dtype=2
      ELSE
         dtype=1
      END IF
   ELSE IF (INDEX(pfmt0,'UNIT') /= 0) THEN
      IF (INDEX(pfmt0,'KJ') /= 0) THEN
         eunit=2
      ELSE IF (INDEX(pfmt0,'KC') /= 0) THEN
         eunit=3
      ELSE 
         eunit=1
      END IF
   ELSE IF (INDEX(pfmt0,'METASTEP') /= 0) THEN
      READ(5,*,err=333,end=333)mtdstp
   ELSE IF (INDEX(pfmt0,'CUBE') /= 0) THEN
      cube_cv=.true.
      IF (INDEX(pfmt0,'ONLY')   /= 0) cube_only=.true.
      IF (INDEX(pfmt0,'RADIUS') /= 0) READ(*,*,err=333,end=333) cube_rad
      IF (INDEX(pfmt0,'CVMDCK') /= 0) cube_cvmdck=.true.
      IF (INDEX(pfmt0,'COLVAR') /= 0) cube_colvar=.true.
      IF (INDEX(pfmt0,'FULL')   /= 0) cube_full=.true.
      IF (.not.(cube_cvmdck.or.cube_colvar.or.cube_full))cube_full=.true.
   ELSE IF (INDEX(pfmt0,'PROJECT') /= 0) THEN
      proj=.true.
      IF(reduce) STOP ' PROJECT and REDUCE can not be used together! '
      READ(*,*,err=333,end=333)nproj
      ALLOCATE(cv_proj(nproj))
      DO i=1,nproj
         READ(*,*,err=333,end=333)cv_proj(i)
      END DO
   ELSE IF (INDEX(pfmt0,'REDUCE') /= 0) THEN
      reduce=.true.
      IF(proj) STOP ' PROJECT and REDUCE can not be used together! '
      READ(*,*,err=333,end=333)nproj
      ALLOCATE(cv_proj(nproj))
      DO i=1,nproj
         READ(*,*,err=333,end=333)cv_proj(i)
      END DO
   ELSE IF (INDEX(pfmt0,'TORSION') /= 0) THEN
        READ(*,*,err=333,end=333)ntorsion
        READ(*,*,err=333,end=333)int_tmp(1:ntorsion)
   ELSE
      PRINT *, 'unknown keyword: ', pfmt0
   END IF
END DO read_input

IF ((.not.tube).and.shft) THEN
  WRITE(*,'(A)') '!!!!WARNING!!!! SHIFT is only possible with TUBE in CPMD'
  WRITE(*,'(A)') '!!!!WARNING!!!! check your input! results can be wrong!'
ENDIF
IF (proj) THEN
  WRITE(*,'(A)') '!!!!WARNING!!!! NOT FULLY TESTED'
  WRITE(*,'(A)') '!!!!WARNING!!!! NOT FULLY TESTED'
ENDIF


!rearrange according to projection
ALLOCATE(dims(ncoll))
IF (proj.or.reduce) THEN
   kk=1
   j=1
   DO iii=1,ncoll
      IF (toproject(iii,nproj,cv_proj)) THEN
         dims(iii)=ncoll-nproj+kk
         kk=kk+1
      ELSE
         dims(iii)=j
         j=j+1
      END IF
   END DO
   write(*,'(t2,a)')'!WARNING! Dimensions are rearranged as below:'
   DO iii=1,ncoll
      WRITE(*,'(t2,i8,a,i8)')iii,'     ->',dims(iii)
   END DO
ELSE 
   DO i=1,ncoll
      dims(i)=i
   END DO
ENDIF

if(ntorsion.gt.0)then
  do i=1,ntorsion
    iii=int_tmp(i)
    is_torsion(dims(iii))=.true.
  end do
end if
DEALLOCATE(int_tmp)

IF (proj.or.reduce) THEN
   ALLOCATE(tmp(3*nd))
   DO iii=1,nd
      tmp(dims(iii))=smin(iii)
      tmp(dims(iii)+nd)=smax(iii)
      tmp(dims(iii)+2*nd)=dg(iii)
   END DO
   DO i=1,nd
      smin(i)=tmp(i)
      smax(i)=tmp(i+nd)
      dg(i)  =tmp(i+2*nd)  
   END DO
   DEALLOCATE(tmp)
ENDIF
!Calculate total number of grid points
l=1
DO i=1,nd
  n(i) = NINT((ABS(smax(i)-smin(i)))/dg(i)) + 1 
  IF((reduce).AND.(i.GT.(nd-nproj)))then
    n(i)=1
    smax(i)=smin(i)
    WRITE(*,'(t2,a,i2,a,f8.4)')  'CV #',i,' is reduced at value  =', smin(i)
  END IF
  l = l*n(i)
END DO


OPEN(11,file='colvar_mtd',status='old',form='formatted',iostat=ios)
IF (ios /= 0) STOP '!!! Cannot open file "colvar_mtd" for reading'
OPEN(12,file='parvar_mtd',status='old',form='formatted',iostat=ios)
IF (ios /= 0) STOP '!!! Cannot open file "parvar_mtd" for reading'
OPEN(13,file='V.final.out',  status='unknown',form='formatted',iostat=ios)
IF (ios /= 0) STOP '!!! Cannot open file "V.final.out" for writing'

IF (prnt_dyn) OPEN(14,file='V.dynamic.out',status='unknown',form='formatted')

IF (walker) then
   OPEN(17,file='enevar_mtd',status='old',form='formatted',iostat=ios)
   IF (ios /= 0) STOP '!!! PRINT WALKER requires the file "enevar_mtd"'
   OPEN(15,file='walker.out',status='unknown',form='formatted')
ENDIF

WRITE(*,*)
!
! OpenMP info and setup.
!$OMP parallel
!$ nthreads = OMP_GET_NUM_THREADS()
!$ CALL OMP_SET_NESTED(.false.)
!$OMP master
!$ WRITE(*,'(t2,a,i4,a,/)')  'OpenMP version using ', nthreads, ' threads.'
!$OMP end master
!$OMP end parallel
!
WRITE(*,'(t2,a,i8)')  'Number of coll. coordinates     =', ncoll
IF (nd /= ncoll)&
WRITE(*,'(t2,a,i8)')  'Number of CVs to reconstruct    =', nd
IF (proj)&
WRITE(*,'(t2,a,i8)')  'Number of CVs to project        =', nd
IF (reduce)&
WRITE(*,'(t2,a,i8)')  'Number of CVs to reduce         =', nd

pfmt4='*'
IF (nd<10) WRITE(pfmt4,'(a6,i1,a3)')'(t2,a,',nd,'i8)'
IF (nd /= ncoll) WRITE(*,pfmt4) 'CVs selected                    =',scoll(1:nd)
IF (proj) THEN
   pfmt4='*'
   IF (nd<10) WRITE(pfmt4,'(a6,i1,a3)')'(t2,a,',nproj,'i8)'
   WRITE(*,pfmt4) 'CVs to project                  =',cv_proj(1:nproj)
END IF
IF (nd==1)pfmt4='(t2,a,i8)'
IF ((nd>1).and.(nd-1<10)) WRITE(pfmt4,'(a9,i1,a12)')'(t2,a,i8,',nd-1,'("   x",i8))'
WRITE(*,pfmt4)    'Grid Size                       =', n(1:nd)    
WRITE(*,pfmt4)    'Torsions:'
do it=1,nd
  if(is_torsion(it))then 
   write(*,'(a,i5,a)') '                          ', it, ': True '
  else
   write(*,'(a,i5,a)') '                          ', it, ': False '
  end if
end do

!find # of mtd steps
it=0
ALLOCATE(tmp(ncoll))
DO 
  READ(11,*,iostat=ios)tmp(1:ncoll), la
  IF (ios/=0) EXIT
  it=it+1
END DO
DEALLOCATE(tmp)

WRITE(*,'(t2,a,i8)') 'Number of Gaussians found       =', it
IF (mtdstp > it) THEN
   WRITE(*,'(t2,a,i8,a,i8,a)') '!!Warning: METASTEP=',mtdstp,' > number of Gaussians=', it,'!!'
   WRITE(*,'(t2,a,i8)') '!!Warning: METASTEP is set to ',it
ELSE IF (mtdstp.gt.0) THEN
   it=mtdstp
ENDIF
WRITE(*,'(t2,a,i8)') 'Number of Gaussians used        =', it
REWIND(11) !rewinding colvar_mtd file

!Get the printing formats for V.dynamic.out
pfmt1=' ' 
ndyn=1
IF (prnt_dyn) THEN
   ndyn = (it+dyn_frq-1)/dyn_frq
   IF (dtype==1) i = ndyn+nd-nproj
   IF (dtype==2) i = ndyn

   IF (i<10) THEN
      WRITE(pfmt1,'(a1,i1,a6)')'(',i,'f16.8)'
   ELSE IF (i<100) THEN
      WRITE(pfmt1,'(a1,i2,a6)')'(',i,'f16.8)'
   ELSE IF (i<1000) THEN
      WRITE(pfmt1,'(a1,i3,a6)')'(',i,'f16.8)'
   ELSE IF (i<10000) THEN
      WRITE(pfmt1,'(a1,i4,a6)')'(',i,'f16.8)'
   ELSE IF (i<100000) THEN
      WRITE(pfmt1,'(a1,i5,a6)')'(',i,'f16.8)'
   END IF

END IF

pfmt0='*'
IF (dtype==1) THEN
   i=nd+1
   IF (proj) i = nd-nproj+1
   IF (i<10) WRITE(pfmt0,'(a,i1,a6)')'(',i,'f16.8)'
ENDIF
IF (dtype==2) pfmt0='(f16.8)'

IF (walker) THEN
   pfmt2='*'
   i=nd+1
   IF (proj) i = nd-nproj
   IF (i<10) WRITE(pfmt2,'(a,i1,a6)')'(',i,'f16.8)'
END IF

ALLOCATE(so(it,ncoll))
ALLOCATE(ds(it))
ALLOCATE(w(it))
ALLOCATE(scal0(ncoll))
ALLOCATE(scal(it,nd))

!Read and store the information from colvar_mtd and parvar_mtd and others
DO i=1,it
   READ(11,*,iostat=ios)j, so(i,1:ncoll), scal0(1:ncoll)
   IF (walker) THEN
      READ(17,*) j,ed(1:5)
      WRITE(15,pfmt2) (so(i,scoll(iii)),iii=1,nd),ed(4)+ed(5)
   END IF
   DO iii=1,nd
      IF (i==1) THEN
         mincol(iii)=so(i,scoll(iii))
         maxcol(iii)=so(i,scoll(iii))
      ELSE
         mincol(iii)=MIN(mincol(iii),so(i,scoll(iii)))
         maxcol(iii)=MAX(maxcol(iii),so(i,scoll(iii)))
      END IF
      scal(i,iii)=scal0(scoll(iii))
   END DO
   IF (ios<0) THEN
      WRITE(*,*) 'premature end of file colvar_mtd'
      WRITE(*,*) 'i =', i
      STOP
   ELSE IF (ios >0) THEN
      WRITE(*,*) 'error reading file colvar_mtd'
      WRITE(*,*) 'i =', i
      STOP
   END IF

   READ(12,*,iostat=ios) j, sp, ds(i), w(i)
   IF (i==1) THEN
      maxw=w(i)
      minw=w(i)
   ELSE
      maxw=MAX(maxw,w(i))
      minw=MIN(minw,w(i))
   END IF

   IF (ios<0) THEN
      WRITE(*,*) 'premature end of file parvar_mtd'
      WRITE(*,*) 'i =', i
      STOP
   ELSEIF (ios >0) THEN
      WRITE(*,*) 'error reading file parvar_mtd'
      WRITE(*,*) 'i =', i
      STOP
   END IF
END DO

ALLOCATE(tmp(nd*2))
DO j=1,it
   DO iii=1,nd  
      tmp(dims(iii))=so(j,iii) 
      tmp(dims(iii)+nd)=scal(j,iii)
   END DO
   so(j,1:nd)=tmp(1:nd)
   scal(j,1:nd)=tmp(ncoll+1:2*ncoll)
END DO
DEALLOCATE(tmp)

!print some more details
DO iii=1,nd
  WRITE(*,'(t2,a,i4,a,f8.4)') 'Min Value of Coll. Cord. #',scoll(iii), &
                              '  =', mincol(iii)
  WRITE(*,'(t2,a,i4,a,f8.4)') 'Max Value of Coll. Cord. #',scoll(iii), &
                              '  =', maxcol(iii)
END DO
WRITE(*,'(t2,a,f8.4,1x,a)') 'Min Value of Hill Height        =',minw*econv(eunit),trim(eunit_string(eunit))
WRITE(*,'(t2,a,f8.4,1x,a)') 'Max Value of Hill Height        =',maxw*econv(eunit),trim(eunit_string(eunit))


WRITE(*,*)
WRITE(*,'(t2,a)'   )    'Total potential written to      = V.final.out'
IF (prnt_dyn) WRITE(*,'(t2,a,i6,a)')  &
     'Dynamic potential written to    = V.dynamic.out  every',dyn_frq,' steps' 
IF (walker)   WRITE(*,'(t2,a,i6,a)')  &
     'Walker positions  written to    = walker.out     every',dyn_frq,' steps'

!find the normalization factor for projection
IF (proj) THEN
   norm_proj=1.0_dp
   lproj=1
   DO i=1,nproj
      norm_proj=norm_proj*DBLE(n(dims(cv_proj(i))))*dg(dims(cv_proj(i)))
      lproj=lproj*n(dims(cv_proj(i)))
   END DO
   DO i=1,nproj
      norm_proj=norm_proj/dg(dims(cv_proj(i)))
   END DO
   IF (norm_proj<1.0d-10) STOP '** Error: projection normalization too small **'
END IF

g(1:nd)=smin(1:nd)
IF (cube_cv) THEN
   IF ((nd-nproj)/=3) STOP '*** Generating Cube file NOT possible..because #CV /= 3 ***'
300 CONTINUE
   IF (cube_cvmdck) THEN
      WRITE(*,'(t2,a)') 'Final potential for cvmdck_mtd  as volumetric data in V.cvmdck.cube'
      OPEN(16,file='V.cvmdck.cube',status='unknown',form='formatted')
   ELSE IF (cube_full) THEN
      WRITE(*,'(t2,a)') 'Final potential for entire space as volumetric data in V.final.cube'
      OPEN(16,file='V.final.cube',status='unknown',form='formatted')
   ELSE IF (cube_colvar) THEN
      WRITE(*,'(t2,a)') 'Final potential for colvar_mtd   as volumetric data in V.colvar.cube'
      OPEN(16,file='V.colvar.cube',status='unknown',form='formatted')
   END IF

   WRITE(16,'(a)') 'CV cube file'
   WRITE(16,'(a,a)') 'generated by Vreco_CPMD.x ', version
   WRITE(16,'(i5,3f12.6)') 1, smin(1), smin(2), smin(3)  
   WRITE(16,'(i5,3f12.6)') n(1), dg(1), zero,  zero     ! grid x
   WRITE(16,'(i5,3f12.6)') n(2), zero,  dg(2), zero     ! grid y
   WRITE(16,'(i5,3f12.6)') n(3), zero,  zero,  dg(3)    ! grid z
   WRITE(16,'(i5,4f12.6)') 1, zero, zero,  zero,  zero     ! "fake" atom coordinates (one h-atom at (0,0,0)
   
   ALLOCATE(cube_value(n(1),n(2),n(3)))
   ALLOCATE(tmp(3*nd))
   cube_value = zero
   il=0
   IF (cube_cvmdck) OPEN(18,file='cvmdck_mtd',status='old',form='formatted')
   DO_cube: DO il=1,l
      IF (cube_cvmdck) THEN
         read(18,*,iostat=ios)i,j,tmp(1:2*nd) 
         IF (ios /= 0) THEN
            cube_cvmdck=.false.
            exit DO_cube
            close(18)
         END IF
         s(1:nd)=tmp(nd+1:2*nd+1)
         ix=NINT((s(1)-smin(1))/dg(1))+1
         iy=NINT((s(2)-smin(2))/dg(2))+1
         iz=NINT((s(3)-smin(3))/dg(3))+1
         CALL get_potential(so,scoll,is_torsion,scal,ds,s,w,.false.,1,tmp,v)
         cube_value(ix,iy,iz)=MAX(cube_value(ix,iy,iz),v)
      ELSE IF (cube_full) THEN
         IF (il == 1)s(1:3)=smin(1:3)
         IF (il == 1)g(1:3)=smin(1:3)
         IF (il == l) THEN
            cube_full=.false.
         END IF
         ix=NINT((s(1)-smin(1))/dg(1))+1
         iy=NINT((s(2)-smin(2))/dg(2))+1
         iz=NINT((s(3)-smin(3))/dg(3))+1
         CALL get_potential(so,scoll,is_torsion,scal,ds,s,w,.false.,1,tmp,v)
         cube_value(ix,iy,iz)=MAX(MAX(cube_value(ix,iy,iz),v),1.0d-6)
         g(nd)=g(nd)+dg(nd)
         ll=1
         DO_dim0:DO k=nd,1,-1 
            ll=ll*n(k)
            IF (MOD(il,ll) == 0) THEN
               g(k-1)=g(k-1)+dg(k-1)
               DO kk=k,nd
                  g(kk)=smin(kk)
               END DO
            ELSE
               EXIT DO_dim0
            END IF
         END DO DO_dim0
         s(1:nd)=g(1:nd)
      ELSE IF (cube_colvar) THEN
         IF (il == it) THEN
            cube_colvar=.false.
         END IF
         s(1:3)=so(il,scoll(1:3))
         ix=NINT((s(1)-smin(1))/dg(1))+1
         iy=NINT((s(2)-smin(2))/dg(2))+1
         iz=NINT((s(3)-smin(3))/dg(3))+1
         CALL get_potential(so,scoll,is_torsion,scal,ds,s,w,.false.,1,tmp,v)
         cube_value(ix,iy,iz)=MAX(cube_value(ix,iy,iz),v)
      END IF
   END DO DO_cube
   DEALLOCATE(tmp)
!   il=0
!   ix=1
!   iy=1
!   iz=1
!   k=0
!   DO il=1,l
!      IF ((il>1).and.(MOD(il-1,n(3))) == 0) THEN
!         iz=1
!         iy=iy+1
!      END IF
!      IF ((il>1).and.(MOD(il-1,n(2)*n(3))) == 0) THEN
!         iz=1
!         iy=1
!         ix=ix+1
!      END IF
!      k=k+1
!      cube_tmp(k)=cube_value(ix,iy,iz)
!!      iz=iz+1
!      IF(iz==n(3)) THEN
!        WRITE(pfmt4,'(a,i1,a6)')'(',k,'e13.5)'
!        WRITE(16,pfmt4) -cube_tmp(1:k)*econv(eunit)
!      ELSE IF(MOD(iz,6)==0) THEN
!        WRITE(16,'(6e13.5)') -cube_tmp(1:6)*econv(eunit)
!        k=0
!      END IF
!      iz=iz+1

!Nn check
!do ix=1,n(1)
  !do iy=1,n(2)
    !do k=1,n(3)/6
      !WRITE(16,'(6e13.5)')(-cube_value(ix,iy,iz)*econv(eunit), iz=(k-1)*6+1,k*6)
    !end do
   !end do
!end do
do ix=1,n(1)
  do iy=1,n(2)
      WRITE(16,'(6e13.5)')(-cube_value(ix,iy,iz)*econv(eunit), iz=1,n(3))
   end do
end do

!      IF (il == l) THEN
!        WRITE(pfmt4,'(a,i1,a6)')'(',k,'e13.5)'
!        WRITE(16,pfmt4) -cube_tmp(1:k)*econv(eunit)
!     ELSE IF (MOD(il,6) == 0) THEN
!        WRITE(16,'(6e13.5)') -cube_tmp(1:6)*econv(eunit)
!        k=0
!     END IF

!   END DO
   
   CLOSE(16)
   DEALLOCATE(cube_value)
   IF (cube_colvar.or.cube_full.or.cube_cvmdck) THEN
      j=1
      g(1:nd)=smin(1:nd)
      GOTO 300
   END IF
   WRITE(*,'(/,t2,a,/)')  '================    Cube File Completed    =================== '
END IF

IF (cube_only) STOP

ALLOCATE(tmp(1:ndyn))

vmax    = zero
v_proj  = zero
iproj   = 0
g(1:nd) = smin(1:nd)

grids: DO il=1,l ! loop over all grid points
   
   s(1:nd) = g(1:nd)

   CALL get_potential(so,scoll,is_torsion,scal,ds,s,w,prnt_dyn,ndyn,tmp,v)
   
   IF (proj) THEN

      v_proj = v_proj+v
      IF (MOD(il,lproj) == 0) THEN
         iproj = iproj+1 
         v_proj = v_proj/norm_proj !normalize
         vmax = MAX(vmax,v_proj)
         IF (dtype == 1) THEN
            !TODO
            !to correct pfmt
            WRITE(13,pfmt0) s(1:nd-nproj), v_proj*(-1.d0)*econv(eunit)
            !TODO
            !to correct tmp
            !TODO
            !to correct pfmt
            IF (prnt_dyn) WRITE(14,pfmt1) s(1:nd-nproj), tmp(1:ndyn)*(-1.d0)*econv(eunit)
         ELSE 
            WRITE(13,pfmt0) v_proj*(-1.0)*econv(eunit)
            IF (prnt_dyn) WRITE(14,pfmt1) tmp(1:ndyn)*(-1.d0)*econv(eunit)
         END IF
         v_proj = zero
      END IF

      IF (prnt_dyn.and.(MOD(il,lproj*n(nd-nproj))==0)) WRITE(14,*)' '
      IF (MOD(il,lproj*n(nd-nproj))==0) WRITE(13,*)' '

   ELSE

      vmax = MAX(vmax,v)
      IF (dtype == 1) THEN
         WRITE(13,pfmt0) s(1:nd-nproj), v*(-1.d0)*econv(eunit)
         IF (prnt_dyn) WRITE(14,pfmt1) s(1:nd-nproj), tmp(1:ndyn)*(-1.d0)*econv(eunit)
      ELSE 
         WRITE(13,pfmt0) v*(-1.0)*econv(eunit)
         IF (prnt_dyn) WRITE(14,pfmt1) tmp(1:ndyn)*(-1.d0)*econv(eunit)
      END IF

      IF (prnt_dyn.and.(MOD(il,n(nd-nproj))==0)) WRITE(14,*)' '
      IF (MOD(il,n(nd-nproj))==0) WRITE(13,*)' '

   END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To get the grid point of N dimension
   g(nd) = g(nd)+dg(nd)
   ll = 1
   DO_dim: DO k=nd,2,-1
      ll = ll*n(k)
      IF (MOD(il,ll) == 0) THEN
         g(k-1) = g(k-1)+dg(k-1)
         DO kk=k,nd
            g(kk) = smin(kk)
         END DO
      ELSE
         EXIT DO_dim
      END IF
   END DO DO_dim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END DO grids

IF (ALLOCATED(tmp)) DEALLOCATE(tmp)

WRITE(*,'(t2,a,f16.8,1x,a)')'Max Reconstructed Potential     =', vmax*econv(eunit), &
                                                             trim(eunit_string(eunit))
WRITE(*,'(/,t2,a,/)')  '================ Reconstruction  Completed =================== '
WRITE(*,*)

CLOSE(11)
CLOSE(12)
CLOSE(13)
IF (prnt_dyn) CLOSE(14)
IF (walker) THEN
  CLOSE(17)
  CLOSE(15)
END IF

DEALLOCATE(scoll)
DEALLOCATE(smin)
DEALLOCATE(smax)
DEALLOCATE(dg)
DEALLOCATE(n)
DEALLOCATE(dif1)
DEALLOCATE(dif2)
DEALLOCATE(s)
DEALLOCATE(s1)
DEALLOCATE(s2)
DEALLOCATE(so)
DEALLOCATE(w)
DEALLOCATE(ds)
DEALLOCATE(g)
DEALLOCATE(scal)
DEALLOCATE(scal0)
DEALLOCATE(maxcol)
DEALLOCATE(mincol)
DEALLOCATE(dims)
IF (ALLOCATED(cv_proj)) DEALLOCATE(cv_proj)
STOP

!error
333 CONTINUE
WRITE(*,'(/,a,/)')'xxxxxxxx ERROR IN READING INPUT xxxxxxxx'
WRITE(*,'(/,a,a,/)')'Last line read is:',pfmt0
STOP
!
END PROGRAM Vreco
