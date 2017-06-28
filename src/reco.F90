!****************************************************************************************!
!
!****************************************************************************************!
PROGRAM RecoPar
IMPLICIT NONE

REAL*8, ALLOCATABLE :: s0(:,:),w(:),ds(:,:),grid0(:,:), v(:)
INTEGER :: ncv, nmtd

REAL*8 :: Vreco,temp, deltat, alpha
LOGICAL :: parent
INTEGER :: i,j, ngrid1,ngrid2,ngrid3

CALL MPI_Start
CALL Set_Parent(parent)

IF(parent)THEN
   READ(*,*)ncv, nmtd, temp, deltat
  alpha=(temp+deltat)/deltat
END IF

CALL IBcast(ncv,1)
CALL IBcast(nmtd,1)
CALL RBcast(alpha,1)

ALLOCATE(grid0(3,ncv))
ALLOCATE(s0(ncv,nmtd))
ALLOCATE(w(nmtd))
ALLOCATE(ds(ncv,nmtd))

IF(parent)THEN
  DO i=1,ncv
    READ(*,*)grid0(1:3,i)
  END DO
END IF
CALL RBcast(grid0,3*ncv)


IF(parent)THEN
  OPEN(1,file='HILLS',STATUS='OLD')
  DO i=1,nmtd
    READ(1,*)j,s0(1:ncv,i),w(i),ds(1:ncv,i)
  END DO 
  CLOSE(1)
END IF
CALL RBcast(s0,ncv*nmtd)
CALL RBcast(w,nmtd)
CALL RBcast(ds,ncv*nmtd)


IF(ncv.EQ.1)THEN
  ngrid1=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
  ALLOCATE(v(ngrid1))
  v(1:ngrid1)=0.0
  CALL sumhills1d(ncv,nmtd,grid0,s0,w,ds,temp,deltat)
  ngrid2=1
  ngrid3=1
ELSE IF(ncv.EQ.2)THEN
  ngrid1=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
  ngrid2=nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1
  ALLOCATE(v(ngrid1*ngrid2))
  v(1:ngrid1*ngrid2)=0.0
  CALL sumhills2d(ncv,nmtd,grid0,s0,w,ds)
  ngrid3=1
ELSE IF(ncv.EQ.3)THEN
  ngrid1=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
  ngrid2=nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1
  ngrid3=nint((grid0(2,3)-grid0(1,3))/grid0(3,3))+1
  ALLOCATE(v(ngrid1*ngrid2*ngrid3))
  v(1:ngrid1*ngrid2*ngrid3)=0.0
  CALL sumhills3d(ncv,nmtd,grid0,s0,w,ds)
ELSE
  PRINT *, 'ncv =', ncv, ' not implemented'
  STOP
END IF

CALL Print_Potential(parent,ncv,ngrid1*ngrid2*ngrid3,grid0,v)

CALL MPI_Stop
END PROGRAM RecoPar
!****************************************************************************************!


FUNCTION Vmeta(ncv,nmtd,w,ds,s0,s)
!Computes Metadynamics potential
IMPLICIT NONE
REAL*8 :: Vmeta
INTEGER :: ncv,nmtd
REAL*8 :: w(*), ds(ncv,*), s0(ncv,*), s(*)
!
REAL*8 :: dd(ncv),gx
INTEGER :: imtd,i

Vmeta=0.D0
DO imtd=1,nmtd
  DO i=1,ncv
    dd(i)=(s(i)-s0(i,imtd))/ds(i,imtd)
  END DO
  gx=DOT_PRODUCT(DD(1:3),DD(1:3))/2.D0
  Vmeta=Vmeta+w(imtd)*DEXP(gx)
END DO

END 
!****************************************************************************************!

SUBROUTINE sumhills1d(ncv,nmtd,grid0,s0,w,ds,alpha,v)
!Sum Hills 1D
IMPLICIT NONE
INTEGER :: ncv, nmtd
REAL*8 :: grid0(3,ncv), s0(ncv,nmtd), w(nmtd), ds(ncv,nmtd), alpha,v(*)
!
REAL*8 :: Vmeta,grid(3,1), s(ncv)
INTEGER :: ig1, gleng1, rank, ii
!
CALL DistributeGrids(ncv,grid0,grid,rank)
gleng1=nint((grid(2,1)-grid(1,1))/grid(3,1))+1
ii=rank
DO ig1=1,gleng1
 s(1)=grid(1,1)+DFLOAT(ig1-1)*grid(3,1)
 v(ii)=Vmeta(ncv,nmtd,w,ds,s0,s)
 ii=ii+1
END DO
!
END
!****************************************************************************************!


SUBROUTINE sumhills2d(ncv,nmtd,grid0,s0,w,ds,alpha,v)
!Sum Hills 2D
IMPLICIT NONE
INTEGER :: ncv, nmtd
REAL*8 :: grid0(3,ncv), s0(ncv,nmtd), w(nmtd), ds(ncv,nmtd), alpha,v(*)
!
REAL*8 :: Vmeta,grid(3,2), s(ncv)
INTEGER :: ig1, ig2, gleng1, gleng2, rank, ii
!
CALL DistributeGrids(ncv,grid0,grid,rank)
gleng1=nint((grid(2,1)-grid(1,1))/grid(3,1))+1
gleng2=nint((grid(2,2)-grid(1,2))/grid(3,2))+1
ii=rank
DO ig1=1,gleng1
  s(1)=grid(1,1)+DFLOAT(ig1-1)*grid(3,1)
  DO ig2=1,gleng2
    s(2)=grid(1,2)+DFLOAT(ig2-1)*grid(3,2)
    v(ii)=Vmeta(ncv,nmtd,w,ds,s0,s)
    ii=ii+1
  END DO
END DO
!
END
!****************************************************************************************!

SUBROUTINE sumhills3d(ncv,nmtd,grid0,s0,w,ds,alpha,v)
!Sum Hills 3D
IMPLICIT NONE
INTEGER :: ncv, nmtd
REAL*8 :: grid0(3,ncv), s0(ncv,nmtd), w(nmtd), ds(ncv,nmtd), alpha,v(*)
!
REAL*8 :: Vmeta,grid(3,3), s(ncv)
INTEGER :: ig1, ig2, ig3, gleng1, gleng2, gleng3, rank, ii
!
CALL DistributeGrids(ncv,grid0,grid,rank)
gleng1=nint((grid(2,1)-grid(1,1))/grid(3,1))+1
gleng2=nint((grid(2,2)-grid(1,2))/grid(3,2))+1
gleng3=nint((grid(2,3)-grid(1,3))/grid(3,3))+1
ii=rank
DO ig1=1,gleng1
  s(1)=grid(1,1)+DFLOAT(ig1-1)*grid(3,1)
  DO ig2=1,gleng2
    s(2)=grid(1,2)+DFLOAT(ig2-1)*grid(3,2)
    DO ig3=1,gleng3
      s(3)=grid(1,3)+DFLOAT(ig3-1)*grid(3,3)
      v(ii)=Vmeta(ncv,nmtd,w,ds,s0,s)
      ii=ii+1
    END DO
  END DO 
END DO
!
END
!****************************************************************************************!

SUBROUTINE DistributeGrids(ncv,grid0,grid,rank)
!Distribute X grid over processors by mapping grid0 to grid 
IMPLICIT NONE
INTEGER :: ncv,rank
REAL*8 :: grid0(3,ncv), grid(3,ncv)
!
INTEGER :: i,ncpu,icpu,ngrids,ngrids_m

CALL Get_ncpu(ncpu)
CALL Get_cpuid(icpu)

DO i=1,ncv
  grid(1:3,i)=grid0(1:3,i)
END DO 
rank=0

!Distribute X grids
WRITE(*,'(2A7,4A16)') 'CPU','CV', 'GRID LEN.', 'GRID MIN', 'GRID MAX', 'GRID BIN'
IF(ncpu.GT.1)THEN
  ngrids=nint(grid0(2,1)-grid0(1,1)/grid0(3,1))
  ngrids=ngrids/ncpu
  IF(icpu.eq.ncpu-1)THEN
    ngrids_m=ngrids+mod(ngrids+1,ncpu)
    grid(1,1)=DFLOAT(icpu*ngrids)*grid0(1,1)
    grid(2,1)=DFLOAT(icpu*ngrids+ngrids_m-1)*grid0(1,1)
  ELSE
    ngrids_m=ngrids
    grid(1,1)=DFLOAT(icpu*ngrids)*grid0(1,1)
    grid(2,1)=DFLOAT(icpu*ngrids+ngrids-1)*grid0(1,1)
  END IF
  CALL Sync_procs
  DO i=1,ncv
    WRITE(*,'(2I7,I16,2F16.6)') icpu, ngrids_m, grid(1,i), grid(2,i), grid(3,i)
    CALL Sync_procs
  END DO
  rank=ngrids*icpu
END IF 
END 
!****************************************************************************************!

SUBROUTINE Print_Potential(parent,ncv,ngrid,grid,v)
!Print Potential
IMPLICIT NONE
LOGICAL :: parent
INTEGER :: ncv, ngrid
REAL*8 :: v(*), grid(3,ncv)
INTEGER :: gleng1, gleng2, gleng3, ig1, ig2, ig3, ii
REAL*8  :: s(3)

CALL GlobSumR(v,ngrid)

IF(.NOT.parent)RETURN

OPEN(11,FILE='V_sum.dat',FORM='FORMATTED')
IF(ncv.EQ.1)THEN
   gleng1=NINT((grid(2,1)-grid(1,1))/grid(3,1))+1
   DO ig1=1,gleng1
     s(1)=grid(1,1)+DFLOAT(ig1-1)*grid(3,1)
     WRITE(11,'(F16.6,E16.8)')s(1),v(ig1)
   END DO
ELSE IF(ncv.EQ.2)THEN
   gleng1=nint((grid(2,1)-grid(1,1))/grid(3,1))+1
   gleng2=nint((grid(2,2)-grid(1,2))/grid(3,2))+1
   ii=1
   DO ig1=1,gleng1
     s(1)=grid(1,1)+DFLOAT(ig1-1)*grid(3,1)
     DO ig2=1,gleng2
       s(2)=grid(1,2)+DFLOAT(ig2-1)*grid(3,2)
       WRITE(11,'(2F16.6,E16.8)')s(1),s(2), v(ii)
       ii=ii+1
     END DO
     WRITE(11,*)
   END DO
ELSE
   STOP 'ncv 3 print not implemented '
   ! CALL Write_cube(grid0,ncv) !FIXME
END IF
CLOSE(11)

END 
!****************************************************************************************!

SUBROUTINE MPI_Start()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
#if defined (_PARALLEL)
  call MPI_INIT(i_err)
#endif
END
!****************************************************************************************!

SUBROUTINE MPI_Stop()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
#if defined (_PARALLEL)
call MPI_FINALIZE(i_err)
#endif
END
!****************************************************************************************!

SUBROUTINE get_ncpu(ncpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: ncpu, i_err
ncpu=1
#if defined (_PARALLEL)
call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,i_err)
#endif
END
!****************************************************************************************!

SUBROUTINE get_cpuid(icpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: icpu, i_err
icpu=0
#if defined (_PARALLEL)
call MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
#endif
END
!****************************************************************************************!

SUBROUTINE IBcast(myint,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: leng, myint(*), i_err
#if defined (_PARALLEL)
CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
#endif
END
!****************************************************************************************!

SUBROUTINE RBcast(myreal,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: myreal(*)
INTEGER :: leng, i_err
#if defined (_PARALLEL)
CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
#endif
END
!****************************************************************************************!

SUBROUTINE Sync_Procs
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER i_err
#if defined (_PARALLEL)
call MPI_Barrier(MPI_COMM_WORLD,i_err)
#endif
END
!****************************************************************************************!

SUBROUTINE Set_Parent(parent)
IMPLICIT NONE
INCLUDE 'mpif.h'
LOGICAL :: parent
INTEGER :: icpu, i_err
parent=.false.
icpu=0
#if defined (_PARALLEL)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
#endif
IF(icpu.eq.0)parent=.true.
END
!****************************************************************************************!

SUBROUTINE GlobSumR(myreal,leng)
IMPLICIT NONE
REAL*8 :: myreal(*)
INTEGER :: leng,i_err
INCLUDE 'mpif.h'
#if defined (_PARALLEL)
call MPI_Allreduce(myreal,myreal,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
#endif
END
!****************************************************************************************!

SUBROUTINE GlobSumI(myint,leng)
IMPLICIT NONE
INTEGER :: myint(*)
INTEGER :: leng,i_err
INCLUDE 'mpif.h'
#if defined (_PARALLEL)
call MPI_Allreduce(myint,myint,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
#endif
END
!****************************************************************************************!
