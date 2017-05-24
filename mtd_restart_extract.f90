program mtdrest_extract
!Reads MTD_RESTART file and recreates the colvar_mtd and parvar_mtd
!files, which are necessary for the reconstruction of the free energy surface.
!Useful,if the *_mtd files got corrupted, or check the biasing hills used in 
!CPMD for actual calculations.
!Output files are in colvar_mtd.reco and parvar_mtd.reco 
!
!NN 2006-10-23 
!written by nisanth.nair@theochem.rub.de

implicit none
real*8, allocatable :: CV_PATH(:,:),CSCL_VAL(:,:),CV_DISP(:),&
                     CV_VEL(:),HLLW_VAL(:),HLLH_VAL(:)
integer :: i, j, ncv, n_mtd
character :: errinfo*10,chnum*10,lineform1*100

open(1,file='MTD_RESTART',form='unformatted',status='old',err=999)
write(*,*)' ./MTD_RESTART file found'
errinfo='n_mtd'
read(1,err=888)n_mtd
write(*,*)' Number of MTD steps =',n_mtd
errinfo='ncv'
read(1,err=888)ncv
write(*,*)' Number of CVs found =',ncv

write(chnum,'(I5)')ncv
lineform1='(1X,I7,'//trim(chnum)//'E14.6,'//trim(chnum)//'E13.3)'
allocate(CV_PATH(n_mtd,ncv),CSCL_VAL(n_mtd,ncv),CV_DISP(ncv),CV_VEL(ncv),&
         HLLW_VAL(n_mtd),HLLH_VAL(n_mtd))

errinfo='cv_path'
do i=1,ncv
 read(1,err=888)(CV_PATH(j,i),j=1,n_mtd)
 read(1,err=888)(CSCL_VAL(j,i),j=1,n_mtd)
end do

errinfo='cv_disp'
read(1,err=888)(CV_DISP(i),i=1,ncv)
errinfo='cv_vel'
read(1,err=888)(CV_VEL(i),i=1,ncv)

errinfo='hllw_val'
read(1,err=888) (HLLW_VAL(j),j=1,n_mtd)
errinfo='hllh_val'
READ(1,err=888) (HLLH_VAL(j),j=1,n_mtd)

close(1)

print *, " writing colvar_mtd.reco, parvar_mtd.reco" 
open(2,file='colvar_mtd.reco',status='unknown',form='formatted')
open(3,file='parvar_mtd.reco',status='unknown',form='formatted')
do j=1,n_mtd
  write(2,lineform1)j, (CV_PATH(j,i),i=1,ncv), (CSCL_VAL(j,i),i=1,ncv)
  write(3,'(1X,I7,3F14.6)')j, 0.0d0, HLLW_VAL(j), HLLH_VAL(j)
end do
close(2)
close(3)
deallocate(CV_PATH,CSCL_VAL,CV_DISP,CV_VEL,HLLW_VAL,HLLH_VAL)
stop

888 continue
    print *, errinfo
    stop 'ERROR: while reading'
999 continue
    print *, errinfo
    stop 'ERROR: ./MTD_RESTART file not found'
end program
