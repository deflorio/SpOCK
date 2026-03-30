!#######################################################################
! MSIS® (NRL-SOF-014-1) SOFTWARE
! NRLMSIS® empirical atmospheric model software. Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! nrlmsis2.1_license.txt, which is part of this software package. BY
! USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
! CONDITIONS OF THE LICENSE.  
!#######################################################################

!!! ===========================================================================
!!! NRLMSIS 2.1:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! ===========================================================================

!==================================================================================================
! MSISTEST: Test program for NRLMSIS 2.1
!==================================================================================================
program msistest

  use msis_init, only          : msisinit

  implicit none

  integer, parameter          :: nrec = 200

  integer                     :: iyd, mass
  real(4)                     :: sec, alt, glat, glong, stl, f107a, f107, ap(7), apd
  real(4)                     :: d(10),t(2)
  
  integer                     :: i
  character(128)              :: dummy

  !Initialize model
  call msisinit(parmpath='',parmfile='msis21.parm')

  !Open input and output files, loop through records, and call model
  open(77,file='msis2.1_test_in.txt',status='old')
  open(78,file='msis2.1_test_out.txt',status='replace')
  read(77,*) dummy
  write(78,'(9a7,10a13,a8)') &
    'iyd','sec','alt','glat','glong','stl','f107a','f107','Ap','He','O','N2','O2','Ar','rho','H','N','O*','NO','T'
  do i = 1,200
	  read(77,*) iyd,sec,alt,glat,glong,stl,f107a,f107,apd
    ap(1) = apd
    call gtd8d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
    write(78,'(2i7,3f7.1,f7.2,3f7.1,10e13.4,f8.2)')  &
      iyd,int(sec),alt,glat,glong,stl,f107a,f107,ap(1),d(1:10),t(2)
  enddo
  close(77)
  close(78)

  stop

end program msistest