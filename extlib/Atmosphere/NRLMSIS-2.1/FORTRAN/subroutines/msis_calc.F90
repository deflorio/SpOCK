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
!!!
!!! MSISCALC: Interface with re-ordered input arguments and output arrays.
!
!     PREREQUISITES:
!       Must first run MSISINIT to load parameters and set switches. The 
!       MSISCALC subroutine checks for initialization and does a default
!       initialization if necessary. This self-initialization will be removed
!       in future versions.
!
!     CALLING SEQUENCE:
!       CALL MSISCALC(DAY, UTSEC, Z, LAT, LON, SFLUXAVG, SFLUX, AP, TN, DN, [TEX])
!  
!     INPUT VARIABLES:
!       DAY       Day of year (1.0 to 365.0 or 366.0)
!       UTSEC     Universal time (seconds)
!       Z         Geodetic altitude (km) (default) or Geopotential height (km)
!       LAT       Geodetic latitude (deg)
!       LON       Geodetic longitude (deg)
!       SFLUXAVG  81 day average, centered on input time, of F10.7 solar
!                 activity index
!       SFLUX     Daily F10.7 for previous day
!       AP        Geomagnetic activity index array:
!                   (1) Daily Ap
!                   (2) 3 hr ap index for current time
!                   (3) 3 hr ap index for 3 hrs before current time
!                   (4) 3 hr ap index for 6 hrs before current time
!                   (5) 3 hr ap index for 9 hrs before current time
!                   (6) Average of eight 3 hr ap indices from 12 to 33 hrs
!                       prior to current time
!                   (7) Average of eight 3 hr ap indices from 36 to 57 hrs
!                       prior to current time
!                 AP(2:7) are only used when switch_legacy(9) = -1.0 in MSISINIT
!
!     NOTES ON INPUT VARIABLES: 
!       - The day-of-year dependence of the model only uses the DAY argument. If
!         a continuous day-of-year dependence is desired, this argument should
!         include the fractional day (e.g., DAY = <day of year> + UTSEC/86400.0
!       - If lzalt_type = .true. (default) in the MSISINIT call, then Z is
!         treated as geodetic altitude.
!         If lzalt_type = .false., then Z is treated as geopotential height.
!       - F107 and F107A values are the 10.7 cm radio flux at the Sun-Earth
!         distance, not the radio flux at 1 AU. 
!
!     OUTPUT VARIABLES:
!       TN     Temperature at altitude (K)
!       DN(1)  Total mass density (kg/m3)
!       DN(2)  N2 number density (m-3)
!       DN(3)  O2 number density (m-3)
!       DN(4)  O number density (m-3)
!       DN(5)  He number density (m-3)
!       DN(6)  H number density (m-3)
!       DN(7)  Ar number density (m-3)
!       DN(8)  N number density (m-3)
!       DN(9)  Anomalous oxygen number density (m-3)
!       DN(10) NO number density (m-3)
!       TEX    Exospheric temperature (K) (optional argument)
!
!     NOTES ON OUTPUT VARIABLES: 
!       - Missing density values are returned as 9.999e-38
!       - Species included in mass density calculation are set in MSISINIT
!
!!! =========================================================================

!**************************************************************************************************
! MSIS_CALC Module: Contains main MSIS entry point
!**************************************************************************************************
module msis_calc

contains

  !==================================================================================================
  ! MSISCALC: The main MSIS subroutine entry point
  !==================================================================================================
  subroutine msiscalc(day,utsec,z,lat,lon,sfluxavg,sflux,ap,tn,dn,tex)

    use msis_constants, only    : rp, dmissing, lnp0, Mbarg0divkB, kB, nspec, nodesTN, nd, zetaF, zetaB, &
                                  Hgamma, zetagamma, maxnbf
    use msis_init, only         : msisinit, initflag, zaltflag, specflag, massflag, masswgt, etaTN
    use msis_gfn, only          : globe
    use msis_tfn, only          : tnparm, tfnparm, tfnx
    use msis_dfn, only          : dnparm, dfnparm, dfnx
    use msis_utils, only        : alt2gph, bspline, dilog

    implicit none

    real(kind=rp), intent(in)  :: day
    real(kind=rp), intent(in)  :: utsec
    real(kind=rp), intent(in)  :: z
    real(kind=rp), intent(in)  :: lat
    real(kind=rp), intent(in)  :: lon
    real(kind=rp), intent(in)  :: sfluxavg,sflux,ap(1:7)
    real(kind=rp), intent(out) :: tn, dn(1:10)
    real(kind=rp), intent(out), optional :: tex
  
    real(kind=rp), save        :: lastday = -9999.0
    real(kind=rp), save        :: lastutsec = -9999.0
    real(kind=rp), save        :: lastlat = -9999.0
    real(kind=rp), save        :: lastlon = -9999.0
    real(kind=rp), save        :: lastz = -9999.0
    real(kind=rp), save        :: lastsflux = -9999.0
    real(kind=rp), save        :: lastsfluxavg = -9999.0
    real(kind=rp), save        :: lastap(1:7) = -9999.0
    real(kind=rp), save        :: gf(0:maxnbf-1)
    real(kind=rp), save        :: Sz(-5:0,2:6)
    integer, save              :: iz
    type(tnparm), save         :: tpro
    type(dnparm), save         :: dpro(1:nspec-1)

    real(8)                    :: zaltd, latd
    real(kind=rp)              :: zeta, lndtotz, Vz, Wz, HRfact, lnPz, delz
    integer                    :: i, j, kmax, ispec

    ! Check if model has been initialized; if not, perform default initialization
    if (.not. initflag) call msisinit()

    ! Calculate geopotential height, if necessary
    if(zaltflag) then
      zaltd = dble(z)
      latd = dble(lat)
      zeta = alt2gph(latd,zaltd)
    else
      zeta = z
    endif

    ! If only altitude changes then update the local spline weights
    if (zeta .lt. zetaB) then
      if (zeta .ne. lastz) then
        if (zeta .lt. zetaF) then
          kmax = 5
        else
          kmax = 6
        endif
        call bspline(zeta,nodesTN,nd+2,kmax,etaTN,Sz,iz)
        lastz = zeta
      endif
    endif

    ! If location, time, or solar/geomagnetic conditions change then recompute the profile parameters
    if ((day .ne. lastday)     .or. (utsec .ne. lastutsec)       .or. &
        (lat .ne. lastlat)     .or. (lon .ne. lastlon)           .or. &
        (sflux .ne. lastsflux) .or. (sfluxavg .ne. lastsfluxavg) .or. &
        any(ap .ne. lastap)) then
      call globe(day,utsec,lat,lon,sfluxavg,sflux,ap,gf)
      call tfnparm(gf,tpro)
      do ispec = 2, nspec-1
        if (specflag(ispec)) call dfnparm(ispec,gf,tpro,dpro(ispec))
      enddo
      lastday = day
      lastutsec = utsec
      lastlat = lat
      lastlon = lon
      lastsflux = sflux
      lastsfluxavg = sfluxavg
      lastap = ap
    endif

    ! Exospheric temperature
    if (present(tex)) then
      tex = tpro%tex
    endif

    ! Temperature at altitude
    tn = tfnx(zeta,iz,Sz(-3:0,4),tpro)

    ! Temperature integration terms at altitude, total number density
    delz = zeta - zetaB
    if (zeta .lt. zetaF) then
      i = max(iz-4,0)
      if (iz .lt. 4) then
        j = -iz
      else
        j = -4
      endif
      Vz = dot_product(tpro%beta(i:iz),Sz(j:0,5)) + tpro%cVS
      Wz = 0.0_rp
      lnPz = lnP0 - Mbarg0divkB*(Vz - tpro%Vzeta0)
      lndtotz = lnPz - log(kB*tn)
    else
      if (zeta .lt. zetaB) then
        Vz = dot_product(tpro%beta(iz-4:iz),Sz(-4:0,5)) + tpro%cVS
        Wz = dot_product(tpro%gamma(iz-5:iz),Sz(-5:0,6)) + tpro%cVS*delz + tpro%cWS
      else
        Vz = (delz + log(tn/tpro%tex)/tpro%sigma)/tpro%tex + tpro%cVB
        Wz = (0.5_rp*delz*delz + dilog(tpro%b*exp(-tpro%sigma*delz))/tpro%sigmasq)/tpro%tex &
              + tpro%cVB*delz + tpro%cWB
      endif
    endif
        
    ! Species number densities at altitude
    HRfact = 0.5_rp * (1.0_rp + tanh(Hgamma*(zeta - zetagamma)))  !Reduction factor for chemical/dynamical correction scale height below zetagamma
    do ispec = 2, nspec-1
      if (specflag(ispec)) then
        dn(ispec) = dfnx(zeta,tn,lndtotz,Vz,Wz,HRfact,tpro,dpro(ispec))
      else
        dn(ispec) = dmissing
      endif
    enddo

    ! Mass density
    if (specflag(1)) then
      dn(1) = dot_product(dn,masswgt)
    else
      dn(1) = dmissing
    endif

    return

  end subroutine msiscalc

end module msis_calc
