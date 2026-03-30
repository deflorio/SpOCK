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

!**************************************************************************************************
! MSIS_UTILS Module: Contains the following auxiliary subroutines:
!  alt2gph:  Converts geodetic altitude to geopotential height
!  gph2alt:  Converts geopotential height to geodetic altitude
!  bspline:  Computes B-splines using input nodes and up to specified order
!  dilog:    Computes dilogarithm function (expansion truncated at order 3, error < 1E-5)
!**************************************************************************************************

module msis_utils

contains

  !==================================================================================================
  ! ALT2GPH: Altitude to Geopotential Height
  ! References:
  !   DMA Technical Report TR8350.2 (1987),
  !     http://earth-info.nga.mil/GandG/publications/historic/historic.html
  !   Featherstone, W. E., and S. J. Claessens (2008), Closed-form transformation between
  !     geodetic and ellipsoidal coordinates, Studia Geophysica et Geodaetica, 52, 1-18
  !   Jekeli, C. (2009), Potential theory and static gravity field of the Earth, in 
  !     Treatise on Geophysics, ed. T. Herring, vol 3, 11-42
  !   NIMA Technical Report TR8350.2 (2000, 3rd edition, Amendment1), 
  !     http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350_2.html
  !==================================================================================================
  real(8) function alt2gph(lat,alt)

    implicit none

    ! Input variables
    real(8), intent(in) :: lat    !Geodetic latitude (deg)
    real(8), intent(in) :: alt    !Geodetic altitude (km)

    real(8), parameter  :: deg2rad = 0.017453292519943295d0

    ! WGS84 Defining parameters
    real(8), parameter  :: a = 6378.1370d0 * 1d3 !Semi-major axis of reference ellipsoid (m)
    real(8), parameter  :: finv = 298.257223563d0 ! 1/f = Reciprocal of flattening
    real(8), parameter  :: w = 7292115d-11 !Angular velocity of Earth rotation (rad/s)
    real(8), parameter  :: GM = 398600.4418 * 1d9 !Gravitational constant x Earth mass (m^3/s^2)

    ! WGS84 Derived parameters
    real(8), parameter  :: asq = a*a
    real(8), parameter  :: wsq = w*w
    real(8), parameter  :: f = 1.0d0 / finv
    real(8), parameter  :: esq = 2*f - f*f
    real(8), parameter  :: e = sqrt(esq)  !Ellipsoid eccentricity
    real(8), parameter  :: Elin = a*e     !Linear eccentricity of ellipsoid
    real(8), parameter  :: Elinsq = Elin*Elin
    real(8), parameter  :: epr = e / (1-f)  !Second eccentricity
    real(8), parameter  :: q0 = ((1.0d0 + 3.0d0/(epr*epr))*atan(epr) - 3.0d0/epr)/2.0d0  !DMA Technical Report tr8350.2, Eq. 3-25
    real(8), parameter  :: U0 = -GM*atan(epr)/Elin - wsq*asq/3d0 !Theoretical potential of reference ellipsoid (m^2/s^2)
    real(8), parameter  :: g0 = 9.80665d0 !Standard gravity (m/s^2), CGPM 1901; WMO
    real(8), parameter  :: GMdivElin = GM / Elin
  
    ! Parameters for centrifugal potential taper
    real(8), parameter  :: x0sq = 2d7**2   !Axial distance squared at which tapering begins (m^2)
    real(8), parameter  :: Hsq = 1.2d7**2  !Relaxation scale length of taper (m^2)

    ! Working variables
    real(8)             :: altm, sinsqlat, v, xsq, zsq
    real(8)             :: rsqminElinsq, usq, cossqdelta, epru, atanepru, q, U, Vc

    ! Compute Cartesian and ellipsoidal coordinates
    altm = alt * 1000.0d0
    sinsqlat = sin(lat*deg2rad)**2
    v = a / sqrt(1-esq*sinsqlat)           !Radius of curvature of the reference ellipsoid, Featherstone eq. 4
    xsq = (v + altm)**2 * (1 - sinsqlat)   !Squared x-coordinate of geocentric system, Featherstone eq. 1
    zsq = (v*(1-esq) + altm)**2 * sinsqlat !Squared z-coordinate of geocentric system, Featherstone eq. 3
    rsqminElinsq = xsq + zsq - Elinsq
    usq = rsqminElinsq/2.0d0 + sqrt(rsqminElinsq**2 / 4.0d0 + Elinsq*zsq)  !Ellipsoidal distance coordinate, Featherstone eq. 19 
    cossqdelta = zsq / usq                 !Ellipsoidal polar angle, Featherstone eq. 21

    ! Compute gravitational potential
    epru = Elin / sqrt(usq)                !Second eccentricity at ellipsoidal coordinate u
    atanepru = atan(epru)
    q = ((1+3.0d0/(epru*epru))*atanepru - 3.0d0/epru)/2.0d0   !Jekeli, eq. 114
    U = -GMdivElin * atanepru - wsq * ( asq * q * (cossqdelta - 1/3.0d0) / q0 ) / 2.0d0   !Jekeli, eq. 113

    ! Compute centrifugal potential and adjust total potential
    if (xsq .le. x0sq) then
      Vc = (wsq/2.0d0) * xsq
    else
      Vc = (wsq/2.0d0) * (Hsq*tanh((xsq-x0sq)/Hsq) + x0sq) !Centrifugal potential taper
    endif
    U = U - Vc
  
    ! Compute geopotential height
    alt2gph = (U - U0) / g0 / 1000.0d0

    return

  end function alt2gph

  !==================================================================================================
  ! GPH2ALT: Geopotential Height to Altitude
  !==================================================================================================
  real(8) function gph2alt(theta,gph)

    implicit none

    real(8), intent(in)  :: theta
    real(8), intent(in)  :: gph

    integer, parameter   :: maxn = 10
    real(8), parameter   :: epsilon = 0.0005

    real(8)              :: x,dx,y,dydz
    integer              :: n

    x = gph
    n = 0
    dx = epsilon + epsilon
    do while ((abs(dx) .gt. epsilon) .and. (n .lt. 10))
      y = alt2gph(theta,x)
      dydz = (alt2gph(theta,x+dx) - y)/dx
      dx = (gph - y)/dydz
      x = x + dx
      n = n + 1
    end do

    gph2alt = x

  end function gph2alt

  !==================================================================================================
  ! BSPLINE: Returns array of nonzero b-spline values, for all orders up to specified order (max 6)
  !==================================================================================================
  subroutine bspline(x,nodes,nd,kmax,eta,S,i)

    use msis_constants, only:  rp

    implicit none

    ! Input variables
    real(kind=rp), intent(in)         :: x              !Location at which splines are to be evaluated
    real(kind=rp),dimension(0:),intent(in)  :: nodes    !Spline node locations
    integer, intent(in)               :: nd             !Number of spline nodes minus one (0:nd)
    integer, intent(in)               :: kmax           !Maximum order (up to 6 allowed) of evaluated splines
    real(kind=rp), intent(in)         :: eta(0:30,2:6)  !Precomputed weights for recursion (reciprocals of node differences)
    ! Ouput variables
    real(kind=rp), intent(out)        :: S(-5:0,2:6)    !b-spline values (spline index relative to i (-5:0), spline order (2:6))
    integer, intent(out)              :: i              !Index of last nonzero b-spline

    ! Working variables
    integer                           :: j, k, l
    integer                           :: low, high
    real(kind=rp)                     :: w(-4:0)        !Weights for recursion relation
 
    ! Initialize to zero
    S(:,:) = 0.0_rp

    ! Find index of last (rightmost) nonzero spline
    if (x .ge. nodes(nd)) then
      i = nd
      return
    endif
    if (x .le. nodes(0)) then
      i = -1
      return
    endif
    low = 0
    high = nd
    i = (low + high)/2
    do while (x .lt. nodes(i) .or. x .ge. nodes(i + 1))
        if (x .lt. nodes(i)) then
          high = i
        else
          low = i
        endif
        i = (low + high)/2
    end do

    ! Initialize with linear splines
    S(0,2) = (x - nodes(i)) * eta(i,2)
    if (i .gt. 0) S(-1,2) = 1 - S(0,2)
    if (i .ge. nd-1) S(0,2) = 0.0_rp   !Reset out-of-bounds spline to zero

    ! k = 3 (quadratic splines)
    w(:) = 0.0_rp
    w(0) = (x - nodes(i)) * eta(i,3)
    if (i .ne. 0) w(-1) = (x - nodes(i-1)) * eta(i-1,3)
    if (i .lt. (nd-2)) S(0,3) = w(0)*S(0,2)
    if ( ((i-1) .ge. 0) .and. ((i-1) .lt. (nd-2)) ) &
        S(-1,3) = w(-1) * S(-1,2) + (1.0_rp - w(0))*S(0,2)
    if ((i-2) .ge. 0) S(-2,3) = (1.0_rp - w(-1))*S(-1,2)
    
    ! k = 4 (cubic splines)
    do l = 0, -2, -1
      j = i + l
      if (j .lt. 0) exit  !Skip out-of-bounds splines
      w(l) = (x - nodes(j)) * eta(j,4)
    enddo
    if (i .lt. (nd-3)) S(0,4) = w(0)*S(0,3)
    do l = -1, -2, -1
        if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-3)) ) &
            S(l,4) = w(l)*S(l,3) + (1.0_rp - w(l+1))*S(l+1,3)
    enddo
    if ((i-3) .ge. 0) S(-3,4) = (1.0_rp - w(-2))*S(-2,3)
  
    ! k = 5
    do l = 0, -3, -1
      j = i + l
      if (j .lt. 0) exit  !Skip out-of-bounds splines
      w(l) = (x - nodes(j)) * eta(j,5)
    enddo
    if (i .lt. (nd-4)) S(0,5) = w(0)*S(0,4)
    do l = -1, -3, -1
        if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-4)) ) &
            S(l,5) = w(l)*S(l,4) + (1.0_rp - w(l+1))*S(l+1,4)
    enddo
    if ((i-4) .ge. 0) S(-4,5) = (1.0_rp - w(-3))*S(-3,4)
    if (kmax .eq. 5) return  !Exit if only 5th order spline is needed

    ! k = 6
    do l = 0, -4, -1
      j = i + l
      if (j .lt. 0) exit  !Skip out-of-bounds splines
      w(l) = (x - nodes(j)) * eta(j,6)
    enddo
    if (i .lt. (nd-5)) S(0,6) = w(0)*S(0,5)
    do l = -1, -4, -1
      if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-5)) ) &
          S(l,6) = w(l)*S(l,5) + (1.0_rp - w(l+1))*S(l+1,5)
    enddo
    if ((i-5) .ge. 0) S(-5,6) = (1.0_rp - w(-4))*S(-4,5)

    return

  end subroutine bspline

  !==================================================================================================
  ! DILOG: Calculate dilogarithm in the domain [0,1)
  ! Retains terms up to order 3 in the expansion, which results in relative errors less than 1E-5.
  ! Reference: 
  !   Ginsberg, E. S., and D. Zaborowski (1975), The Dilogarithm function of a real argument, 
  !   Commun. ACM, 18, 200–202.
  !==================================================================================================
  real(kind=rp) function dilog(x0)

    use msis_constants, only     : rp, pi

    implicit none

    real(kind=rp), intent(in)   :: x0
    real(kind=rp), parameter    :: pi2_6 = pi*pi / 6.0_rp
    real(kind=rp)               :: x, xx, x4, lnx

    x = x0
    if (x .gt. 0.5_rp) then
      lnx = log(x)
      x = 1.0_rp - x          !Reflect argument into [0,0.5] range
      xx = x*x
      x4 = 4.0_rp*x
      dilog = pi2_6 - lnx*log(x) &
              - (4.0_rp*xx*(23.0_rp/16.0_rp + x/36.0_rp + xx/576.0_rp + xx*x/3600.0_rp) &
                  + x4 + 3.0_rp*(1.0_rp - xx)*lnx) / (1.0_rp + x4 + xx)
    else
      xx = x*x
      x4 = 4.0_rp*x
      dilog = (4.0_rp*xx*(23.0_rp/16.0_rp + x/36.0_rp + xx/576.0_rp + xx*x/3600.0_rp) &
                + x4 + 3.0_rp*(1.0_rp - xx)*log(1.0_rp - x)) / (1.0_rp + x4 + xx)
    endif

    return

  end function dilog

end module msis_utils